#!/bin/bash -l

# qsub options
#$ -l h_rt=48:00:00
#$ -l mem_per_core=12G
#$ -pe omp 16
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup --------------------------------------------------------
# functions
mesg () { echo -e "[MSG] $@"; }
err () { echo -e "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}

# default values and help message
KDB="$(dirname "$0")/pipeline/kraken2db"
VODKA="$(dirname "$0")/pipeline/vodka/VODKA"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") -o OUTPUT -x R1 -y R2 -b BOWTIE2 -f FASTA [-k KRAKEN2 -v VODKA]

arguments (options):
  -o output directory
  -x R1 FASTQ file
  -y R2 FASTQ file
  -b bowtie2 index
  -f virus reference FASTA
  -k kraken2 database (default: $KDB)
  -v VODKA installation (default: $VODKA)
  -h show this message and exit
"

# parsing arguments
while getopts ":ho:x:y:b:f:k:v:" opt 
do 
  case ${opt} in 
    o ) ODIR="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    b ) BOWTIE2="${OPTARG}"
      ;;
    f ) FASTA="${OPTARG}"
      ;;
    k ) KDB="${OPTARG}"
      ;;
    v ) VODKA="${OPTARG}"
      ;;
    h ) echo "${HELP}" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## print job info for output log --------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID"
echo "=========================================================="
echo ""

## check inputs -------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# check kraken2 is in $PATH
if [ -z "$(which kraken2 2> /dev/null)"]
then
  err "kraken2 not found"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided."
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# R1 FASTQ file
if [ -z "$R1" ]
then
  err "No R1 FASTQ file provided"
elif [ -f "$R1" ]
then
  mesg "Valid R1 FASTQ file: $R1"
else
  err "Invalid R1 FASTQ file: $R1"
fi

# R2 FASTQ file
if [ -z "$R2" ]
then
  err "No R2 FASTQ file provided"
elif [ -f "$R2" ]
then
  mesg "Valid R2 FASTQ file: $R2"
else
  err "Invalid R2 FASTQ file: $R2"
fi

# bowtie2 index
if [ -z "$BOWTIE2" ]
then
  err "No bowtie2 index provided"
elif [ "$(ls -1 ${BOWTIE2}*.bt2 | wc -l)" -eq "6" ]
then
  mesg "Valid bowtie2 index: $BOWTIE2"
else
  err "Invalid bowtie2 index: $BOWTIE2"
fi

# fasta reference sequence
if [ -z "$FASTA" ]
then
  err "No reference FASTA provided"
elif [ -f "$FASTA" ]
then
  mesg "Valid reference FASTA: $FASTA"
else
  err "Invalid reference FASTA: $FASTA"
fi

# kraken2 database
if [ -d "$KDB" ]
then
  mesg "Kraken2 database directory: $KDB"
else
  err "Invalid kraken2 database directory: $KDB"
fi

# vodka executable
if [ -f "$VODKA" ]
then
  mesg "VODKA executable: $VODKA"
else
  err "Invalid VODKA executable: $VODKA"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## kraken2 classification ---------------------------------------
mesg "STEP 1: KRAKEN2 CONTAMINANT CHECK"

# make sure gcc has been loaded
module load gcc

# build and execute command
CMD="kraken2 --threads 16 --paired --db '$KDB' --output '$ODIR/kraken2-output.tsv' --report '$ODIR/kraken2-report.tsv' --use-names --gzip-compressed '$R1' '$R2'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "kraken2"
echo ""

## bowtie2 alignment --------------------------------------------
mesg "STEP 2: BOWTIE2 ALIGNMENT"

# load module
module load bowtie2
module load samtools

# build and execute alignment
mesg "Align to virus reference sequence"
CMD="bowtie2 --threads 16 -x '$BOWTIE2' -1 '$R1' '$R2' > '$ODIR/alignment.sam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "bowtie2 alignment"

# compress and sort
mesg "Compress and sort alignment"
samtools view --threads 16 -b -h "$ODIR/alignment.sam" > "$ODIR/alignment.bam"
samtools sort --threads 16 --reference "$FASTA" -o "$ODIR/alignment-sorted.bam" --output-fmt BAM "$ODIR/alignment.bam"

# clean up intermediate files
rm "$ODIR/alignment.sam"
rm "$ODIR/alignment.bam"
mv "$ODIR/alignment-sorted.bam" "$ODIR/alignment.bam"
echo ""

## coverage -----------------------------------------------------
mesg "STEP 3: CALCULATE COVERAGE"

# build command and run
CMD="samtools depth -a -H -d 0 -o '$ODIR/coverage.tsv' '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "coverage"
echo ""

## lofreq and SNV calling ---------------------------------------


## print package versions ---------------------------------------
mesg "FIN. PIPELINE COMPLETED SUCCESSFULLY."
kraken2 --version
module list
echo ""

