#!/bin/bash -l

# qsub options
#$ -l h_rt=24:00:00
#$ -pe omp 8
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

# pre-set variables
LOFREQ="pipeline/lofreq/lofreq"
KDB="pipeline/kraken2db"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -i INDEX -f FASTA -o ODIR -x R1 [-y R2]
Please submit the job from the pipeline directory!

arguments:
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory (e.g., sample ID)
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -h print this message and exit
"

# parsing arguments
while getopts ":hi:f:o:s:x:y:" opt 
do 
  case ${opt} in 
    i ) IDX="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
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
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs -------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# if no kraken, fail
if [ -z "$(which kraken2 2> /dev/null)" ]
then
  err "No Kraken2 found. Please install and try again."
fi

# if no bowtie2, load module
if [ -z "$(bowtie2 --version 2> /dev/null)" ]
then
  module load bowtie2
  checkcmd "Loading bowtie2"
fi

# if no samtools, load module
if [ -z "$(samtools version 2> /dev/null)" ]
then
  module load samtools
  checkcmd "Loading samtools"
fi

# double-check that lofreq exists
if [ -z "$($LOFREQ version 2> /dev/null)" ]
then
  err "LoFreq error: $LOFREQ"
fi

# bowtie2 index should have 6 files
if [ -z "$IDX" ]
then
  err "No bowtie2 index provided"
elif [ "$(ls -1 ${IDX}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Bowtie2 index: $IDX"
else
  err "Invalid bowtie2 index; run setup.sh first: $IDX"
fi

# kraken2 db shoudl have 4 files
if [ "$(ls -1 $KDB 2> /dev/null | wc -l)" -ne 11 ]
then
  err "Invalid Kraken2 database; rerun setup.sh"
fi

# reference sequence
if [ -z "$REFSEQ" ]
then
  err "No reference FASTA provided"
elif [ -f "$REFSEQ" ]
then
  mesg "Reference FASTA: $REFSEQ"
else
  err "Invalid reference FASTA: $REFSEQ"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided"
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# R1/R0 FASTQ file
if [ -z "$R1" ]
then
  err "No FASTQ file provided (-x)"
elif [ -f "$R1" ]
then
  mesg "First FASTQ file: $R1"
else
  err "Invalid first FASTQ file: $R1"
fi

# check for R2 FASTQ file
if [ -z "$R2" ]
then 
  mesg "Only 1 FASTQ file detected. Running as unpaired."
elif [ -f "$R2" ]
then
  mesg "Second FASTQ file: $R2"
else
  err "Invalid second FASTQ file: $R2"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## contaminant check --------------------------------------------
mesg "STEP 1: CONTAMINANT CHECK"

# build command based on whether has paired reads
CMD="kraken2 --db '$KDB' --threads 8 --output - --report '$ODIR/metagenomics.tsv' --use-names --gzip-compressed"
if [ -z "$R2" ] # unpaired
then
  CMD="$CMD '$R1'"
else # paired
  CMD="$CMD --paired '$R1' '$R2'"
fi

# run classification
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Kraken2 metagenomic classification"
echo ""

## alignment ----------------------------------------------------
mesg "STEP 2: ALIGN TO GENOME"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 -x '$IDX' -U '$R1' > '$ODIR/alignment.sam'"
else # paired
  CMD="bowtie2 --threads 8 -x '$IDX' -1 '$R1' -2 '$R2' > '$ODIR/alignment.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Alignment"

# compress SAM to BAM
CMD="samtools view --threads 8 -b -h '$ODIR/alignment.sam' > '$ODIR/alignment-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$ODIR/alignment.sam"
echo ""

## process BAM --------------------------------------------------
mesg "STEP 3: PROCESS BAM"

# sort BAM
CMD="samtools sort --threads 8 '$ODIR/alignment-raw.bam' > '$ODIR/alignment-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$ODIR/alignment-raw.bam"

# score indels to get final BAM
CMD="$LOFREQ indelqual --dindel --ref '$REFSEQ' '$ODIR/alignment-sorted.bam' > '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$ODIR/alignment-sorted.bam"

# index final BAM
CMD="samtools index '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage -------------------------------------------
mesg "STEP 4: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="samtools depth --threads 8 -a -H '$ODIR/alignment.bam' > '$ODIR/coverage.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus -------------------------------------------
mesg "STEP 5: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="samtools consensus --threads 8 --use-qual --min-depth 10 --call-fract 0.5 --output '$ODIR/consensus-tmp.fa' '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"

# update consensus header
mesg "Updating consensus header"
echo ">$ID" > "$ODIR/consensus.fa"
cat "$ODIR/consensus-tmp.fa" | grep "^[^>]" >> "$ODIR/consensus.fa"
rm "$ODIR/consensus-tmp.fa"
echo ""

## quantify SNVs ------------------------------------------------
mesg "STEP 6: QUANTIFY SNVs"

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="$LOFREQ call-parallel --pp-threads 8 --call-indels --min-cov 10 --ref '$REFSEQ' '$ODIR/alignment.bam' > '$ODIR/snvs.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"
echo ""

## package version ----------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
kraken2 --version
echo ""
echo "LoFreq"
$LOFREQ version
echo ""
samtools --version
echo "" 
