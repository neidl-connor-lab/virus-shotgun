#!/bin/bash -l

# qsub options
#$ -l h_rt=24:00:00
#$ -pe omp 4
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup -----------------------------------------------------------------------
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
FORMAT="pipeline/format.r"
MINCOV="10"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -i INDEX -f FASTA -o ODIR -s SAMPLE -x R1 [-y R2] [-t THLD]
Please submit the job from the pipeline directory!

arguments:
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -t [OPTIONAL] minimum aligned read depth (default: $MINCOV)
  -h print this message and exit
"

# parsing arguments
while getopts ":hi:f:o:s:x:y:t:" opt 
do 
  case ${opt} in 
    i ) IDX="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    s ) SAMPLE="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    t ) MINCOV="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## print job info for output log -----------------------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# if no kraken, fail
if [ -z "$(which kraken2 2> /dev/null)" ]
then
  err "No Kraken2 found. Please install and try again."
fi

# double-check that lofreq exists
if [ -z "$($LOFREQ version 2> /dev/null)" ]
then
  err "LoFreq error: $LOFREQ"
fi

# double-check the formatting script exists
if [ ! -f "$FORMAT" ]
then 
  err "VCF formatting script not detected: $FORMAT"
fi

# load bowtie2 module
module load bowtie2/2.4.2
checkcmd "Loading bowtie2"

# load samtools module 
module load htslib/1.18
module load samtools/1.18
checkcmd "Loading samtools"

# load R
module load R/4.0.2
checkcmd "Loading R"

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

# kraken2 db should have 12 files
if [ "$(ls -1 $KDB 2> /dev/null | wc -l)" -ne 12 ]
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

# check sample-specific output directory
VAR="$ODIR/$SAMPLE"
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
elif [ ! -d "$VAR" ]
then
  mesg "Outputting files to: $VAR"
  mkdir -p "$VAR"
else
  mesg "Outputting files to: $VAR"
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

# minimum coverage
if [ -z "$MINCOV" ]
then
  err "No minimum coverage specified"
else
  mesg "Minimum aligned read depth: $MINCOV"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## contaminant check -----------------------------------------------------------
mesg "STEP 1: CONTAMINANT CHECK"

# build command based on whether has paired reads
CMD="kraken2 --db '$KDB' --threads 4 --output - --report '$VAR/metagenomics.tsv' --use-names"
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
mesg "STEP 2: ALIGN TO VIRAL GENOME"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 4 -x '$IDX' -U '$R1' 1> '$VAR/alignment.sam' 2> '$VAR/bowtie2.log'"
else # paired
  CMD="bowtie2 --threads 4 -x '$IDX' -1 '$R1' -2 '$R2' 1> '$VAR/alignment.sam' 2> '$VAR/bowtie2.log'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Alignment"

# compress SAM to BAM
CMD="samtools view --threads 4 -b -h '$VAR/alignment.sam' > '$VAR/alignment-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$VAR/alignment.sam"
echo ""

## process BAM --------------------------------------------------
mesg "STEP 3: PROCESS BAM"

# sort BAM
CMD="samtools sort --threads 4 '$VAR/alignment-raw.bam' > '$VAR/alignment-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$VAR/alignment-raw.bam"

# score indels to get final BAM
CMD="$LOFREQ indelqual --dindel --ref '$REFSEQ' '$VAR/alignment-sorted.bam' > '$VAR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$VAR/alignment-sorted.bam"

# index final BAM
CMD="samtools index '$VAR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage -------------------------------------------
mesg "STEP 4: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="samtools depth --threads 4 -a -H '$VAR/alignment.bam' > '$VAR/coverage.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus -------------------------------------------
mesg "STEP 5: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="samtools consensus --threads 4 --use-qual --min-depth $MINCOV --call-fract 0.5 --output '$VAR/consensus-tmp.fa' '$VAR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"

# update consensus header
mesg "Updating consensus header"
echo ">$SAMPLE" > "$VAR/consensus.fa"
cat "$VAR/consensus-tmp.fa" | grep "^[^>]" >> "$VAR/consensus.fa"
rm "$VAR/consensus-tmp.fa"
echo ""

## quantify SNVs ------------------------------------------------
mesg "STEP 6: QUANTIFY SNVs"

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="$LOFREQ call-parallel --pp-threads 4 --call-indels --min-cov $MINCOV --ref '$REFSEQ' '$VAR/alignment.bam' > '$VAR/snvs.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"
echo ""

# format VCF
CMD="Rscript $FORMAT --vcf '$VAR/snvs.vcf' --ofile '$VAR/snvs.csv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "SNV formatting"
echo ""

## package version ----------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
kraken2 --version
echo ""
echo "LoFreq"
$LOFREQ version
echo ""
