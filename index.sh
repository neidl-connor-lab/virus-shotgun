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
RESOURCES="pipeline"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") -f FASTA -b BOWTIE

arguments (options):
  -f viral genome FASTA
  -b bowtie2 index directory and name
  -h show this message and exit"

# parsing arguments
while getopts ":hf:b:" opt
do
  case ${opt} in
    f ) FASTA="${OPTARG}"
      ;;
    b ) BOWTIE2="${OPTARG}"
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
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs -------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# viral genome FASTA
if [ -z "$FASTA" ]
then
  err "No reference genome FASTA provided!"
elif [ -f "$FASTA" ]
then
  mesg "Valid reference genome: $FASTA"
else
  err "Invalid reference genome: $FASTA"
fi

# check bowtie2 was provided
if [ -z "$BOWTIE2" ]
then
  err "No bowtie2 prefix provided"
fi

# bowtie2 directory
if [ -d "$(dirname $BOWTIE2)" ]
then
  mesg "Valid bowtie2 index directory: $(dirname $BOWTIE2)"
else
  mesg "Creating bowtie2 index directory: $(dirname $BOWTIE2)"
  mkdir -p "$(dirname $BOWTIE2)"
fi

# bowtie2 name
mesg "Index name: $(basename $BOWTIE2)"

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## bowtie index -------------------------------------------------

# load bowtie2 
module load bowtie2

# set up and run command
mesg "Creating bowtie2 index from reference"
CMD="bowtie2-build --threads 16 --quiet '$FASTA' '$BOWTIE2'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "bowtie2 index"
echo ""

## all done -----------------------------------------------------
mesg "Bowtie2 index complete!"
module list
echo ""
