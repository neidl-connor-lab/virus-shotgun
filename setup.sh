#!/bin/bash -l

# qsub options
#$ -l h_rt=12:00:00
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
LOFREQ="pipeline/lofreq"
KDIR="pipeline/kraken2db"
KLINK="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20221209.tar.gz"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index path and prefix
  -h show this message and exit"

# parsing arguments
while getopts ":hf:b:" opt
do
  case ${opt} in
    f ) FASTA="${OPTARG}"
      ;;
    b ) BOWTIE="${OPTARG}"
      ;;
    h ) echo "${HELP}" && exit 0
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
# load kraken2
module load kraken2/2.1.2
checkcmd "Loading kraken2"

# load bowtie2
module load bowtie2/2.4.2
checkcmd "Loading bowtie2" 

# load samtools module 
module load htslib/1.18
module load samtools/1.18
checkcmd "Loading samtools"

# check reference FASTA
if [ -z "$FASTA" ]
then
  err "No reference FASTA provided"
elif [ -f "$FASTA" ]
then
  mesg "Valid reference FASTA file: $FASTA"
else
  err "Invalid reference FASTA file: $FASTA"
fi

# check that bowtie prefix was provided
if [ -z "$BOWTIE" ]
then
  err "No bowtie2 index name provided"
else
  mesg "Bowtie2 index prefix: $BOWTIE"
fi

# make bowtie index output directory if necessary
if [ ! -d "$(dirname $BOWTIE)" ]
then
  mesg "Creating bowtie2 index directory: $(dirname $BOWTIE)"
  mkdir -p "$(dirname $BOWTIE)"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## unpack LoFreq if necessary --------------------------------------------------
# check for lofreq tarball
if [ -f "$LOFREQ.tar.bz2" ]
then 
  mesg "LoFreq tarball detected. Expanding."
  tar -xf "$LOFREQ.tar.bz2" -C "$(dirname $LOFREQ)"
  checkcmd "LoFreq decompression"
  # remove tarball
  rm "$LOFREQ.tar.bz2"
fi

# check that the lofreq executable works
if [ -z "$($LOFREQ/lofreq version 2> /dev/null)" ]
then
  err "Problem with LoFreq: $LOFREQ"
else
  mesg "LoFreq is ready to go!"
fi
echo ""

## set up Kraken2 if necessary -------------------------------------------------
# check for kraken dir; should contain 4 files
if [ -d "$KDIR" ] && [ "$(ls -1 $KDIR 2> /dev/null | wc -l)" -eq 11 ]
then
  mesg "Kraken2 database is ready to go!"
# if any problem with the folder or it doesn't exists, make a new db
else
  mesg "Setting up the Kraken2 database. This will take a while!"
  # remove the directory if it exists
  if [ -d "$KDIR" ]
  then
    rm -r "$KDIR"
  fi
  
  # move into the pipeline directory and check for tarball
  # if it's there, remove it and pull a fresh one
  cd "$(dirname $KDIR)"
  if [ -f "$(basename $KLINK)" ] # tarball pulled
  then
    rm "$(basename $KLINK)"
  fi
  
  # pull fresh tarball
  mesg "Pulling database tarball..."
  wget --quiet "$KLINK"
  checkcmd "Pulling tarball"
  
  # expand the tarball
  mesg "Expanding database tarball..."
  mkdir "$(basename $KDIR)"
  tar -xf "$(basename $KLINK)" -C "$(basename $KDIR)"
  checkcmd "Expanding tarball"
  
  # clean up and move back up
  mesg "Cleaning up..."
  kraken2-build --clean --db "$(basename $KDIR)"
  checkcmd "Cleanup"
  rm "$(basename $KLINK)"
  cd ..
fi
echo ""

## make index ------------------------------------------------------------------
mesg "Creating Bowtie2 index: $BOWTIE"

# build command and execute
CMD="bowtie2-build --quiet '$FASTA' '$BOWTIE'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Bowtie2 index"
echo ""

## all done --------------------------------------------------------------------
mesg "Setup complete!"
module list
