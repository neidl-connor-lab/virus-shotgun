#!/bin/bash -l

# qsub options
#$ -l h_rt=12:00:00
#$ -l mem_total=1000G
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
BDIR="pipeline/indices"
LOFREQ="pipeline/lofreq"
KDIR="pipeline/kraken2db"
KLINK="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20221209.tar.gz"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index name
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
# check for Kraken2 install
if [ -z "$(which kraken2 2> /dev/null)" ]
then
  err "No Kraken2 detected. Please install before using."
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

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## unpack LoFreq if necessary -----------------------------------
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

## set up Kraken2 if necessary ----------------------------------
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
  cd "$(dirname $KDIR)"
  if [ -f "$(basename $KLINK)" ] # tarball pulled
  then
    # remove tarball and pull fresh
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
  kraken2-build --clean --db "$KDIR"
  rm "$(basename $KLINK)"
  cd ..
fi
echo ""

## make index ---------------------------------------------------
mesg "Creating Bowtie2 index: $BDIR/$BOWTIE"

# create directory for bowtie2 if it doesn't already exist
if [ ! -d "$BDIR" ]
then
  mkdir -p "$BDIR"
fi

# build command and execute
CMD="bowtie2-build --quiet '$FASTA' '$BDIR/$BOWTIE'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Bowtie2 index"
echo ""

## all done -----------------------------------------------------
mesg "Setup complete!"
module list
