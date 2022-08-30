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
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") -r RESOURCES

arguments:
  -r resources directory (default: $RESOURCES)
  -h show this message and exit"

# parsing arguments
while getopts ":hr:" opt
do
  case ${opt} in
    r ) RESOURCES="${OPTARG}"
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

# create directory for resources
if [ -z "$RESOURCES" ]
then
  err "No resources directory provided"
elif [ -d "$RESOURCES" ]
then
  mesg "Valide resources directory: $RESOURCES"
else
  mesg "Creating resource directory: $RESOURCES"
  mkdir -p "$RESOURCES"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## install Kraken2 ----------------------------------------------
mesg "STEP 1: SET UP KRAKEN2"

# if Kraken2 is not installed, install it
if [ -z "$(which kraken2 2> /dev/null)" ]
then
  mesg "Kraken2 installation not detected; installing: $HOME/kraken2"
  
  # download and extract Kraken2 v2.1.2
  wget --quiet "https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz"
  tar -xf "v2.1.2.tar.gz"
  
  # clean up tarball
  rm "v2.1.2.tar.gz"
  
  # install to ~/kraken2; main scripts at ~/bin
  cd "kraken2-2.1.2"
  mkdir -p "$HOME/kraken2"
  ./install_kraken2.sh $HOME/kraken2
  cp $HOME/kraken2/kraken2{,-build,-inspect} $HOME/bin
  cd ..
  rm -r "kraken2-2.1.2"
  
  # check again for Kraken2 in $PATH
  # fail if not found
  if [ -z "$(which kraken2 2> /dev/null)" ]
  then
    err "Attempted Kraken2 installation failed"
  fi
  mesg "Kraken2 installed!"
else
  mesg "Kraken2 detected: $(which kraken2)"
fi

# install Kraken2 standard database
KDB="$RESOURCES/kraken2db"
mesg "Building Kraken2 database: $KDB"
mkdir -p "$KDB"
CMD="kraken2-build --threads 16 --standard --db $KDB"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Building Kraken2 database"

# clean up Kraken2's mess
mesg "Clean up installation"
CMD="kraken2-build --threads 16 --clean --db $KDB"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Cleaning Kraken2 database"

# done with Kraken2 check
mesg "Kraken2 setup complete!"
echo ""

## install vodka ------------------------------------------------
mesg "STEP 2: SET UP VODKA"

# install weblogo
mesg "Installing weblogo"
CMD="pip3 install --user weblogo"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "weblogo installation"

# download VODKA and set up
mesg "Installing VODKA to $RESOURCES"
wget --quiet "https://github.com/itmat/VODKA/archive/refs/tags/v0.1f.tar.gz"
tar -xf "v0.1f.tar.gz"
rm "v0.1f.tar.gz"
mv "VODKA-0.1f/scripts" "$RESOURCES/vodka"
rm -r "VODKA-0.1f"

# done with VODKA setup
mesg "VODKA setup complete!"
echo ""

## all done -----------------------------------------------------
mesg "Setup complete!"
module list
echo ""

