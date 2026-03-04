#!/bin/bash

# This script is used to hadd the output of the postprocessor
FOLDER=$1

# Check if the folder exists
if [ ! -d $FOLDER ]; then
    echo "The folder $FOLDER does not exist. Exiting."
    exit 1
fi

THIS_DIR=$PWD
cd $FOLDER
# remote_hadd 'analysisOutput_*.root' combined.root
echo "INFO  : Hadding files in $FOLDER"
echo "INFO  : I'm in folder $PWD"

voms-proxy-init -rfc -noregen -voms=dune:/dune/Role=Analysis -valid 120:00
TMPFILE=$(mktemp)
find *proc_bsmtrigger*.root | sed 's#/pnfs/dune/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/#' > $TMPFILE
# find *anaOut*.root | sed 's#/pnfs/dune/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/#' > $TMPFILE
hadd -f combined.root @${TMPFILE}


cd $THIS_DIR

