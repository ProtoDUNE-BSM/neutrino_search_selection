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
# find *prod_protodunehd*.root | sed 's#/pnfs/dune/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/#' > $TMPFILE
# find *anaOut*.root | sed 's#/pnfs/dune/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/#' > $TMPFILE
hadd -f combined.root @${TMPFILE}


# create a copy of the combined file in the persistent storage, following the same structure as the input folder
SOURCE_BASE="/pnfs/dune/scratch/users/dpullia/"
DEST_BASE="/pnfs/dune/persistent/users/dpullia/"
REL_PATH="${FOLDER#$SOURCE_BASE}"
DEST_DIR="${DEST_BASE}${REL_PATH}"

echo "INFO  : Copying combined.root to persistent storage: $DEST_DIR"
CURRENT_DIR="$DEST_BASE"
IFS='/' read -ra PARTS <<< "$REL_PATH"
for PART in "${PARTS[@]}"; do
    [ -z "$PART" ] && continue
    CURRENT_DIR="${CURRENT_DIR%/}/$PART"
    if [ ! -d "$CURRENT_DIR" ]; then
        echo "INFO  : Creating directory $CURRENT_DIR"
        mkdir "$CURRENT_DIR"
        if [ $? -ne 0 ]; then
            echo "ERROR : Failed to create directory $CURRENT_DIR. Exiting."
            cd $THIS_DIR
            exit 1
        fi
    fi
done
cp "$FOLDER/combined.root" "$DEST_DIR/combined.root"
if [ $? -eq 0 ]; then
    echo "INFO  : Successfully copied combined.root to $DEST_DIR"
else
    echo "ERROR : Failed to copy combined.root to $DEST_DIR"
fi

cd $THIS_DIR

