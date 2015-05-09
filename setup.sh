#!/bin/bash
export KINFIT2PATH=/afs/desy.de/user/t/tlange/workspace/HHKinFit2
export LD_LIBRARY_PATH=.:$KINFIT2PATH:$LD_LIBRARY_PATH
echo "KinFit path set to $KINFIT2PATH"

if [[ $HOST == *"nafhh"* ]]
then
    module load git
    module load root
    echo "ROOT and Git initialized via module"
else
    echo "Root and Git not initialized. Please make sure that git and root is initialized"
fi