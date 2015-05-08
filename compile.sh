#!/bin/bash
if [ -z "$KINFITPATH" ];
then
    echo "KINFITPATH not specified."
else
    echo "creating executable"
#g++ mainTorbenTest.C `root-config --cflags --glibs` -I $KINFITPATH/include -L $KINFITPATH -lHHKinFit  -o runHHTauTauEventGenerator
#g++ controlplots.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o createcontrolplots
#g++ KinFitwithEventGenerator.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit  -o KinFitwithEventGenerator
    g++ backgroundtest.C `root-config --cflags --glibs` -I $KINFITPATH/include -L $KINFITPATH -lHHKinFit  -o BackgroundTest
#g++ -std=c++11 compareKinFits.C `root-config --cflags --glibs` -I ./include -I ../HHKinFit/interface -L . -L ../HHKinFit -lHHKinFit2 -lHHKinFit -o compareKinFits
    g++ newmain.C `root-config --cflags --glibs` -I $KINFITPATH/include -L $KINFITPATH -lHHKinFit  -o multiplefit
 
fi