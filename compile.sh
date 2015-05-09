#!/bin/bash
if [ -z "$KINFIT2PATH" ];
then
    echo "KINFIT2PATH not specified."
else
    echo "creating executable"
#g++ mainTorbenTest.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -L $KINFIT2PATH -lHHKinFit2  -o runHHTauTauEventGenerator
#g++ controlplots.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit2  -o createcontrolplots
#g++ KinFitwithEventGenerator.C `root-config --cflags --glibs` -I ./include -L . -lHHKinFit2  -o KinFitwithEventGenerator
    g++ backgroundtest.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -L $KINFIT2PATH -lHHKinFit2  -o BackgroundTest
    g++ newmain.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -L $KINFIT2PATH -lHHKinFit2  -o multiplefit
 
fi