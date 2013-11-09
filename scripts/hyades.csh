#!/bin/csh

set nProc = 2

echo -n " time stamp: "; date; echo ""
./bin/mpiMcmc

echo -n " time stamp: "; date; echo ""
./bin/sampleWDMass

echo -n " time stamp: "; date
