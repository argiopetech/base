#!/bin/csh

set nProc = 2

echo -n " time stamp: "; date; echo ""
mpirun -np $nProc bin/mpiMcmc

echo -n " time stamp: "; date; echo ""
mpirun -np $nProc bin/sampleWDMass

echo -n " time stamp: "; date
