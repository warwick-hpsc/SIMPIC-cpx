#!/bin/bash
#######run the following command to display a full description of memory issues:
#cuda-memcheck ./runsimpic-1.sh
####or
#cuda-memcheck ./runsimpic-1.sh |more

NT=200  # of time steps
PPC=200 # of particles per cell
CPP=500  # of cell per process
DTFACTOR=0.001 #defines a fraction of dt vs. time steps from plasma frequency;
               #must be positive
DIAGSTEPS=1 #steps saving diags
LHSV=25000 #applied voltage on left-hand side; RHS is grounded; 
PROGNAME=a.out

CMDLINE="mpirun -np 1 $PROGNAME -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV -diagsteps $DIAGSTEPS"
echo $CMDLINE
$CMDLINE
