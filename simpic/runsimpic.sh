#!/bin/bash

#SBATCH --job-name=cupcfd_mpi     # Job name
#SBATCH --ntasks=768                # Run all processes on a single node	
#SBATCH --cpus-per-task=1
#SBATCH --nodes=6
#SBATCH --time=00:20:00
#SBATCH --output=simpic_%j.log
#SBATCH --partition=standard
#SBATCH --qos=short

export CC=cc
export CXX=CC
export FC=ftn

cd /work/e647/e647/awp/SIMPIC-cpx/simpic/mpi


NT=50000  # of time steps
PPC=100 # of particles per cell
CPP=666 # of cell per process
DTFACTOR=0.000001 #defines a fraction of dt vs. time steps from plasma frequency;
               #must be positive
NPROC=4 # of processors
LHSV=20000 #applied voltage on left-hand side; RHS is grounded;
PROGNAME=a.out

CMDLINE="-ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"


srun --unbuffered --distribution=block:block --hint=nomultithread ./a.out $CMDLINE
