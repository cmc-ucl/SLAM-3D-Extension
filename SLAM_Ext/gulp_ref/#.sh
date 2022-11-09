#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N #
#$ -e #.err
#$ -o #.out
#$ -pe smp # Keywords:

module purge
module load gulp/gnu-4.9/5.0
/opt/openmpi/gfortran/1.8.4/bin/mpirun -n $NSLOTS gulp-mpi < # > #.got

