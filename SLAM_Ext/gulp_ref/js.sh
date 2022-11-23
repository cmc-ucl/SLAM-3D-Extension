#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N task_!
#$ -e task_!.err
#$ -o task_!.out
#$ -pe smp 1

module purge
module load gulp/gnu-4.9/5.0
/opt/openmpi/gfortran/1.8.4/bin/mpirun -n $NSLOTS gulp-mpi < single.in > single.in.got

