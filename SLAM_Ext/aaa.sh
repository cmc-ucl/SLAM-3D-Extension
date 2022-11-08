#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N aaa
#$ -e aaa.err
#$ -o aaa.out
#$ -pe smp 1


/home/uccawkj/Work/SLAM-3D-Extension/SLAM_Ext/out.x > test.out
