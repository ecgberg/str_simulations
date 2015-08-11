#!/bin/sh
#
#$ -N run_STRs_1
#
#$ -l h_vmem=4G
#
#$ -l h_rt=5:59:59
#
#$ -m eab
#
#$ -M eglassbe@stanford.edu
#
#$ -w e
#
#$ -e /home/eglassbe/strs/jobs/str_1.err
#
#$ -cwd
#
#$ -o /home/eglassbe/strs/jobs/str_1.out
#

module load python/2.7

mu=.0005
p=0.84
sigsq_g=20

beta=0.5

path='/home/eglassbe/str_simulations/str_res/'

python str_sims.py -m $mu $p $sigsq_g -b $beta -p $path
