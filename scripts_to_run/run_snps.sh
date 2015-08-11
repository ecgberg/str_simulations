#!/bin/sh
#
#$ -N run_SNPs_1
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
#$ -e /home/eglassbe/strs/jobs/snp_1.err
#
#$ -cwd
#
#$ -o /home/eglassbe/strs/jobs/snp_1.out
#

module load python/2.7

mu=.0005

beta=1

path='/home/eglassbe/str_simulations/str_res/'

python str_sims.py -m $mu -b $beta -p $path
