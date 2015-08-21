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
#$ -e /home/eglassbe/str_simulations/jobs/snp_1.err
#
#$ -cwd
#
#$ -o /home/eglassbe/str_simulations/jobs/snp_1.out
#

module load python/2.7

mu=.0005

beta=( 1 0.9 0.8 0.7 0.6 0.5 0.4 0.2 0.3 0.1 )

runnum=1

path='/home/eglassbe/str_simulations/str_res/run_'$runnum'/'

python str_sims.py -m $mu -b ${beta[0]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[1]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[2]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[3]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[4]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[5]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[6]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[7]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[8]} -p $path -w 10
python str_sims.py -m $mu -b ${beta[9]} -p $path -w 10
