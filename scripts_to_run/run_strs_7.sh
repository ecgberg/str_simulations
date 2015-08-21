#!/bin/sh
#
#$ -N run_STRs_7
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
#$ -e /home/eglassbe/str_simulations/jobs/str_7.err
#
#$ -cwd
#
#$ -o /home/eglassbe/str_simulations/jobs/str_7.out
#

module load python/2.7

mu=.0005
p=0.84
sigsq_g=20

beta=0.5

path='/home/eglassbe/str_simulations/str_res/'

python str_sims.py -m $mu $p $sigsq_g -b $beta -p $path

beta=( 0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.1 0.15 0.05 )

runnum=7

path='/home/eglassbe/str_simulations/str_res/run_'$runnum'/'

python str_sims.py -m $mu $p $sigsq_g -b ${beta[0]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[1]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[2]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[3]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[4]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[5]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[6]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[7]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[8]} -p $path -w 10
python str_sims.py -m $mu $p $sigsq_g -b ${beta[9]} -p $path -w 10
