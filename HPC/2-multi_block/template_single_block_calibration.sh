#!/bin/sh
#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_single_block
#PBS -q med-bio
#PBS -J 1-1100

cd /work/bbodinie/stability_selection/Scripts/2-multi_block
module load anaconda3/personal

simul_study_id={simul_study_id_input}
topology={topology_input}
do_exp={do_exp_input}
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX
PFER_thr={PFER_thr_input}

Rscript single_block_calibration.R $simul_study_id $topology $do_exp $params_id $seed $PFER_thr
