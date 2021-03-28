#!/bin/sh
#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_comp
#PBS -q med-bio
#PBS -J 1-1100

cd /work/bbodinie/stability_selection/Scripts/2-variable_selection
module load anaconda3/personal

simul_study_id={simul_study_id_input}
do_exp={do_exp_input}
sd_pred_error=2
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX
PFER_thr={PFER_thr_input}

Rscript comparison.R $simul_study_id $do_exp $sd_pred_error $params_id $seed $PFER_thr
