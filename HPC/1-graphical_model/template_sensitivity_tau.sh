#!/bin/sh
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_tau
#PBS -q med-bio
#PBS -J 1-1100

cd /work/bbodinie/stability_selection/Scripts/1-graphical_model
module load anaconda3/personal
source activate selection

simul_study_id={simul_study_id_input}
topology={topology_input}
do_exp={do_exp_input}
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX

Rscript sensitivity_tau.R $simul_study_id $topology $do_exp $params_id $seed
