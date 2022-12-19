#!/bin/sh
#PBS -l walltime=120:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N simul_p
#PBS -q med-bio
#PBS -J 1-1100

cd /work/bbodinie/stability_selection/Scripts/3-regression_model
module load anaconda3/personal
source activate selection

simul_study_id=1
seed=$PBS_ARRAY_INDEX

Rscript computation_times_p.R $simul_study_id $seed
