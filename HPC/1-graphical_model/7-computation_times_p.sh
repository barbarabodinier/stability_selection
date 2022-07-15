#!/bin/sh
#PBS -l walltime=120:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_p
#PBS -q med-bio
#PBS -J 1-1100

cd /work/bbodinie/stability_selection/Scripts/1-graphical_model
module load anaconda3/personal

simul_study_id=1
topology="random"
seed=$PBS_ARRAY_INDEX

Rscript computation_times_p.R $simul_study_id $topology $seed
