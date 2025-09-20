#!/bin/bash
#SBATCH --chdir=/home/lakhatib/3country/final_scripts/adonis
#SBATCH --output=/home/lakhatib/3country/final_scripts/adonis/slurm_out/%x_%A_%a.out
#SBATCH --mail-user="lakhatib@ucsd.edu"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

source ~/miniforge3/bin/activate qiime2-2023.5


python alpha_beta_calc.py