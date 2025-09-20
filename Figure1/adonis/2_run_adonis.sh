#!/bin/bash
#SBATCH --chdir=/home/lakhatib/3country/final_scripts/adonis
#SBATCH --output=/home/lakhatib/3country/final_scripts/adonis/slurm_out/adonis_%A_%a.out
#SBATCH --mail-user="lakhatib@ucsd.edu"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --array=0-3

source ~/miniforge3/bin/activate qiime2-2023.5

DMS=("uwuf" "wuf" "rpca" "prpca")
DM="${DMS[$SLURM_ARRAY_TASK_ID]}"

python adonis_run.py --dm "$DM"
