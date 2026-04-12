#!/bin/bash
#SBATCH --job-name=is110_hmm
#SBATCH --output=/home/kuangh/logs/is110_hmm_%j.out
#SBATCH --error=/home/kuangh/logs/is110_hmm_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --qos=standard
#SBATCH --partition=standard

conda run -n opfi python /home/kuangh/tools/Cycle/scripts/is110_hmm_only.py
