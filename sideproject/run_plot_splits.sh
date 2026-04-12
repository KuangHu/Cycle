#!/bin/bash
#SBATCH --job-name=plot_splits
#SBATCH --output=/home/kuangh/logs/plot_splits_%j.out
#SBATCH --error=/home/kuangh/logs/plot_splits_%j.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --qos=standard
#SBATCH --partition=standard

conda run -n opfi python /home/kuangh/tools/Cycle/sideproject/plot_split_example.py
