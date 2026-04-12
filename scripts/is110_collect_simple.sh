#!/bin/bash
#SBATCH --job-name=is110_collect
#SBATCH --output=/home/kuangh/logs/is110_collect_%j.out
#SBATCH --error=/home/kuangh/logs/is110_collect_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --qos=standard
#SBATCH --partition=standard

conda run -n opfi python /home/kuangh/tools/Cycle/scripts/is110_collect_simple.py
