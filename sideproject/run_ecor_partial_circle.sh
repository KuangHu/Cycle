#!/bin/bash
#SBATCH --job-name=ecor_pc
#SBATCH --output=/home/kuangh/logs/ecor_pc_%j.out
#SBATCH --error=/home/kuangh/logs/ecor_pc_%j.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --qos=standard
#SBATCH --partition=standard

conda run -n opfi python /home/kuangh/tools/Cycle/sideproject/run_ecor_partial_circle.py
