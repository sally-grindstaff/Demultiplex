#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --time=1-1:00:00
#SBATCH --nodes=1

conda activate bgmp_py39

/usr/bin/time -v ./qual_scores.py

exit