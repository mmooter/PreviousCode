#!/bin/bash
#SBATCH --job-name=sim4cc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=150gb
#SBATCH --output=/work/pauley/shared/projects/Neandertal/Logs/sim4cc.log
#SBATCH --error=/work/pauley/shared/projects/Neandertal/Logs/sim4cc.err
python /work/pauley/shared/projects/Scripts/sim4cc.py
