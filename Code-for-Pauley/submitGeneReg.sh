#!/bin/bash
#SBATCH --job-name=geneReg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=150gb
#SBATCH --output=/work/pauley/shared/projects/Denisova/Logs/geneReg.log
#SBATCH --error=/work/pauley/shared/projects/Denisova/Logs/geneReg.err
python /work/pauley/shared/projects/Scripts/createGeneRegionsFromSim4cc.py
