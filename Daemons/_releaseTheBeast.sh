#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=AI
#SBATCH --job-name=DaemonDB
#SBATCH --comment=""
#SBATCH --output=Daemons.out
#SBATCH --error=Daemons.err
## --mem=0

source /home/$USER/.bashrc

srun beelzebul.ru