#!/bin/bash
#SBATCH --job-name=gogogo    # Job name
#SBATCH --mail-user=chongmin@andrew.cmu.edu     # Where to send mail
#SBATCH --time=24:00:00
#SBATCH -N 1                 # Run on a single CPU
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --ntasks-per-node 20
#SBATCH --output=gogogo.log   # Standard output and error log

. $HOME/.bashrc
conda activate eclip
