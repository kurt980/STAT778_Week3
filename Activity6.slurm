#!/bin/bash
#SBATCH --job-name=stat778_activity6
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=0-02:00:00
#SBATCH --output=rscript_%a.out
#SBATCH --error=rscript_%a.err

module load gnu12 r/4.4.0-od

Rscript Activity6.R --repl 1000 --n 100 ${SLURM_ARRAY_TASK_ID}
