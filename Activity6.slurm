#!/bin/bash
#SBATCH --job-name=stat778_activity6
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=0-02:00:00
#SBATCH --output=rscript_%a.out
#SBATCH --error=rscript_%a.err

module purge
module load hosts/hopper
module load gnu12/12.3.0
module load openmpi4/4.1.6-a2
module load r/4.4.0-od

which Rscript
Rscript --version

Rscript Activity6.R --repl 1000 --n 100 ${SLURM_ARRAY_TASK_ID}
