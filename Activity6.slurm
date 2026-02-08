#!/bin/bash
#SBATCH --partition=normal # on hopper, you should always submit your non-interactive jobs to the normal partition
#SBATCH --cpus-per-task=[number of CPUS requested for each task]
#SBATCH --job-name [a short descriptive name of the job]
#SBATCH --mem-per-cpu=[memory required per CPU]
#SBATCH --time=[time the job needs in the format days-hours:minutes]
#SBATCH --output=[path for the output log]
#SBATCH --error=[path for the error log]

# load modules and other preparations
module load gnu12 r/4.4.0-od

# run the script
Rscript Activity6.R --repl 1000 --n 100 ${SLURM_ARRAY_TASK_ID}
