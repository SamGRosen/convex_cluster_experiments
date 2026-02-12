#!/bin/bash
#SBATCH --output=./slurm-output/%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-400
#SBATCH --account=dctrl-sgr26
#SBATCH --exclude=dcc-core-[30-40]


module load R/4.4.0
module load GLPK/5.0
module load OpenBLAS/3.23

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/

pwd

echo $SLURM_JOB_ID
echo $SLURM_ARRAY_JOB_ID
echo $SLURM_ARRAY_TASK_ID

Rscript /hpc/home/sgr26/convex_clustering_experiments/trial.R $SLURM_ARRAY_JOB_ID \
  $SLURM_ARRAY_TASK_ID
