#!/bin/bash -l
#SBATCH -J 2d-pl-SignedDist_10epochs
#SBATCH --mail-type=end
#SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=0-1:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o %x-%j.log

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR


# Activate python environment
conda activate fem-env

# Run Python script with the parameters
python -u SeparateTraining_SignedDistance.py

echo "==Ending run at $(date)"

