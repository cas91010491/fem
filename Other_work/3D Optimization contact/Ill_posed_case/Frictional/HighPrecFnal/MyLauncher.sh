#!/bin/bash -l
#SBATCH -J HighPrecFnal
#SBATCH --mail-type=end
#SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=1-23:59:59
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o %x-%j.log


echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR


### load modules 


### activate python env 
conda activate n386

### run python script 
python FrictionRigidSphere.py



echo "==Ending run at $(date)"
