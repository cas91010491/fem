#!/bin/bash -l
#SBATCH -J SLSQP-BFGS_perf
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


### activate python env 
conda activate n386

### run python script 
# Define the values for each parameter
methods=("BFGS" "SLSQP" "L-BFGS-B")
actions=("d" "f" "df")

# Loop through each combination
for method in "${methods[@]}"; do
  for action in "${actions[@]}"; do
    # Run Python script with the current combination of parameters
    python IllConditioned.py "$method" "$action"
  done
done

echo "==Ending run at $(date)"
