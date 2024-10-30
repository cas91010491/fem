#!/bin/bash -l
#SBATCH -J 3d-LBFGS-500
#SBATCH --mail-type=end
#SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=1-23:59:59
#SBATCH -p batch
#SBATCH --array=0-5
#SBATCH --qos=normal
#SBATCH -o %x-%j.log

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# Define possible values for each argument
min_methods=("LBFGS500")
meshes=(5 10 15)
plastics=(0 1)  #  0:elastic, 1:plastic

# Calculate total number of combinations
total_combinations=$(( ${#min_methods[@]} * ${#meshes[@]} * ${#plastics[@]} ))

# Ensure the job array range matches the total combinations
if (( SLURM_ARRAY_TASK_ID >= total_combinations )); then
    echo "Task ID ${SLURM_ARRAY_TASK_ID} is out of range. Exiting."
    exit 1
fi

# Determine the index values for each parameter based on the task ID
method_index=$(( SLURM_ARRAY_TASK_ID / (${#meshes[@]} * ${#plastics[@]}) ))
mesh_index=$(( (SLURM_ARRAY_TASK_ID / ${#plastics[@]}) % ${#meshes[@]} ))
plastic_index=$(( SLURM_ARRAY_TASK_ID % ${#plastics[@]} ))

# Set each parameter based on the calculated indices
min_method=${min_methods[$method_index]}
mesh=${meshes[$mesh_index]}
plastic=${plastics[$plastic_index]}

# Activate python environment
conda activate fem-env

# Run Python script with the parameters
python -u ContactPotato_slideX.py --min_method ${min_method} --mesh ${mesh} --plastic ${plastic}

echo "==Ending run at $(date)"

