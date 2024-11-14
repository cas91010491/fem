#!/bin/bash -l
#SBATCH -J Net_Parallel            # Job name
#SBATCH --mail-type=end
###SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --mem=0
#SBATCH --time=1-23:59:59
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o %x-%j.log

# Load MATLAB module (if available)
module load toolchain/intel
module load math/MATLAB/2021a

# Preload libraries
export LD_PRELOAD="/usr/lib64/libcrypto.so.1.1 /usr/lib64/libssl.so.1.1"
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Run MATLAB script and redirect output to log files
srun -n $SLURM_NTASKS --cpu-bind=no matlab -nodisplay -nodesktop -r "run('LarsNet.m'); exit" > Recov_Output.log 2> Recov_errors.log

# You can add additional commands or post-processing steps here
echo "==Ending run at $(date)"

