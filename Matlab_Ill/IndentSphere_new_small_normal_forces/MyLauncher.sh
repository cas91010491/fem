#!/bin/bash -l
#SBATCH -J HexaIndent            # Job name
#SBATCH --mail-type=end
#SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH -o HexaIndent.out        # Output file for standard output
#SBATCH -e HexaIndent.err        # Error file for standard error
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=1-23:59:59
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o %x-%j.log

# Load MATLAB module (if available)
module load math/MATLAB/2021a

# Run MATLAB script and redirect output to log files
matlab -nodisplay -nodesktop -r "run('IndentSphereQuarter_rec.m'); exit" > Recov_Output.log 2> Recov_errors.log

# You can add additional commands or post-processing steps here
echo "==Ending run at $(date)"


