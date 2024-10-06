#!/bin/bash -l
#SBATCH -J IllCondMethsForD
#SBATCH --mail-type=end
#SBATCH --mail-user=diego.hurtado@uni.lu
#SBATCH --time=0-23:59:59
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1


# Just useful check
print_error_and_exit() { echo "***ERROR*** $*"; exit 1; }
hash parallel 2>/dev/null && test $? -eq 0 || print_error_and_exit "Parallel is not installed on the system"

# Increase the user process limit.
ulimit -u 10000

echo "Node: ${SLURM_NODELIST}"
echo "Executing ${SLURM_NTASKS} independant tasks at the same time"
export TIMESTAMP=$(date +"%Y%m%dT%H%M%S")


# the --exclusive to srun makes srun use distinct CPUs for each job step
# -N1 -n1 single task with ${SLURM_CPUS_PER_TASK} cores
SRUN="srun  --exclusive -n1 -c ${SLURM_CPUS_PER_TASK:=1} --cpu-bind=cores"

HOSTNAME=$(hostname)
LOGS="logs.${TIMESTAMP}"
RESUME=""
TASKFILE="IllConditioned.py"
NTASKS=""


#=======================
# Get Optional arguments
#=======================
while [ $# -ge 1 ]; do
    case $1 in
        -r | --resume)           shift; LOGS=$1; RESUME=" --resume ";;
        -n | --ntasks)           shift; NTASKS="$1"                            ;;
        -* | --*)                echo "[Warning] Invalid option $1"          ;;
        *)                       break                                       ;;
    esac
    shift
done

#=======================
# Get Mandatory  Options
#=======================

if [[ "$#" < 1 ]]; then
    print_error_and_exit "No tasks defined"
else
    TASKFILE="$1"
    TASKFILE_DIR=$(cd "$(dirname "${TASKFILE}")" && pwd)
    TASKFILE_NAME="$(basename "${TASKFILE}")"
fi

echo "Starting parallel worker initialisation on $(hostname)"

#=======================
# Set logs for resuming
#=======================

LOGS_DIR="${TASKFILE_DIR}/${LOGS}"
TASKFILE_MAINLOG="${LOGS_DIR}/${TASKFILE_NAME//sh/log}"
PARALLEL="parallel --delay 0.2 -j ${SLURM_NTASKS} --joblog ${TASKFILE_MAINLOG} ${RESUME}"
METHS=(CG BFGS Newton-CG TNC COBYLA SLSQP dogleg trust-ncg trust-exact trust-krylov)
DISPS=(d f df)

echo "Create logs directory if not existing"
mkdir -p ${LOGS_DIR}

if [[ -z ${NTASKS} ]]; then
  parallel "${SRUN} python IllConditioned.py {} > ${LOGS_DIR}/IllConditioned.log.{#}" ::: "${METHS[@]}" ::: "${DISPS[@]}"
else
  parallel "${SRUN} python IllConditioned.py {} > ${LOGS_DIR}/IllConditioned.log.{#}" ::: "${METHS[@]}" ::: "${DISPS[@]}" ::: $(seq 1 ${NTASKS})
fi
