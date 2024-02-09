#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=def-fmaps
#SBATCH --job-name=mitgcmuv_all
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxime.benoit-gagne@takuvik.ulaval.ca
#SBATCH --array=1-5

# run in
# /home/benoitga/projects/def-fmaps/benoitga/gud_groups/gud_1d_35+16/input_noradtrans

echo "Starting task ${SLURM_ARRAY_TASK_ID}"

# set constants
RBCS_FILES="rbcs_files_exp1_3.txt"
RUN_DIRS="run_dirs_exp1_3.txt"

# test current directory
if [[ $PWD != */input_noradtrans ]]
then
    echo 'ERROR: This file should be run in input_noradtrans'
    exit 1
fi

# run mitgcmuv
run_dir=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RUN_DIRS})
rbcs_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RBCS_FILES})
./mitgcmuv_one_exp1_3.sh "${run_dir}" "${rbcs_file}"
