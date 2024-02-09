#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2048M
#SBATCH --account=def-fmaps
#SBATCH --job-name=mitgcmuv_all
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxime.benoit-gagne@takuvik.ulaval.ca
#SBATCH --array=1-122

# run in
# /home/benoitga/projects/def-fmaps/benoitga/gud_groups/gud_1d_35+16/input_noradtrans

echo "Starting task ${SLURM_ARRAY_TASK_ID}"

# set constants
CSV_DIR="run_20230426_0000_EXP2_3_122_sims_par"
RBCS_FILES="rbcs_files.txt"
RUN_DIR_PREFIX="run_20230426_0000_EXP2_3"
RUN_DIRS="run_dirs_exp2.txt"

# test current directory
if [[ $PWD != */input_noradtrans ]]
then
    echo 'ERROR: This file should be run in input_noradtrans'
    exit 1
fi

# run mitgcmuv
run_dir=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RUN_DIRS})
rbcs_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RBCS_FILES})
./mitgcmuv_one_exp2.sh "${run_dir}" "${rbcs_file}"

# gather CSV files
cd ..
mkdir -p "${CSV_DIR}"
for f in "${RUN_DIR_PREFIX}"*/comp/*.csv
do
    cp "$f" "${CSV_DIR}"
done
