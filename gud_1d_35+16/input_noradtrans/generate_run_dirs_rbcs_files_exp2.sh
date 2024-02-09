#!/bin/bash

# set constants
declare -a NO3_LIST=($(seq -f %02g 0 2 20))
RBCS_FILES="rbcs_files.txt"
RUN_DIRS="run_dirs_exp2.txt"
RUN_DIR_PREFIX="run_20230426_0000_EXP2_3"
RUN_DIR_SUFFIX="_122_sims_par"
declare -a SIOH4_LIST=($(seq -f %02g 0 2 20))

rm -f "${RUN_DIRS}"
rm -f "${RBCS_FILES}"

# function to write files
write_files () {
    run_dir="${RUN_DIR_PREFIX}no3_${no3}_si_${sioh4}${RUN_DIR_SUFFIX}"
    echo "${run_dir}" >> "${RUN_DIRS}"
    rbcs_file="data.rbcs.no3.${no3}.si.${sioh4}"
    echo "${rbcs_file}" >> "${RBCS_FILES}"
}

# standard simulation
for no3 in "standard"
do
    for sioh4 in "standard"
    do
	write_files
    done
done

# loop over nutrients
for no3 in "${NO3_LIST[@]}"
do
    for sioh4 in "${SIOH4_LIST[@]}"
    do
	write_files
    done
done
