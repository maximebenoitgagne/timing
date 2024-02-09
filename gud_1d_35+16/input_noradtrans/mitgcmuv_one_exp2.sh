#!/bin/bash
# ./mitgcmuv_one_exp2.sh run_dir rbcs_file
# run_dir: run directory
# rbcs_file: data.rbcs

# run in
# /home/benoitga/projects/def-fmaps/benoitga/gud_groups/gud_1d_35+16/input_noradtrans

# set constants
date=$(date '+%Y%m%d')
BUILD_DIR="build_code_20220325_0000_EXP2_1rbcsno3times1_00_defineGUD_READ_PARUICE"
DIAGS_DIR="diags_${date}_0001"
MITGCMUV_NAME="mitgcmuv_20220325_0000_EXP2_1rbcsno3times1_00_defineGUD_READ_PARUICE"
PREVIOUS_DIR="run_20230324_0000_EXP2_1rbcsno3times1_00sioh4times1_00"

# test current directory
if [[ $PWD != */input_noradtrans ]]
then
    echo 'ERROR: This file should be run in input_noradtrans'
    exit 1
fi

# change current directory
cd ..
# read arguments
run_dir=$1
##rbcs_file=$2
# load modules for mitgcmuv
##module load nixpkgs/16.09
##module load gcc/5.4.0
##module load netcdf-fortran/4.4.4
# mitgcmuv
##cp -r input_noradtrans ${run_dir}
##cp ${BUILD_DIR}/mitgcmuv ${run_dir}/${MITGCMUV_NAME}
cd ${run_dir}
##ln -s ${MITGCMUV_NAME} mitgcmuv
##ln -fs ${rbcs_file} data.rbcs
##tar -xf input/NEMO_GE_2016.tar -C input/
##./mitgcmuv
##rm -r input/NEMO_GE_2016
# move output
##mv ${DIAGS_DIR}/*.nc .
# load modules for groups.py
module --force purge
module load StdEnv/2020
module load mii/1.1.2
module load python/3.7 scipy-stack
source ../../gud_groups/bin/activate
# groups.py
cp ../../output/${PREVIOUS_DIR}/*.py .
cp ../../output/${PREVIOUS_DIR}/*.mplstyle .
python groups.py
# copy CSV files
##mkdir comp
cp *.run_*.csv ../${run_dir}/comp
