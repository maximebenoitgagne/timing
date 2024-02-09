#!/bin/bash
# ./mitgcmuv_one_exp1.sh run_dir gud_file
# run_dir: run directory
# gud_file: data.gud

# run in
# /home/benoitga/projects/def-fmaps/benoitga/gud_groups/gud_1d_35+16/input_noradtrans

# set constants
date=$(date '+%Y%m%d')
BUILD_DIR="build_code_20220325_0000_EXP2_1rbcsno3times1_00_defineGUD_READ_PARUICE"
DIAGS_DIR="diags_${date}_0001"
MITGCMUV_NAME="mitgcmuv_20220325_0000_EXP2_1rbcsno3times1_00_defineGUD_READ_PARUICE"
PREVIOUS_DIR="run_20230707_0000_EXP0_translucent_snow"

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
gud_file=$2
# load modules for mitgcmuv
module load nixpkgs/16.09
module load gcc/5.4.0
module load netcdf-fortran/4.4.4
# mitgcmuv
rsync -a --exclude='*.out' input_noradtrans/ ${run_dir}
# copy executable
cp ${BUILD_DIR}/mitgcmuv ${run_dir}/${MITGCMUV_NAME}
# prepare run
cd ${run_dir}
ln -s ${MITGCMUV_NAME} mitgcmuv
ln -fs ${gud_file} data.gud
tar -xf input/NEMO_GE_2016.tar -C input/
# run
./mitgcmuv
# delete a directory that was unarchived
rm -r input/NEMO_GE_2016
# copy post-processing scripts
cp -r ../../output/${PREVIOUS_DIR}/comp ${DIAGS_DIR}
rm -f ${DIAGS_DIR}/comp/*.png
cp ../../output/${PREVIOUS_DIR}/*.ipynb ${DIAGS_DIR}
cp ../../output/${PREVIOUS_DIR}/*.mplstyle ${DIAGS_DIR}
cp ../../output/${PREVIOUS_DIR}/*.py ${DIAGS_DIR}
# move output
cp data data.cal data.diagnostics data.exf data.gchem data.gud data.kpp data.mnc data.off data.pkg data.ptracers data.rbcs data.traits gud_params.txt gud_traits.txt ${DIAGS_DIR} 
