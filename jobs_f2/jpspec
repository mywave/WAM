#!/bin/bash
#
#     Arno Behrens (May 2019)
#
#==>  WAM post-processing pspec
#
#SBATCH --job-name=pspec
#SBATCH --partition=pCluster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00
#SBATCH --mail-type=FAIL
#SBATCH --account=cluster
#SBATCH --output=pspec.o%j
#SBATCH --error=pspec.e%j
#
module load compilers/intel/2019.4.243
module load intelmpi/2019.4.243
module load netcdf
#
set -k
WAMDIR=WAMDIR=/gpfs/home/behrens/WAM_Cycle_6
WORK=/gpfs/work/behrens/WAM_Cycle_6
#
cd ${WORK}/tempsg
cp ${WAMDIR}/const/Coarse_Grid/JAN/Spectra_User .
cp ${WAMDIR}/abs/pspec pspec.exe
#
./pspec.exe
#
mv Spectra_Prot ${WAMDIR}/dayfiles/coarse/pspec_prot_coarse
rm Spectra_User pspec.exe
#
# ===================================================================
#  GRID FILES HAVE BEEN CREATED AND SAVED.
#  END OF JOB PGRID.
# ===================================================================
#
