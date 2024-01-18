#!/bin/bash

#SBATCH --job-name=ISIMIP3b-pos
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --time=05-00:00:00
#SBATCH --mem=180G


###SBATCH --partition=broadwell
###SBATCH --ntasks=28
###SBATCH --reservation=gotm


##SBATCH --ntasks=40
##SBATCH --partition=skylake
##SBATCH --mem=180G


echo EMPIEZA TODO `date`

ml SciPy-bundle/2021.10-foss-2021b
ml netcdf4-python/1.5.8-foss-2021b
ml xarray/2022.6.0-foss-2021b
ml NCO/5.0.3-foss-2021b
ml CDO/2.0.5-gompi-2021b

echo create_ISIMIP_NetCDF_fullprofile_IPSL.py -p ISIMIP3b -t bias-adjusted -m IPSL-CM6A-LR -c historical -fi 1850 -f 2014
./create_ISIMIP_NetCDF_fullprofile_IPSL.py -p ISIMIP3b -t bias-adjusted -m IPSL-CM6A-LR -c historical -fi 1850 -f 2014

echo TERMINA TODO `date`

