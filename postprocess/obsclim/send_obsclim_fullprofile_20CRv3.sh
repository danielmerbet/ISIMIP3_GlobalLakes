#!/bin/bash

#SBATCH --job-name=ISIMIP3a-pos-GSWP3
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --mem=230G
#SBATCH --partition=broadwell
#SBATCH --ntasks=28
#SBATCH --time=05-00:00:00

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

echo create_ISIMIP_NetCDF_fullprofile.py -p ISIMIP3a -t obsclim -m 20CRv3 -c historical -fi 1901 -f 2015
./create_ISIMIP_NetCDF_fullprofile.py -p ISIMIP3a -t obsclim -m 20CRv3 -c historical -fi 1901 -f 2015

echo TERMINA TODO `date`

