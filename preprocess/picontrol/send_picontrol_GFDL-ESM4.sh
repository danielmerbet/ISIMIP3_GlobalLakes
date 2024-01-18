#!/bin/bash

#SBATCH --job-name=ISIMIP3bP
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00

###SBATCH --partition=broadwell
###SBATCH --mem=230G
###SBATCH --ntasks=28


###SBATCH --reservation=gotm

ml SciPy-bundle/2021.10-foss-2021b
ml netcdf4-python/1.5.8-foss-2021b
ml xarray/2022.6.0-foss-2021b

echo ./extract_lakes_xarray_GOTMv1.py --phase ISIMIP3b -t bias-adjusted -m GFDL-ESM4 -c picontrol --base $VSC_DATA_VO/data/dataset/ISIMIP --out $VSC_SCRATCH_VO/vsc10623/ISIMIP3b/Input/bias-adjusted/picontrol/GFDL-ESM4
./extract_lakes_xarray_GOTMv1.py --phase ISIMIP3b -t bias-adjusted -m GFDL-ESM4 -c picontrol --base $VSC_DATA_VO/data/dataset/ISIMIP --out $VSC_SCRATCH_VO/vsc10623/ISIMIP3b/Input/bias-adjusted/picontrol/GFDL-ESM4
