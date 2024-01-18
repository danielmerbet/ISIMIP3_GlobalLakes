#!/bin/bash

#SBATCH --job-name=ISIMIP3b-pre
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00

###SBATCH --ntasks=28
###SBATCH --partition=broadwell
###SBATCH --reservation=gotm

ml SciPy-bundle/2021.10-foss-2021b
ml netcdf4-python/1.5.8-foss-2021b
ml xarray/2022.6.0-foss-2021b

echo ./extract_lakes_xarray_GOTMv1_12h.py --phase ISIMIP3b -t bias-adjusted -m IPSL-CM6A-LR -c ssp370 --base $VSC_DATA_VO/data/dataset/ISIMIP --out $VSC_SCRATCH_VO/vsc10623/ISIMIP3b/Input/bias-adjusted/ssp370/IPSL-CM6A-LR
./extract_lakes_xarray_GOTMv1_12h.py --phase ISIMIP3b -t bias-adjusted -m IPSL-CM6A-LR -c ssp370 --base $VSC_DATA_VO/data/dataset/ISIMIP --out $VSC_SCRATCH_VO/vsc10623/ISIMIP3b/Input/bias-adjusted/ssp370/IPSL-CM6A-LR
