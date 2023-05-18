#!/bin/bash

#SBATCH --job-name=ISIMIP3a-pre-GSWP3
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --partition=broadwell
#SBATCH --mem=230G
#SBATCH --time=05-00:00:00


###SBATCH --reservation=gotm

ml SciPy-bundle/2021.10-foss-2021b
ml netcdf4-python/1.5.8-foss-2021b
ml xarray/2022.6.0-foss-2021b

echo ./extract_lakes_xarray_GOTMv1_counterclim.py --phase ISIMIP3a -t counterclim -m GSWP3-W5E5 -c historical

./extract_lakes_xarray_GOTMv1_counterclim.py --phase ISIMIP3a -t counterclim -m GSWP3-W5E5 -c historical --base $VSC_DATA_VO/data/dataset/ISIMIP --out $VSC_SCRATCH_VO/vsc10623/ISIMIP3a/Input/counterclim/historical/GSWP3-W5E5
