#!/bin/bash
#SBATCH --time=04:00:00

echo  start at `date`
####TO CHANGE
phase=ISIMIP3a
scenario=historical
type=counterclim #obsclim
model=20CRv3-ERA5 #20CRv3 #20CRv3-ERA5  #20CRv3-W5E5 #GSWP3-W5E5
model_l=20crv3-era5 #20crv3 #20crv3-era5   #20crv3-w5e5 #gswp3-w5e5
year_ini=1901
year_end=2021
####TO CHANGE
root_save=$VSC_DATA_VO/vsc10623/${phase}/Output_GOTM/${type}/${scenario}/${model}/
root_gotm=$VSC_SCRATCH_VO/vsc10623/${phase}/GOTM/
root_input=$VSC_SCRATCH_VO/vsc10623/${phase}/Input/${type}/${scenario}/${model}/
root_output=$VSC_SCRATCH_VO/vsc10623/${phase}/Output_GOTM/${type}/${scenario}/${model}/

for folder in {1..100}
do

  cd ${root_gotm}/${folder}/

  mv *xr.nc ${root_output}
  mv *2021.nc ${root_save}

done

echo  end at `date`
