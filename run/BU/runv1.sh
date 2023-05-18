#!/bin/bash

echo empieza 
module load R/4.1.2-foss-2021b
#module load CDO/2.0.5-gompi-2021b
#module load GLM-AED/3.3.0a5-gompi-2021b
module load GOTM/5.4.0-20210127-gompi-2021b-lake-stim

echo  start at `date`
####TO CHANGE
phase=ISIMIP3a
scenario=historical
type=obsclim
model=GSWP3-W5E5
model_l=gswp3-w5e5
year_ini=1901
year_end=2019
####TO CHANGE

root_gotm=$VSC_SCRATCH_VO/vsc10623/${phase}/GOTM/
root_input=$VSC_DATA_VO/vsc10623/${phase}/Input/${type}/${scenario}/${model}/
root_output=$VSC_SCRATCH_VO/vsc10623/${phase}/Output_GOTM/${type}/${scenario}/${model}/
pixel=1
folder=1
echo pixel: ${pixel} in folder ${folder}


cp hypso_gotm/h_${pixel}.dat ${root_gotm}/${folder}/hypsograph.dat
cp ${root_input}/${model_l}_${scenario}_${type}_lake_${pixel}_daily_${year_ini}_${year_end}.csv ${root_gotm}/${folder}/meteo_file.dat
cp gotm.yaml ${root_gotm}/${folder}
cp config_gotm_yaml_spinup.R ${root_gotm}/${folder}

cd ${root_gotm}/${folder}/

sed -i "2s/.*/pixel <- ${pixel}/" config_gotm_yaml_spinup.R
sed -i "3s/.*/year_ini <- ${year_ini}/" config_gotm_yaml_spinup.R
sed -i "4s/.*/year_end <- ${year_end}/" config_gotm_yaml_spinup.R
Rscript config_gotm_yaml_spinup.R
sed -i "288s/.*/   ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}:/" gotm.yaml

gotm

cp ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}.nc ${root_output}

echo final final
