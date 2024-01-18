#!/bin/bash

echo empieza 
module load R/4.1.2-foss-2021b
module load GOTM/5.4.0-20210127-gompi-2021b-lake-stim
ml NCO/5.0.3-foss-2021b
echo  start at `date`
####TO CHANGE
phase=ISIMIP3b
scenario=picontrol #ssp585 #ssp370 #ssp126 #historical
type=bias-adjusted # counterclim #obsclim
model=MRI-ESM2-0 #ISIMIP3b: MRI-ESM2-0 #MPI-ESM1-2-HR #IPSL-CM6A-LR #UKESM1-0-LL #GFDL-ESM4 # ISIMIP3a: 20CRv3 #20CRv3-ERA5  #20CRv3-W5E5 #GSWP3-W5E5
model_l=mri-esm2-0 #ISIMIP3b: mri-esm2-0 #mpi-esm1-2-hr #ipsl-cm6a-lr #ukesm1-0-ll #gfdl-esm4 #20crv3 #20crv3-era5   #20crv3-w5e5 #gswp3-w5e5
year_ini=1841 #picontrol: 1841, historical: 1850, ssp: 2015
year_end=2100 #picontrol: 20100, historical: 2014, ssp: 2100
####TO CHANGE
root_save=$VSC_DATA_VO/vsc10623/${phase}/Output_GOTM/${type}/${scenario}/${model}/
root_gotm=$VSC_SCRATCH_VO/vsc10623/${phase}/GOTM/
root_input=$VSC_SCRATCH_VO/vsc10623/${phase}/Input/${type}/${scenario}/${model}/
root_output=$VSC_SCRATCH_VO/vsc10623/${phase}/Output_GOTM/w5e5/${scenario}/${model}/
folder=11

cp $VSC_SCRATCH/run/gotm.yaml ${root_gotm}/${folder}
cp $VSC_SCRATCH/run/config_gotm_yaml_spinup.R ${root_gotm}/${folder}

for pixel in {4151..4565}
do
  echo pixel: ${pixel} in folder ${folder}


  cp $VSC_SCRATCH/run/hypso_gotm/h_${pixel}.dat ${root_gotm}/${folder}/hypsograph.dat
  cp ${root_input}/${model_l}_${scenario}_${type}_lake_${pixel}_daily_${year_ini}_${year_end}.csv ${root_gotm}/${folder}/meteo_file.dat

  cd ${root_gotm}/${folder}/

  sed -i "2s/.*/pixel <- ${pixel}/" config_gotm_yaml_spinup.R
  sed -i "3s/.*/year_ini <- ${year_ini}/" config_gotm_yaml_spinup.R
  sed -i "4s/.*/year_end <- ${year_end}/" config_gotm_yaml_spinup.R
  Rscript config_gotm_yaml_spinup.R
  sed -i "288s/.*/   ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}:/" gotm.yaml

  gotm
  
  ncrename -v z,z_coord ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}.nc ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}xr.nc
  mv ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}xr.nc ${root_output}
  #mv ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}.nc ${root_save}
  rm ${model_l}_${scenario}_${type}_gotm_${pixel}_daily_${year_ini}_${year_end}.nc

done

echo  end at `date`
