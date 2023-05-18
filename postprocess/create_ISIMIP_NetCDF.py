#!/usr/bin/env python3
# run with Python 3.9

"""

functions to convert individual nc files into final netCDF for ISIMIP3

"""

# import modules 
import argparse
import os
import fnmatch
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
import datetime
import pylake

parser = argparse.ArgumentParser(description='Extract climate data for local lakes sector.')

parser.add_argument('-p', '--phase', dest='phase', required=True,
                    default='ISIMIP3b',
                    help='ISIMIP phase.')
parser.add_argument('-t', '--datatype', dest='datatype', required=True,
                    help='Input data types, e.g. 3a: [obsclim|counterclim], 3b: [bias-corrected]')
parser.add_argument('-m', '--model', dest='model', required=True,
                    help='Input model to process')
parser.add_argument('-c', '--climate-forcing', dest='climforcing', required=True,
                    help='climate forcing model to process, e.g [picontrol|historical|ssp126|ssp370|ssp585].')
parser.add_argument('-b', '--base', dest='basedir',
                    default='/scratch/brussel/vo/000/bvo00012/vsc10623/',
                    help='base directory')
parser.add_argument('-o', '--out', dest='outdir',
                    default='/data/brussel/vo/000/bvo00012/vsc10623/',
                    help='base directory to post-processed data')
parser.add_argument('-f', '--fye', dest='finalyear',
                    default=2019,
                    help='the final year of all the simulation')
parser.add_argument('-fi', '--fiye', dest='firstyear',
                    default=1901,
                    help='the final year of all the simulation')

args = parser.parse_args()


lakes_df = pd.read_csv('/scratch/brussel/106/vsc10623/preprocess/coord_area_depth.csv')

#variable_name='temp'
#variable_name=args.variable
variables=['temp', 'qe', 'qh', 'Hice', 'Tice_surface']
variables_save=['surftemp', 'bottemp', 'latentheatf', 'sensheatf', 'icethick', 'icetemp']
units = ['K', 'K', 'W m-2', 'W m-2', 'm', 'K', 'm', 'm']
longnames = ['Temperature of Lake Surface Water', 'Temperature of Lake Bottom Water', 'Latent Heat Flux at Lake-Atmosphere Interface', 'Sensible Heat Flux at Lake-Atmosphere Interface', 'Ice Thickness', 'Ice Temperature at Upper Surface', 'surface mixed layer depth', 'bottom mixed layer depth']

#level = -1 #is surface and 0 is bottom

date = datetime.date.today().strftime("%d/%m/%Y")
attrs_global = {'creation_date': date,
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'institution':'Catalan Institute for Water Research (ICRA); Vrije Universiteit Brussel (VUB)',
                        'comment' : 'Data prepared for '+ args.phase +' by DMB, data revised by RM; for validation of the data, visit: https://github.com/icra/GOTM_ISIMIP3_validation'}

final_year = args.finalyear
firstyear = args.firstyear

#nc_individual_path = "/scratch/brussel/vo/000/bvo00012/vsc10623/ISIMIP3a/Output_GOTM/obsclim/historical/GSWP3-W5E5"
nc_individual_path = args.basedir + args.phase + '/Output_GOTM/' + args.datatype + '/' + args.climforcing + '/' + args.model + '/'
os.chdir(nc_individual_path)

#args.outdir = ''
outpath = Path(args.outdir + args.phase + '/OutputData/lakes_global/' + args.datatype + '/' + args.climforcing + '/')
outpath_str= args.outdir + args.phase + '/OutputData/lakes_global/'+ args.datatype + '/' + args.climforcing + '/'

# isimip resolution
resolution = 0.5
coord_lat = [[np.arange(-89.75, 90,0.5)], [np.arange(1,361)]]
coord_lon = [[np.arange(-179.75, 180,0.5)], [np.arange(1,721)]]
lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

lon_da = xr.DataArray(lons, 
                      coords = {'lon':lons}, 
                      dims='lon', 
                      attrs={'standard_name':'longitude', 'long_name':'Longitude', 'units':'degrees_east', 'axis':'X'})

lat_da = xr.DataArray(lats,
                      coords = {'lat':lats}, 
                      dims='lat', 
                      attrs={'standard_name':'latitude', 'long_name':'Latitude','units':'degrees_north', 'axis':'Y'})

#climate data path, to extract decades
climate_path='/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/' + args.phase +'/InputData/climate/atmosphere/' + args.datatype + '/global/daily/' + args.climforcing + '/' + args.model + '/'
def get_periods(path):
    files = sorted(fnmatch.filter(os.listdir(path), '*_tas_*'))
    periods = [file.split('daily_')[1].split('.nc')[0] for file in files]
    return periods

for period in get_periods(climate_path):
    print('comienza' + period)
    first_year = period.split('_')[0]
    if int(first_year) > 1991:
        continue
    last_year = period.split('_')[1]
    start_time = first_year + '-01-01'#insertar period[0] aqui
    end_time = last_year + '-12-31'##insertar period[1] aqui
    time_decade = np.arange(np.datetime64(start_time), np.datetime64(end_time)+1)
    len_time = len(time_decade) 

    #create empty values
    values_surftemp = np.full((len_time,360,720), 1.e+20)
    values_bottemp = np.full((len_time,360,720), 1.e+20)
    values_latentheatf = np.full((len_time,360,720), 1.e+20)
    values_sensheatf = np.full((len_time,360,720), 1.e+20)
    values_icethick = np.full((len_time,360,720), 1.e+20)
    values_icetemp = np.full((len_time,360,720), 1.e+20)
    #values_thermodepth = np.full((len_time,360,720), 1.e+20)
    #values_strat = np.full((len_time,360,720), 1.e+20)

    #time_da = xr.DataArray(time_decade,
    #                       coords = {'time':time_decade}, 
    #                       dims='time', 
    #                       attrs={'units':'days since 1661-01-01','calendar':'proleptic_gregorian', 'standard_name':'time', 'long_name':'Time', 'axis':'T'})

    for index, lake in lakes_df.iterrows():
        print('lago número' + str(index+1))
        lake_file = args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_gotm_' + str(index+1) + '_daily_' + firstyear + '_' + final_year
        #os.system('rm ' + lake_file + 'xr.nc')
        #os.system('ncrename -v z,z_coord '+ lake_file + '.nc ' + lake_file + 'xr.nc')
        lake_ds = xr.load_dataset(nc_individual_path + '/' + lake_file + 'xr.nc')
        time_np = lake_ds['time'].values
        z_np = lake_ds['z'].values
        start_decade = np.where(time_np==np.datetime64(start_time))
        end_decade = np.where(time_np==np.datetime64(end_time))
        lat_pos_values = coord_lat[1][0][np.where(coord_lat[0][0]==lake[1])]
        lon_pos_values = coord_lon[1][0][np.where(coord_lon[0][0]==lake[2])]
        
        for variable_name in variables:
            if variable_name=='temp':
                
                #surface temperature 
                lake_np = lake_ds.sel(z=z_np[-1])[variable_name].values
                lake_np_decade = lake_np[int(start_decade[0]):int(end_decade[0])+1]
                lake_np_decade = lake_np_decade.squeeze(axis=1) #remove one dimension=1 to use with values
                if np.max(lake_np_decade)>100:
                    print("OJO valor mayor a 100º en lago" + str(index+1))
                    #break
                else:
                    values_surftemp[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade + 273.15
                
                #bottom temperature
                lake_np = lake_ds.sel(z=z_np[0])[variable_name].values #set level as 0
                lake_np_decade = lake_np[int(start_decade[0]):int(end_decade[0])+1]
                lake_np_decade = lake_np_decade.squeeze(axis=1) #remove one dimension=1 to use with values
                if np.max(lake_np_decade)>100:
                    print("OJO valor mayor a 100º en lago" + str(index+1))
                    #break
                else:
                    values_bottemp[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade + 273.15

                #thermocline depth
                #lake_np = lake_ds[variable_name].values
                #lake_np = lake_np[:,:,0,0]
                #print(lake_np.shape)
                #print(z_np.shape)
                #thermoD, thermoInd = pylake.thermocline(lake_np, z_np, time_np)
                #print(thermoD)
                #print(thermoD.shape)
                #print(thermoInd)
                #strat=thermoD
                #print(strat)
                #strat[np.where(strat==np.nan)]=0
                #print(strat)
            

            else:
                lake_np = lake_ds.sel(z=z_np[-1])[variable_name].values
                lake_np_decade = lake_np[int(start_decade[0]):int(end_decade[0])+1]
                lake_np_decade = lake_np_decade.squeeze(axis=1) #remove one dimension=1 to use with values
                if variable_name=='qe':
                    values_latentheatf[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade * (-1) #positive is upward
                if variable_name=='qh':
                    values_sensheatf[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade * (-1) #positive i upward
                if variable_name=='Hice':
                    values_icethick[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade
                if variable_name=='Tice_surface':
                    values_icetemp[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade + 273.15
                #if variable_name=='mld_surf':
                #    values_mldsurf[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade
                #if variable_name=='mld_bott':
                #    values_mldbott[:, lat_pos_values-1, lon_pos_values-1] = lake_np_decade
        #if index==11:
        #    break
    #break
    #surface temperature
    #if variable_name=='temp':
    attrs_variable = {'standard_name': variables_save[0], 'long_name': longnames[0], 'units': units[0], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_surftemp, 
                             coords = {'time':time_decade,'lat':lats,'lon':lons},
                             dims=('time','lat','lon'),
                             attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                               'lon' : lon_da,
                               variables_save[0] : values_da}, 
                               attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00' #ojo esto funciona, si quiero cambiar algo
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[0] + '_global_daily_'+ first_year + '_' + last_year + '.nc'
                            
    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[0]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

    #bottom temperature
    attrs_variable = {'standard_name': variables_save[1], 'long_name': longnames[1], 'units': units[1], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_bottemp, 
                             coords = {'time':time_decade,'lat':lats,'lon':lons},
                             dims=('time','lat','lon'),
                             attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                                   'lon' : lon_da,
                                   variables_save[1] : values_da}, 
                                   attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00' #ojo esto funciona, si quiero cambiar algo
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[1] + '_global_daily_'+ first_year + '_' + last_year + '.nc'
                                
    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[1]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

   #latent heat
    #if variable_name=='qe':     
    attrs_variable = {'standard_name': variables_save[2], 'long_name': longnames[2], 'units': units[2], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_latentheatf, 
                             coords = {'time':time_decade,'lat':lats,'lon':lons},
                             dims=('time','lat','lon'),
                             attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                               'lon' : lon_da,
                                variables_save[2] : values_da}, 
                                attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[2] + '_global_daily_'+ first_year + '_' + last_year + '.nc'
                                
    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[2]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

    #sensible heat
    #if variable_name=='qh':     
    attrs_variable = {'standard_name': variables_save[3], 'long_name': longnames[3], 'units': units[3], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_sensheatf, 
                                coords = {'time':time_decade,'lat':lats,'lon':lons},
                                dims=('time','lat','lon'),
                                attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                                   'lon' : lon_da,
                                   variables_save[3] : values_da}, 
                                   attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[3] + '_global_daily_'+ first_year + '_' + last_year + '.nc'

    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[3]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

    #ice thickness
    #if variable_name=='Hice':     
    attrs_variable = {'standard_name': variables_save[4], 'long_name': longnames[4], 'units': units[4], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_icethick, 
                             coords = {'time':time_decade,'lat':lats,'lon':lons},
                             dims=('time','lat','lon'),
                             attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                                   'lon' : lon_da,
                                   variables_save[4] : values_da}, 
                                   attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[4] + '_global_daily_'+ first_year + '_' + last_year + '.nc'

    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[4]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

    #ice temperature
    #if variable_name=='Tice_surface':     
    attrs_variable = {'standard_name': variables_save[5], 'long_name': longnames[5], 'units': units[5], '_FillValue' : 1e20, 'missing_value':1e20}

    values_da = xr.DataArray(values_icetemp, 
                             coords = {'time':time_decade,'lat':lats,'lon':lons},
                             dims=('time','lat','lon'),
                             attrs = attrs_variable)

    ds = xr.Dataset(data_vars={'lat' : lat_da,
                                   'lon' : lon_da,
                                   variables_save[5] : values_da}, 
                                   attrs=attrs_global)

    ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[5] + '_global_daily_'+ first_year + '_' + last_year + '.nc'
                                
    ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save[5]:{'dtype':'float32'}})
    ds.close()

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)
    
    #surface mixed layer depth
    #if variable_name=='mld_surf':     
    #attrs_variable = {'standard_name': variables_save[6], 'long_name': longnames[6], 'units': units[6], '_FillValue' : 1e20, 'missing_value':1e20}

    #values_da = xr.DataArray(values_mldsurf, 
    #                         coords = {'time':time_decade,'lat':lats,'lon':lons},
    #                         dims=('time','lat','lon'),
    #                         attrs = attrs_variable)

    #ds = xr.Dataset(data_vars={'lat' : lat_da,
    #                               'lon' : lon_da,
    #                               variables_save[6] : values_da}, 
    #                               attrs=attrs_global)

    #ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    #output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[6] + '_global_daily_'+ first_year + '_' + last_year + '.nc'
                                
    #ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w')
    #ds.close()

    #os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    #os.system('rm '+ outpath_str + '/' + output_netcdf)
    #os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    #os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    #os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, ' +  outpath_str + '/' +  output_netcdf)

    #bottom mixed layer depth
    #if variable_name=='mld_bott':     
    #attrs_variable = {'standard_name': variables_save[7], 'long_name': longnames[7], 'units': units[7], '_FillValue' : 1e20, 'missing_value':1e20}

    #values_da = xr.DataArray(values_mldbott, 
    #                         coords = {'time':time_decade,'lat':lats,'lon':lons},
    #                         dims=('time','lat','lon'),
    #                         attrs = attrs_variable)

    #ds = xr.Dataset(data_vars={'lat' : lat_da,
    #                               'lon' : lon_da,
    #                               variables_save[7] : values_da}, 
    #                               attrs=attrs_global)

    #ds.time.encoding['units'] = 'days since 1901-01-01 00:00:00'
    
    #output_netcdf = 'gotm_'+ args.model.lower() + '_' + args.datatype + '_histsoc_default_'+ variables_save[7] + '_global_daily_'+ first_year + '_' + last_year + '.nc'

    #ds.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w')
    #ds.close()

    #os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    #os.system('rm '+ outpath_str + '/' + output_netcdf)
    #os.system("ncap2 -h -O -s 'time=double(time)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    #os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    #os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, ' +  outpath_str + '/' +  output_netcdf)

    #break
