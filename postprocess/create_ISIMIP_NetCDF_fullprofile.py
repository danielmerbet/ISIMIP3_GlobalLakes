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


parser = argparse.ArgumentParser(description='Create final output for ISIMIP3.')

parser.add_argument('-p', '--phase', dest='phase', required=True,
                    default='ISIMIP3a',
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

reftime = 'days since 1601-01-01 00:00:00'

variables='temp'
variables_save='watertemp'
units = 'K'
longnames = 'Temperature of Lake Water'

n_level = 50

#level = -1 #is surface and 0 is bottom

attrs_variable = {'long_name': longnames,'units': units, '_FillValue' : 1e20, 'missing_value':1e20}
attrs_variable_depth = {'standard_name':'depth_below_surface','long_name': 'Depth of Vertical Layer Center Below Surface','units':'m', 'positive' : 'down', 'axis':'Z',  '_FillValue' : 1e20, 'missing_value':1e20}
date = datetime.date.today().strftime("%d/%m/%Y")
attrs_global = {'creation_date': date,
                        'contact' : 'Daniel Mercado-Bettín - ICRA (dmercado@icra.cat); Rafael Marcé - ICRA (rmarce@icra.cat)',
                        'institution':'Catalan Institute for Water Research (ICRA); Vrije Universiteit Brussel (VUB)',
                        'comment' : 'Data prepared for '+ args.phase +' by DMB, data revised by RM; for validation of the data, visit: https://github.com/icra/GOTM_ISIMIP3_validation' }

final_year =args.finalyear
firstyear = args.firstyear

#nc_individual_path = "/scratch/brussel/vo/000/bvo00012/vsc10623/ISIMIP3a/Output_GOTM/obsclim/historical/GSWP3-W5E5"
nc_individual_path = args.basedir + args.phase + '/Output_GOTM/' + args.datatype + '/' + args.climforcing + '/' + args.model + '/'
os.chdir(nc_individual_path)

#args.outdir = ''
outpath = Path(args.outdir + args.phase + '/OutputData/lakes_global/' + args.datatype + '/' + args.climforcing + '/')
outpath_str= args.outdir + args.phase + '/OutputData/lakes_global/' + args.datatype + '/' + args.climforcing + '/'

# isimip resolution hardcoded
resolution = 0.5
coord_lat = [[np.arange(-89.75, 90,0.5)], [np.arange(1,361)]]
coord_lon = [[np.arange(-179.75, 180,0.5)], [np.arange(1,721)]]

#to ssave netcdf
lons = np.arange(-180+resolution/2,180+resolution/2,resolution)
lats = np.arange(-90+resolution/2,90+resolution/2,resolution)
#levlak = np.arange(1,11)
levlak = np.arange(1,51)

levlak_da = xr.DataArray(levlak,
                      coords = {'levlak':levlak}, 
                      dims='levlak', 
                      attrs={'standard_name': 'water_layer', 'long_name': 'Vertical Water Layer Index','units':'-', 'axis':'Z', 'positive':'down'})

lon_da = xr.DataArray(lons, 
                      coords = {'lon':lons}, 
                      dims='lon', 
                      attrs={'standard_name': 'longitude', 'long_name': 'Longitude', 'units':'degrees_east', 'axis':'X'})

lat_da = xr.DataArray(lats,
                      coords = {'lat':lats}, 
                      dims='lat', 
                      attrs={'standard_name':'latitude', 'long_name':'Latitude','units':'degrees_north', 'axis':'Y'})


#climate data path, to extract decades
climate_path='/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/' + args.phase +'/InputData/climate/atmosphere/bias-adjusted/global/daily/' + args.climforcing + '/' + args.model + '/'
def get_periods(path):
    files = sorted(fnmatch.filter(os.listdir(path), '*_tas_*'))
    periods = [file.split('daily_')[1].split('.nc')[0] for file in files]
    return periods

if args.climforcing=='historical':
    total_periods = ['1850_1850']
    for i in range(1850, 2010, 10):
        total_periods.append(str(i+1) + '_' + str(i+10))
    total_periods.append('2011_2014')

if args.climforcing=='ssp126' or args.climforcing=='ssp370' or args.climforcing=='ssp585':
    total_periods = ['2015_2020']
    for i in range(2020, 2100, 10):
        total_periods.append(str(i+1) + '_' + str(i+10))

#for period in get_periods(climate_path):
for period in total_periods:
    print('comienza' + period)
    # delete if file exists
    #if os.path.isfile(filename_netcdf):
       # os.system('rm '+filename_netcdf)
    #inputs
    first_year = period.split('_')[0]
    #if int(first_year) < 1951:
    #    continue
    last_year = period.split('_')[1]
    start_time = first_year + '-01-01'#insertar period[0] aqui
    end_time = last_year + '-12-31'##insertar period[1] aqui
    time_decade = np.arange(np.datetime64(start_time), np.datetime64(end_time)+1)
    len_time = len(time_decade) 


    #create empty values for watertemp (time, levlak, lat, lon)
    values = np.full((len_time,n_level,360,720), 1.e+20)
    
    #create values for depth (levlak, lat, lon)
    values_depth = np.full((n_level,360,720), 1.e+20)
    #values.shape

    #time_da = xr.DataArray(time_decade,
    #                       coords = {'time':time_decade},
    #                       dims='time',
    #                       attrs={'units':'days since 1661-01-01', 'calendar':'proleptic_gregorian', 'standard_name':'time', 'long_name':'Time', 'axis':'T'})
        
    for index, lake in lakes_df.iterrows():
        print('lago número' + str(index+1))
        lake_file = args.model.lower() + '_' + args.climforcing + '_bias-adjusted_gotm_' + str(index+1) + '_daily_' + firstyear + '_' + final_year
        #os.system('rm ' + lake_file + 'xr.nc')
        #os.system('ncrename -v z,z_coord '+ lake_file + '.nc ' + lake_file + 'xr.nc')
        lake_ds = xr.load_dataset(nc_individual_path + '/' + lake_file + 'xr.nc')
        time_np = lake_ds['time'].values
        z_np = lake_ds['z'].values
        z_np = np.reshape(z_np, (len(z_np), 1))
        start_decade = np.where(time_np==np.datetime64(start_time))
        end_decade = np.where(time_np==np.datetime64(end_time))
        lat_pos_values = coord_lat[1][0][np.where(coord_lat[0][0]==lake[1])]
        lon_pos_values = coord_lon[1][0][np.where(coord_lon[0][0]==lake[2])]

        lake_np = lake_ds[variables].values
        lake_np_decade = lake_np[int(start_decade[0]):int(end_decade[0])+1]
        lake_np_decade = lake_np_decade.squeeze(axis=3) #remove dimension number 4 to be able to use with values
        if np.max(lake_np_decade)>100:
            print("OJO valor mayor a 100º en lago" + str(index+1))
            #break
        else:
            sel_levels = np.round(np.linspace(1,lake_np_decade.shape[1],n_level)).astype(int)-1
            values[:, :,lat_pos_values-1, lon_pos_values-1] = lake_np_decade[:,sel_levels,:] + 273.15
            print(values) 
            values_depth[:,lat_pos_values-1, lon_pos_values-1] = z_np[sel_levels]
        if index==11:
            break


    values_da = xr.DataArray(values, 
                            coords = {'time':time_decade,'levlak':levlak, 'lon':lons,'lat':lats},
                            dims=('time','levlak','lat','lon'),
                            attrs = attrs_variable)

    values_depth_da = xr.DataArray(values_depth,
                            coords = {'levlak':levlak,'lon':lons,'lat':lats},
                            dims=('levlak','lat','lon'),
                            attrs = attrs_variable_depth)

    ds = xr.Dataset(data_vars={'levlak': levlak_da,   
                                'lat' : lat_da,
                                'lon' : lon_da,
                                variables_save : values_da}, 
                                attrs=attrs_global)

    ds.time.encoding['units'] = reftime

    ds_depth = xr.Dataset(data_vars={'levlak':levlak_da,
                                'lat' : lat_da,
                                'lon' : lon_da,
                                'depth' : values_depth_da},
                                 attrs=attrs_global)


    ds_merged = xr.merge([ds, ds_depth])

    output_netcdf = 'gotm_'+ args.model.lower() + '_'  + args.datatype + '_' + args.climforcing + '_histsoc_default_'+ variables_save + '_global_daily_'+ first_year + '_' + last_year + '.nc'
    #'gotm' + '_' + first_year + '_' + last_year + '.nc'
                                
    ds_merged.to_netcdf(outpath / output_netcdf, format='NETCDF4_CLASSIC', mode='w', encoding={variables_save:{'dtype':'float32'}})

    os.system('cdo -f nc4c -z zip_5 -O --history -setmissval,1e+20 -invertlat ' + outpath_str + '/' + output_netcdf + ' ' + outpath_str + '/' +  output_netcdf + 'chito')
    os.system('rm '+ outpath_str + '/' + output_netcdf)
    os.system("ncap2 -h -O -s 'time=double(time)' -s 'levlak=double(levlak)' " + outpath_str + '/' +  output_netcdf + 'chito ' +  outpath_str + '/' +  output_netcdf)
    os.system('rm '+ outpath_str + '/' + output_netcdf + 'chito')
    os.system('ncatted -h -O -a Conventions,global,d,, -a CDI,global,d,, -a CDO,global,d,, -a long_name,time,o,c,Time ' +  outpath_str + '/' +  output_netcdf)

    #break
