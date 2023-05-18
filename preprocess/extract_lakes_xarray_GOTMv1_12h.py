#!/usr/bin/env python3
# run with Python 3.8

import argparse
import os
import fnmatch
import pandas as pd
import xarray as xr
import numpy as np
#from calc_cc_dewp import calc_cc_dewp
from pathlib import Path
import datetime
from math import pi

print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

def calc_cc_dewp(date, airt, relh, swr, lat, lon, elev):
  
  date = pd.to_datetime(date, format='%d/%m/%Y %H:%M:%S')
  date = pd.date_range(start=date[0],end=(date[-1]+datetime.timedelta(days = 1)-datetime.timedelta(hours = 1)),freq='H')
    

  yday = np.array(date.dayofyear)
  
  hour = np.array(date.hour)
  hour[hour==0]=24

  stdmer = range(-180,180, 15)
  Lsm = stdmer[np.argmin([abs(lon-x) for x in list(stdmer)])] # Local standard meridian (degrees)

  Hsc = 1390 # Solar constant (W/m2)
  cd = 0.06 # Dust coefficient
  Rg = 0.045 # Reflectivity of the ground - extended mixed forest


  theta = lat*pi/180 # Latitude in radians

  r = 1 + 0.017 * np.cos((2*pi/365)*(186-yday)) # Relative earth-sun distance

  d = 23.45 * pi/180 * np.cos((2*pi/365)*(172-yday)) # Declination of the sun

  dts = (1/15) * (Lsm-lon) # Fraction of 15-degree increment Llm is east of Lsm
  value = (np.sin(theta)*np.sin(d))
  value = value/(np.cos(theta)*np.cos(d))
  
  value[value > 1] = 1
  value[value < -1] = -1

  tss = (12/pi) * np.arccos(-value) + dts + 12 # Time of sunset
  tsu = -tss + (2 * dts) + 24 # Time of sunrise

  gamma = np.repeat(0, len(tss)) # Correction factor
  dum = np.where(np.logical_and(hour>tsu, hour<tss))
  gamma[dum] = 1

  #Calculate Hb and Htheta
  dum1 = np.where(hour <=12 )
  dum2 = np.where(hour > 12 )
  hb1  = pi/12*(hour-1-dts)
  hb1[dum1] = hb1[dum1]+pi
  hb1[dum2] = hb1[dum2]-pi
  hb  = hb1
  dum3 = np.where(hb1 > 2*pi)
  hb[dum3] = hb[dum3] - 2 * pi
  dum4 = np.where(hb1 < 0)
  hb[dum4] = hb[dum4] + 2 * pi
  #rm(c(dum3, dum4))
  he1  = pi/12*(hour-dts)
  he1[dum1] = he1[dum1]+pi
  he1[dum2] = he1[dum2]-pi
  he  = he1
  dum3 = np.where(he1 > 2*pi)
  he[dum3] = he[dum3] - 2*pi
  dum4 = np.where(he1 < 0)
  he[dum4] = he[dum4] + 2*pi
  #clear dum1 dum2 dum3 dum4

  Ho = Hsc/(r**2)*(np.sin(theta)*np.sin(d)+12/pi*np.cos(theta)*np.cos(d)*(np.sin(he)-np.sin(hb)))*gamma

  # Radiation scattering and absorption #####################################

  w = (he+hb)/2 # Hour angle
  alpha1 = abs(np.sin(theta)*np.sin(d)+np.cos(theta)*np.cos(d)*np.cos(w))
  alpha = np.arctan(alpha1/np.sqrt(1-alpha1**2)) # Solar altitude

  theta_am1 = ((288-0.0065*elev)/288)**5.256
  theta_am2 = np.sin(alpha)+0.15*((alpha*180/pi)+3.855)**(-1.253)
  theta_am = theta_am1/theta_am2 # Optical air mass

  # Dewpoint temperature
  dewt_daily = 243.04*(np.log(relh/100)+((17.625*airt)/(243.04+airt)))/(17.625-np.log(relh/100)-((17.625*airt)/(243.04+airt)))
  dewt_hourly = np.repeat(dewt_daily, 24)  

  Pwc = 0.85*np.exp(0.11+0.0614*dewt_hourly) # Precipitable atmospheric water content

  a2 = np.exp(-(0.465+0.134*Pwc)*(0.179+0.421*np.exp(-0.721*theta_am))*theta_am) # Atmospheric transmission coefficient after scattering and absorption
  a1 = np.exp(-(0.465+0.134*Pwc)*(0.129+0.171*np.exp(-0.88*theta_am))*theta_am)
  at = (a2+0.5*(1-a1-cd))/(1-0.5*Rg*(1-a1-cd)) # attenuation (scattering and absorption)
  #att = mean(at)

  Ho = at*Ho
  #Ho = att*Ho

  dum5 = np.where(Ho<0)
  Ho[dum5] = 1

  
  Ho_daily = np.average(Ho.reshape(-1, 24), axis=1)  
  
  ccsim = np.empty(len(Ho_daily))
  ccsim[:] = np.nan
  for i in range(0,len(Ho_daily)):
    if Ho_daily[i] > swr[i]:
        ccsim[i] = np.sqrt((1 - (swr[i]/Ho_daily[i]))/0.65)
  
  ccsim[ccsim > 1] = 1
  
  #this happens usually in winter in high latitudes, put a mean value of 0.5
  if np.isnan(ccsim[0]):
      ccsim[0] = 0.5
  if np.isnan(ccsim[-1]):
      ccsim[0] = 0.5

  ccsim = pd.DataFrame(ccsim).interpolate()
  ccsim = np.array(ccsim[0])
  
  ccsim[ccsim > 1] = 1
  

  return ccsim, dewt_daily


parser = argparse.ArgumentParser(description='Extract climate data for local lakes sector.')

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
                    default='/p/projects/isimip/isimip',
                    help='base directory')
parser.add_argument('-o', '--out', dest='outdir',
                    default='/p/tmp/buechner/extract_lakes_python',
                    help='base directory to post-processed data')

args = parser.parse_args()
ref_year_isimip = 1901

path = args.basedir + '/' + args.phase + '/InputData/climate/atmosphere/' \
       + args.datatype + '/global/daily/' + '/' + args.climforcing + '/' + args.model

outpath = Path(args.outdir)

lakes_df = pd.read_csv('coord_area_depth.csv')


# return available time periods of our daily data chunked to decades
def get_periods(path):
    files = sorted(fnmatch.filter(os.listdir(path), '*_tas_*'))
    periods = [file.split('daily_')[1].split('.nc')[0] for file in files]
    return periods


first_year = get_periods(path)[0].split(sep='_')[0]
last_year = get_periods(path)[-1].split(sep='_')[1]

for period in get_periods(path):
    print(' Period: ', period)

    # get file names for current period
    tas_file = [file for file in fnmatch.filter(os.listdir(path), '*_tas_*') if period in file][0]
    hurs_file = [file for file in fnmatch.filter(os.listdir(path), '*_hurs_*') if period in file][0]
    pr_file = [file for file in fnmatch.filter(os.listdir(path), '*_pr_*') if period in file][0]
    rsds_file = [file for file in fnmatch.filter(os.listdir(path), '*_rsds_*') if period in file][0]
    rlds_file = [file for file in fnmatch.filter(os.listdir(path), '*_rlds_*') if period in file][0]
    ps_file = [file for file in fnmatch.filter(os.listdir(path), '*_ps_*') if period in file][0]
    sfcwind_file = [file for file in fnmatch.filter(os.listdir(path), '*_sfcwind_*') if period in file][0]

    # open datasets
    print("loading datasets")
    print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    tas_ds = xr.load_dataset(path + '/' + tas_file)
    hurs_ds = xr.load_dataset(path + '/' + hurs_file)
    pr_ds = xr.load_dataset(path + '/' + pr_file)
    rsds_ds = xr.load_dataset(path + '/' + rsds_file)
    rlds_ds = xr.load_dataset(path + '/' + rlds_file)
    ps_ds = xr.load_dataset(path + '/' + ps_file)
    sfcwind_ds = xr.load_dataset(path + '/' + sfcwind_file)
    print("loading finished")
    print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

    # iterate over lakes
    for index, lake in lakes_df.iterrows():
        #print('{0} (lat:{1} lon:{2})'.format(*lake))

        outfile = args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_' + \
            lake[0].replace(' ', '-').lower() + '_daily_' + first_year + '_' + last_year + '.csv'

        # read actual data from NetCDF
        #print('   read tas ...')
        tas_lake = tas_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read hurs ...')
        hurs_lake = hurs_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read pr ...')
        pr_lake = pr_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read rsds ...')
        rsds_lake = rsds_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read rlds ...')
        rlds_lake = rlds_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read ps ...')
        ps_lake = ps_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        #print('   read sfcwind ...')
        sfcwind_lake = sfcwind_ds.sel(lat=lake[1], lon=lake[2], method='nearest')
        
        # write csv files per lake
        #print('   write data ...')
        lake_data = sfcwind_lake.to_dataframe().round(2)
        #lake_data['hour'] = np.repeat('12:00:00', len(lake_data['sfcwind']))
        lake_data = lake_data[['sfcwind']]
        lake_data['sfcwind_e'] = np.repeat('0', len(lake_data['sfcwind']))
        lake_data['ps'] = ps_lake['ps']/100
        lake_data['tas'] = tas_lake['tas']-273.15
        cc, dew = calc_cc_dewp(np.array(hurs_lake['time']), np.array(tas_lake['tas'])-273.15, np.array(hurs_lake['hurs']), np.array(rsds_lake['rsds']),np.array(hurs_lake['lat']),np.array(hurs_lake['lon']),lake[6])
        lake_data['dew'] = dew
        lake_data['cc'] = cc
        lake_data['rsds'] = rsds_lake['rsds']
        lake_data['pr'] = pr_lake['pr']*86400/(1000*24*60*60)
        
        if period.split('_')[0] == first_year:
            lake_data.to_csv(outpath / outfile, header=False, sep="\t")
        else:
            lake_data.to_csv(outpath / outfile, mode='a', header=False, sep="\t")

        # break after first lake for faster debugging
        #break
    #break

print("finished")
print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
