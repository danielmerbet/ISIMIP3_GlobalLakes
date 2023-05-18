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

path = args.basedir + '/' + args.phase + '/OutputGOTM/' \
       + args.datatype + '/' + args.climforcing + '/' + args.model

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

    # get file names for current period for air temperature to be used as template
    tas_file = [file for file in fnmatch.filter(os.listdir(path), '*_tas_*') if period in file][0]

    # open template
    print("loading datasets")
    template_ds = xr.load_dataset(path + '/' + tas_file)

    outfile = 'gotm' + args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_' + \
        '_daily_' + period[0] + '_' + period[1] + '.nc'
    
    #define decade to print...
    range(period[0], period[1], by=1)


    # iterate over lakes
    for index, lake in lakes_df.iterrows():
        print('{0} (lat:{1} lon:{2})'.format(*lake))
        
        output_gotm = args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_' + \
            lake[0].replace(' ', '-').lower() + '_daily_' + first_year + '_' + last_year + '.nc'


        #open lake data
        lake_temporal = xr.load_dataset(path + '/' + output_gotm)
        
        # replace in template
        #print('   read tas ...')
        template_ds.loc(dict(lat=lake[1], lon=lake[2], , time=)) = lake_temporal.sel(time=, var=)
    
    #print final NetCDF
    template_ds.to_netcdf(path=outfile)
        
    break

print("finished")
print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
