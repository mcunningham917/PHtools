#!/usr/bin/env python

import os

# Set to just 'gdal_merge' if it works for you
gdal_merge = 'gdal_merge.py'
# gdal_merge = '/opt/local/bin/gdal_merge.py-3.7'

# Options to be passed to gdal_merge
args = '-o MergedHypsoPeaks.tif'
# args = '-o MergedHypsoPeaks.tif -ot Float32 -v -n -32768'

# Get a list of all the files in the current directory
files = os.listdir('.')
# Assume we want to process all the .tif files beginning with 'Hypso'
tiffs = [file for file in files 
         if (os.path.isfile(file) and file[:5]=='Hypso' and file[-4:]=='.tif')]

# Call gdal_merge on this giant list
os.system(gdal_merge+' '+args+' '+' '.join(tiffs))