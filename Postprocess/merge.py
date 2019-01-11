#!/usr/bin/env python

import os, six

# Set to just 'gdal_merge' if it works for you
gdal_merge = '/opt/local/bin/gdal_merge.py-3.7'
if not os.path.isfile(gdal_merge):
    gdal_merge = 'gdal_merge.py'

# Options to be passed to gdal_merge
merged_filename = 'MergedHypsoPeaks'
# Choose -v for very verbose output
args = ' -ot Float32 -n 32768 '
# args = ' -ot Float32 -v -n -32768 '

# Get a list of all the files in the current directory
files = os.listdir('.')
# Assume we want to process all the .tif files beginning with 'Hypso'
tiffs = [file for file in files 
         if (os.path.isfile(file) and file[:5]=='Hypso' and file[-4:]=='.tif')]
n_tiffs = len(tiffs)
print('Number of GeoTIFFs to merge: {}'.format(n_tiffs))
    
# Call gdal_merge on this giant list
if n_tiffs<1000:
#     os.system(gdal_merge+' -o '+merged_filename+'.tif'+' '+' '.join(tiffs))
    os.system('gdalwarp '+merged_filename+'.tif'+' '+' '.join(tiffs))
else:
    n_batches = n_tiffs//1000
    tmp_filename_list = []
    for batch in range(0,n_batches+1):
        from_tiff = batch*1000
        to_tiff   = min((batch+1)*1000-1,n_tiffs)
        tmp_filename = merged_filename+'_'+str(batch)+'.tif'
        print('Merging files {}-{} to "{}"'.format(from_tiff,to_tiff, tmp_filename))
        sys_call = gdal_merge+' -o '+tmp_filename+args +' '.join(tiffs[from_tiff:to_tiff+1])
        os.system(sys_call)
        tmp_filename_list += [tmp_filename]
    print('Merging temporary files {} to "{}"'.format(tmp_filename_list, merged_filename+'.tif'))
    sys_call = gdal_merge+' -o '+merged_filename+'.tif'+args+' '.join(tmp_filename_list)
    print('System call: "{}"'.format(sys_call))
    os.system(sys_call)
    print('Removing temporary files {}'.format(tmp_filename_list))
    [os.system('rm '+tmp_filename) for tmp_filename in tmp_filename_list]