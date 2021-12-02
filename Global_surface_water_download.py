#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 20:46:27 2021

@author: aparrar1
Code for querying and downloading data from the Global surface water dynamics 1999-2020 dataset.
https://glad.umd.edu/dataset/global-surface-water-dynamics
The dataset is stored in a Google Cloud platform bucket storage. 
Following the filtering data process of the GEDI data which includes the removal of points with a Landsat water persistence value of less than 80%
this code was created to include information from the Global surface water dynamics into ICESat data, to follow the same filtering process. 
According to the GEDI L2B product documentation, the Landsat water persistence flag is defined as "The percent UMD GLAD Landsat observations with classified surface water between 2018 and 2019."
Following this description the annual water percent maps for 2018 and 2019 are used for calculating the percent surface water.
The Google Cloud platform bucket storage can be accessed at:
https://console.cloud.google.com/storage/browser/earthenginepartners-hansen/water/
This code was based on publicly available code and resources:
https://global-surface-water.appspot.com/download
https://storage.googleapis.com/global-surface-water/downloads_ancillary/downloadWaterData_PythonV3_2020.zip

"""
import os
import sys
import urllib.request
import rasterio
from rasterio.merge import merge
import numpy as np
import geopandas as gpd
import pandas as pd
import itertools
sys.path.append('/Users/asparra/Documents/Global_deltas/code')
import Satelite_LIDAR_download_subset as dpf


#Get the lat and lon string values needed to query the Global surface water dynamics dataset from the ROI shapefile
def lon_lat_selection(roi_shp):
    #make sure the roi_shp file exists
    if not os.path.isfile(roi_shp):    
        print('error: ROI shpefile provided does not exist or was not found')
        sys.exit(2)
    #get the lat and lon values of the ROI bounding box
    BBOX=dpf.ROI_to_BBOX(roi_shp) 
    BBOX=BBOX.total_bounds
    lon=[BBOX[0],BBOX[2]] 
    lat=[BBOX[1],BBOX[3]]
    #the Global surface water dynamics dataset is divided in a grid of 10 degrees cells. 
    #Check in which part of the grid the lat and lon fall and change the negative sign to S and W values
    longs_W = [w for w in range(180,0,-10)]
    longs_E = [e for e in range(0,180,10)]
    lats_S = [s for s in range(50,0,-10)]
    lats_N = [n for n in range(0,90,10)]
    lon_values=[]
    for item in lon:
        if item<0:
            letter='W'
            for i in enumerate(longs_W):
                if abs(item)<=i[1] and abs(item)>i[1]-10:
                    filtered_lon=i[1]
                    if filtered_lon>=100:
                        lon_values.append(str(filtered_lon)+letter)
                    else:
                        lon_values.append('0'+str(filtered_lon)+letter)
        else:
            letter='E'
            for i in enumerate(longs_E):
                if item>=i[1] and item<i[1]+10:
                    filtered_lon=i[1]
                    if filtered_lon>=100:
                        lon_values.append(str(filtered_lon)+letter)
                    elif filtered_lon==0:
                        lon_values.append('000'+letter)
                    else:
                        lon_values.append('0'+str(filtered_lon)+letter)
    lat_values=[]
    for item in lat:
        if item<=-10:
            letter='S'
            for i in enumerate(lats_S):
                if abs(item)<=i[1] and abs(item)>i[1]-10:
                    filtered_lat=i[1]
                    lat_values.append(str(filtered_lat)+letter)
        elif item>-10 and item<=0:
            letter='N'
            lat_values.append('00'+letter)        
        else:
            letter='N'
            for i in enumerate(lats_N):
                if item>i[1] and item<=i[1]+10:
                    filtered_lat=i[1]+10
                    lat_values.append(str(filtered_lat)+letter)  
    final_values = list(itertools.product(set(lat_values), set(lon_values)))
    return(final_values)   
    
def water_surface_download(lat, lon, output_directory): #lat and lon values should be strings with 00N and 000E or 00S and 000W format.
    #make sure the output directory ends with a /
    if (output_directory[-1:]!="/"):
        output_directory = output_directory + "/"
    #make sure the output_directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    subdirectory=lat+'_'+lon+'/'
    filename1 = '2018_percent.tif'
    output_file1= '2018_percent'+'_'+ lat+'_'+lon+'.tif'
    filename2 = '2019_percent.tif'
    output_file2= '2019_percent'+'_'+ lat+'_'+lon+'.tif'
    if os.path.exists(output_directory +  output_file1):
        print(output_directory + filename1 + " already exists - skipping")
    else:
        try:
            url = "http://storage.googleapis.com/earthenginepartners-hansen/water/" + subdirectory + filename1
            code = urllib.request.urlopen(url).getcode()
            if (code != 404):
                print("Downloading " + url)
                urllib.request.urlretrieve(url, output_directory + output_file1)
            else:
                print(url + " not found")
        except ValueError as val_error:
            print(val_error)
        except Exception as error:
            print(error)
    if os.path.exists(output_directory +  output_file2):
        print(output_directory + filename2 + " already exists - skipping")
    else:
        try:
            url = "http://storage.googleapis.com/earthenginepartners-hansen/water/" + subdirectory + filename2
            code = urllib.request.urlopen(url).getcode()
            if (code != 404):
                print("Downloading " + url)
                urllib.request.urlretrieve(url, output_directory + output_file2)
            else:
                print(url + " not found")
        except ValueError as val_error:
            print(val_error)
        except Exception as error:
            print(error)
            
def mosaic_construction(input_directory, filename):
#Find all the 2019 and 2018 files that correspond to a ROI, mosaic them and then create a final average raster  
    #make sure the input directory ends with a /
    if (input_directory[-1:]!="/"):
        input_directory = input_directory + "/"
    #make sure the input_directory exists
    if not os.path.exists(input_directory):
        print('The input directory does not exists')
        sys.exit(2)
    else:
        #get the available files in the directory
        available_files=[os.path.join(input_directory, file) for file in os.listdir(input_directory)]
        files_2018=[file for file in available_files if '2018_percent' in file]
        files_2019=[file for file in available_files if '2019_percent' in file]
        #Check that the length of the files matches
        if len(files_2018)!=len(files_2019):
            print('The number of files for 2018 and 2019 are different!')
            sys.exit(2)  
        files_list_2018 = []
        for file in files_2018:
            input_file = rasterio.open(file)
            files_list_2018.append(input_file)
        files_list_2019 = []
        for file in files_2019:
            input_file = rasterio.open(file)
            files_list_2019.append(input_file)
        mosaic2018, out_trans_2018 = merge(files_list_2018)
        mosaic2019, out_trans_2019 = merge(files_list_2019)
        array_list = [mosaic2018, mosaic2019]
        # Perform averaging
        array_out = np.nanmean(array_list, axis=0)
        output_file=input_directory+filename
        new_dataset = rasterio.open(output_file,'w',driver='GTiff',
                            height=mosaic2018.shape[1],width=mosaic2018.shape[2],count=1,dtype=mosaic2018.dtype,
                            crs=files_list_2018[0].crs, transform=out_trans_2018)
        new_dataset.write(array_out)
        new_dataset.close()

def extract_water_persistence(shp_file, raster_file): #future fixes!! add loops for files that are too large
    # Read points from shapefile
    if '.shp' in shp_file or '.json' in shp_file:
        lidar_points = gpd.read_file(shp_file)
        coords = [(x,y) for x, y in zip(lidar_points.geometry.x, lidar_points.geometry.y)]
    else:
        lidar_points = pd.read_csv(shp_file)
        coords = [(x,y) for x, y in zip(lidar_points.Lon, lidar_points.Lat)]
    # Open the Landsat water persistence raster 
    water_persistence = rasterio.open(raster_file)  
    # Sample the raster at every point location and store values in a DataFrame
    out_data = pd.DataFrame() 
    out_data['landsat_WP'] = [x[0] for x in water_persistence.sample(coords)]
    #Add the resulting column into the input shapefile following the GEDI name convention
    lidar_points = lidar_points.assign(landsat_WP=out_data)
    #lidar_points.to_file(shp_file)
    return(lidar_points)
