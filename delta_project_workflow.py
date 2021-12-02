#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 07:56:36 2021

@author: aparrar1
Code for GEDI, ICESat1 GLAS and ICESat2 ATLAS data processing, for the Global Delta project
This code defines functions for creating a workflow for querying, downloading, and subsetting GEDI and ICESat data
This code was based on publicly available code and resources.

The complete data process is done in parallel per file.
Some ROI have very large areas, for these cases downloading all the files that intercept the bounding box results in large demands of CPU space. 
To avoid this, the area of the ROI is first evaluated, if the area is larger than a threshold, the ROI is divided in equal parts. 
An initial search of all the files that intercept each subsection of the ROI is performed and repeated files are deleted. 
Once the final list of available files is complete the workflow is done to each file. 
The file is downloaded, the data extraction is ran, and if there are no point intercepts the file is deleted. 

#A netrc_file needs to be created to be able to download the files. 
#Follow the instructions on the Earth Data login setup code EarthdataLoginSetup.py created by Cole Krehbiel to construct the netrc_file:
#https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_python/browse/EarthdataLoginSetup.py

"""

import os
import sys
import re
import geopandas as gp
import pandas as pd
import numpy as np

sys.path.append('/Users/asparra/Documents/Global_deltas/code')#Change the directory to the folder with the Satelite_LIDAR_download_subset.py file
import Satelite_LIDAR_download_subset as dpf
  
#Set the high level directories, the number of workers, and the file extension (shp, json, or csv) 
output_dir='/Users/asparra/Documents/Global_deltas/data/'
temporary_dir='/Volumes/Fortress_L3/'
NUM_WORKERS = 5
extension='csv'

#get the ROI shapefiles for each delta and loop through each one 
ROI_shapefiles=['/Users/asparra/Documents/Global_deltas/data/Shapefiles/DeltaSHPs/global_map_VistulaBuff5km.shp']

if __name__ == '__main__':
    for i in ROI_shapefiles:
        delta_name=re.split(".shp", i.split("/")[-1])[0]
        #Modify this part of the code depending on the name of your file.
        delta_name=re.sub('global_map_|Buff5km', '', delta_name)
        output_directory=output_dir + delta_name
        #Modify the file_types depending on the dataset you want to download
        file_types = {"GEDI": ['GEDI02_B', '002', 'gedib_query',['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']],
                      "GLAS": ['GLAH14', '034', 'glas_query','Data_40HZ'],
                      "ATLAS03": ['ATL03', '004', 'atlas03_query',['gt1l', 'gt1r', 'gt2l','gt2r','gt3l','gt3r']], 
                      "ATLAS08": ['ATL08', '004', 'atlas08_query',['gt1l', 'gt1r', 'gt2l','gt2r','gt3l','gt3r']]}
        #Check the size of the ROI. If the area is bigger than 2000 km2 divide the BBOX in smaller parts to avoid querying files for empty spaces of the BBOX
        shapefile = gp.read_file(i)
        shapefile = shapefile['geometry'].to_crs({'proj':'cea'}) 
        shp_area=shapefile.area / 10**6
        if shp_area[0]<2000:
            query_type=True
        else:
            query_type=False
        for n in file_types:
            data_name=n
            short_name=file_types[n][0]
            version=file_types[n][1]
            query_name=file_types[n][2]
            data_group=file_types[n][3]
            output_directory2=output_directory + '/' + data_name
            temporary_directory=temporary_dir + data_name
            output_file=output_directory2+'/'+delta_name+'_complete_' + data_name +'.'+extension
            #Check if the complete shapefile for the delta and for each data type already exists. If it does, skip to the next element.
            if os.path.isfile(output_file) or os.path.isfile(re.sub('.'+extension, '_0.'+extension, output_file)):
                continue
            else:
                #Check if the data query file already exists and get the list of overlaping files from that file
                query_file=output_directory2+ '/' + query_name
                if os.path.isfile(query_file):
                   query_list = open(query_file, "r") 
                   files = query_list.readlines()
                   available_files = [re.sub('\n', '', file) for file in files] 
                   if re.search("No available files for this ROI", available_files[0]):
                       continue
                else:   
                    if query_type: 
                        available_files=dpf.cmr_query_wrapper(ROI_file=i, output_directory=output_directory2, short_name=short_name, version=version, file_name=query_name, save_output=False)  
                    else:
                        complete_roi = gp.read_file(i)
                        if complete_roi.crs != 'EPSG:4326':
                            complete_roi  = complete_roi.to_crs("EPSG:4326")
                        bbox = complete_roi.bounds
                        lon=[bbox['maxx'][0],bbox['minx'][0]]  
                        lat=[bbox['maxy'][0],bbox['miny'][0]]
                        if round((lat[0]-lat[1])/0.5)==1:
                            step_lat=2
                        else:
                            step_lat=round((lat[0]-lat[1])/0.5)
                        if round((lon[0]-lon[1])/0.5)==1:
                            step_lon=2
                        else:
                            step_lon=round((lon[0]-lon[1])/0.5)
                        latitudes=np.linspace(lat[0], lat[1], step_lat)
                        longitudes=np.linspace(lon[0], lon[1], step_lon)
                        available_files=[]
                        for x in range(len(latitudes)-1):
                            for z in range(len(longitudes)-1):
                                new_bbox=[latitudes[x], longitudes[z+1], latitudes[x+1],longitudes[z]]
                                ROI=dpf.ROI_to_BBOX(new_bbox)
                                final_ROI = gp.overlay(ROI, complete_roi)
                                if final_ROI.shape[0] != 0:
                                    file_query=dpf.cmr_query_wrapper(ROI_file=new_bbox, output_directory=output_directory2, short_name=short_name, version=version, file_name=query_name, save_output=False) 
                                    available_files.append(file_query)
                        #Remove duplicate files
                        complete_list= [item for sublist in available_files for item in sublist]
                        available_files=list(set(complete_list))        
                #Get the final list of h5 files to download
                available_files=[e for e in available_files if ('.h5' in e or '.H5' in e) and not(".xml" in e or ".MET" in e or ".iso" in e)]
                if len(available_files)==0:
                    if not os.path.exists(output_directory2):
                        os.makedirs(output_directory2)
                    avail_files=open(query_file, "w")
                    avail_files.write('No available files for this ROI')
                    avail_files.close()  
                    continue
                else:
                    #Download and process each file that ovelays the ROI.
                    dpf.parallel_workflow_per_file(available_files, temporary_directory, i, output_directory2, data_group, NUM_WORKERS, extension)
                    #Update the list of available files, to avoid downloading files that do not cover the ROI if the downloading has to be repeated
                    final_files = [item for item in available_files if re.split(".h5|.H5", item.split("/")[-1])[0]+'.'+extension in os.listdir(output_directory2)]
                    #Since some datasets are still producing data, the list of available files for a specific ROI should be updated every once in a while.
                    #check if there are no files that intercept the ROI
                    avail_files=open(query_file, "w")
                    if len(final_files)==0:
                        avail_files.write('No available files for this ROI')
                    else:
                        for element in final_files:
                            avail_files.write(element + "\n")
                    avail_files.close()        
                    #Merge the resulting shapefiles into a single file per delta
                    #Check if there are available files
                    files = os.listdir(output_directory2)
                    paths = [os.path.join(output_directory2, x) for x in files if '.'+extension in x]
                    if len(paths)!=0:
                        if len(paths)<225: #If there are more than 225 files create separate shapefile files per groups of 200 files
                            if extension=='shp'or extension=='json':
                                complete_shapefile = gp.GeoDataFrame(pd.concat([gp.read_file(y) for y in paths]), crs=gp.read_file(paths[0]).crs)
                                complete_shapefile.to_file(output_file)
                            else:
                                complete_file = pd.concat([pd.read_csv(y) for y in paths])
                                complete_file.to_csv(output_file, index=False)
                        else:
                            nums1=list(np.arange(0, len(paths), 201))
                            nums2=list(np.arange(0, len(paths), 200)[1:])
                            nums2.append(len(paths))
                            for number in range(len(nums1)):
                                output_file2=re.sub('.'+extension, ('_'+str(number)+'.'+extension), output_file) 
                                path_list=paths[nums1[number]:nums2[number]]
                                if extension=='shp'or extension=='json':
                                    complete_shapefile = gp.GeoDataFrame(pd.concat([gp.read_file(y) for y in path_list]), crs=gp.read_file(paths[0]).crs)
                                    complete_shapefile.to_file(output_file2)
                                else:
                                    complete_file = pd.concat([pd.read_csv(y) for y in path_list])
                                    complete_file.to_csv(output_file2, index=False)
                                    
 


