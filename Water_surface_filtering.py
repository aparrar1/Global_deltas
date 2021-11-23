#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 12:04:42 2021

@author: aparrar1
Code for GEDI, ICESat GLAS and ICESat2 ATLAS data processing, for the Delta modelling project
Depending on the type of data (GEDI, ICESat1, or ICESat2), different filtering parametes will be used to remove points over land areas
Points with degraded conditions or high uncertaninty will also be removed. 
The Global surface water dynamics 1999-2020 dataset will be used to filter points that fall on water areas for ICESat data

"""

import os
import sys
import re
import geopandas as gpd
sys.path.append('/Users/asparra/Documents/Global_deltas/code')
import Global_surface_water_download as gsw


"""
Filtering data points to select only information corresponding to water surface elevation
"""

def water_surface_filtering(roi_shp, input_file, elevation_threshold=50, overwrite=False, aditional_columns=None):
#Check that the input file and roi shp exists
    if not os.path.isfile(roi_shp):    
        print('error: ROI shpefile provided does not exist or was not found')
        sys.exit(2)
    if not os.path.isfile(input_file):    
        print('error: input file provided does not exist or was not found')
        sys.exit(2)
    if 'GEDI' in input_file:
        #Check if the water_filtered points file already exists
        Name_parts=re.split('_', input_file.split("/")[-1]) 
        filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_GEDI_WaterSurface_points.shp'
        if os.path.isfile(filtered_file) and overwrite==False:
           print('The water surface points shapefile already exists. Change overwrite option to True')
           return(None)
        #Filter the GEDI file using the defined flags
        gedi_file=gpd.read_file(input_file)
        filtered_data=gedi_file[gedi_file.landsat_WP>80] #Filter by the water persistence column
        filtered_data=filtered_data[filtered_data.noDetctMod==1] #Chose only the points with one detected Mode             
        filtered_data=filtered_data[filtered_data.degradeFlg==0] #Remove points with degrade conditions
        #Remove possible outliers above the defined elevation threshold.
        filtered_data=filtered_data[filtered_data.elevLowMod<elevation_threshold]
        #Select only the relevant columns  
        column_subset=['FILE','BEAM','shot_no','Lat','Lon','time','elevLowMod','DEM','l2b_QF','landsatTC','landsat_WP','urbanProp','surfFl','solarElev','geometry']
        if aditional_columns is not None:
            aditional_columns = [x for x in aditional_columns if x not in column_subset] #Make sure the column is not repeated
            aditional_columns = [x for x in aditional_columns if x in list(filtered_data.columns.values)] #Make sure the column is in the data
            [column_subset.append(y) for y in aditional_columns]
        filtered_data=filtered_data.loc[:,column_subset]    
        filtered_data.to_file(filtered_file)
    if 'GLAS' in input_file:
        #Check if the water_filtered points file already exists
        Name_parts=re.split('_', input_file.split("/")[-1]) 
        filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_GLAS_WaterSurface_points.shp'
        if os.path.isfile(filtered_file) and overwrite==False:
           print('The water surface points shapefile already exists. Change overwrite option to True')
           return(None)
        #Filter the ICESat1 GLAS file using the defined flags
        glas_file=gpd.read_file(input_file)
        filtered_data=glas_file[glas_file.satCorrFlg<3] #Filter acording to the saturation correction flag
        filtered_data=filtered_data[filtered_data.elevUseFlg==0] #Filter acording to the elevation use flag                 
        filtered_data=filtered_data[filtered_data.i_numPk==1]  #Chose only the points with one detected peak
        #Remove no data values from the saturation elevation correction, elevation bias correction, and the delta ellipse values
        filtered_data.satElevCor[:]=filtered_data.satElevCor.replace(1.7976931348623157e+23,0)
        filtered_data.ElevBiasCo[:]=filtered_data.ElevBiasCo.replace(1.7976931348623157e+23,0)
        filtered_data.deltaEllip[:]=filtered_data.deltaEllip.replace(1.7976931348623157e+23,0)
        #calculate the final elevation values
        filtered_data.d_elev[:]=filtered_data.d_elev[:]+filtered_data.satElevCor[:]+filtered_data.ElevBiasCo[:]-filtered_data.deltaEllip[:]
        #Remove possible outliers above the defined elevation threshold.
        filtered_data=filtered_data[filtered_data.d_elev<elevation_threshold]
        #Change the SRMT DEM values from TOPEX/Poseidon ellipsoid to WGS84 ellipsoid.
        filtered_data.d_DEM_elv[:]=filtered_data.d_DEM_elv[:]-filtered_data.deltaEllip[:]
        column_subset=['FILE','id','Lat','Lon','time','d_elev','d_DEM_elv','deltaEllip','geometry']
        if aditional_columns is not None:
            aditional_columns = [x for x in aditional_columns if x not in column_subset] #Make sure the column is not repeated
            aditional_columns = [x for x in aditional_columns if x in list(filtered_data.columns.values)] #Make sure the column is in the data
            [column_subset.append(y) for y in aditional_columns]
        filtered_data=filtered_data.loc[:,column_subset]    
        #Save the filtered data to a temporary file
        temp_file=os.path.dirname(input_file)+'/'+'temporary_file.shp'
        filtered_data.to_file(temp_file)
        #Check if the Global surface water dynamic file exists
        output_directory = os.path.split(os.path.dirname(input_file))[0]+'/Landsat_WP/'
        WP_output_file=output_directory+Name_parts[0]+'_WP.tif'
        if not os.path.isfile(WP_output_file): #Only execute this part if the file doesn't exits 
           #Get the Global surface water dynamic rasters
           locations=gsw.lon_lat_selection(roi_shp)
           for i in range(len(locations)):
               gsw.water_surface_download(locations[i][0], locations[i][1], output_directory)
               gsw.mosaic_construction(output_directory, Name_parts[0]+'_WP.tif')
        filtered_data = gsw.extract_water_persistence(temp_file, WP_output_file)
        filtered_data=filtered_data[filtered_data.landsat_WP>80]
        filtered_data.to_file(filtered_file)
        #remove temporary files
        available_files=[os.path.join(output_directory, file) for file in os.listdir(output_directory)]
        remove=[file for file in available_files if '2018_percent' in file]
        [remove.append(file) for file in available_files if '2019_percent' in file]
        remove.append(temp_file)
        remove.append(re.sub('shp', 'dbf',temp_file))
        remove.append(re.sub('shp', 'cpg',temp_file))
        remove.append(re.sub('shp', 'prj',temp_file))
        remove.append(re.sub('shp', 'shx',temp_file))
        [os.remove(file) for file in remove]
    

