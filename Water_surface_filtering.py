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
import pandas as pd
sys.path.append('/Users/asparra/Documents/Global_deltas/code')
import Global_surface_water_download as gsw


"""
Filtering data points to select only information corresponding to water surface elevation
"""

def water_surface_filtering(roi_shp, input_file, elevation_threshold=None, overwrite=False, aditional_columns=None, extension='csv'):
#Check that the input file and roi shp exists
    if not os.path.isfile(roi_shp):    
        print('error: ROI shpefile provided does not exist or was not found')
        sys.exit(2)
    if not os.path.isfile(input_file):    
        print('error: input file provided does not exist or was not found')
        sys.exit(2)  
    Name_parts=re.split('_', input_file.split("/")[-1]) 
    if 'GEDI' in input_file:
        #Check if the water_filtered points file already exists
        filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_GEDI_WaterSurface_points.'+extension
        if os.path.isfile(filtered_file) and overwrite==False:
           print('The water surface points shapefile already exists. Change overwrite option to True')
           return(None)
        #Filter the GEDI file using the defined flags
        if extension=='shp' or extension=='json':
            gedi_file=gpd.read_file(input_file)
        else: 
            gedi_file=pd.read_csv(input_file)
        filtered_data=gedi_file[gedi_file.landsat_WP>80] #Filter by the water persistence column
        filtered_data=filtered_data[filtered_data.noDetctMod==1] #Chose only the points with one detected Mode             
        filtered_data=filtered_data[filtered_data.degradeFlg==0] #Remove points with degrade conditions
        #Remove possible outliers above the defined elevation threshold.
        if elevation_threshold!= None:
            filtered_data=filtered_data[filtered_data.elevLowMod<elevation_threshold]
        #Select only the relevant columns  
        column_subset=['FILE','BEAM','shot_no','Lat','Lon','time','elevLowMod','DEM','l2b_QF','landsatTC','landsat_WP','urbanProp','surfFl','solarElev','geometry']
        if aditional_columns is not None:
            aditional_columns = [x for x in aditional_columns if x not in column_subset] #Make sure the column is not repeated
            aditional_columns = [x for x in aditional_columns if x in list(filtered_data.columns.values)] #Make sure the column is in the data
            [column_subset.append(y) for y in aditional_columns]
        filtered_data=filtered_data.loc[:,column_subset]    
    else: 
        if 'GLAS' in input_file:
            #Check if the water_filtered points file already exists
            filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_GLAS_WaterSurface_points.'+extension
            if os.path.isfile(filtered_file) and overwrite==False:
                print('The water surface points shapefile already exists. Change overwrite option to True')
                return(None)
            #Filter the ICESat1 GLAS file using the defined flags
            if extension=='shp' or extension=='json':
                glas_file=gpd.read_file(input_file)
            else:
                glas_file=pd.read_csv(input_file)
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
            if elevation_threshold!= None:
                filtered_data=filtered_data[filtered_data.d_elev<elevation_threshold]
            #Change the SRMT DEM values from TOPEX/Poseidon ellipsoid to WGS84 ellipsoid.
            filtered_data.d_DEM_elv[:]=filtered_data.d_DEM_elv[:]-filtered_data.deltaEllip[:]
            column_subset=['FILE','id','Lat','Lon','time','d_elev','d_DEM_elv','deltaEllip','geometry']     
        if 'ATLAS03' in input_file:
            #Check if the water_filtered points file already exists
            filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_ATLAS03_WaterSurface_points.'+extension
            if os.path.isfile(filtered_file) and overwrite==False:
                print('The water surface points shapefile already exists. Change overwrite option to True')
                return(None)
            #Filter the ICESat1 GLAS file using the defined flags
            if extension=='shp' or extension=='json':
                atlas_file=gpd.read_file(input_file)
            else:
                atlas_file=pd.read_csv(input_file)
            filtered_data=atlas_file[atlas_file.podppdFlg==0] #Filter acording to the Composite POD/PPD flag that indicates the quality of input geolocation products for the specific ATL03 segment.
            filtered_data=filtered_data[filtered_data.sc_orient<2]#filter tracks adquired when the spacecraft was in transition. Science quality is potentially degraded while in transition mode.
            #Remove possible outliers above the defined elevation threshold.
            if elevation_threshold!= None:
                filtered_data=filtered_data[filtered_data.h_ph<elevation_threshold]
            column_subset=['FILE','group','chanel', 'count', 'pulse','ph_ndx','Lat','Lon','time','h_ph','dem_h','solarElev','sc_orient','geometry']
        if 'ATLAS08' in input_file:
            #Check if the water_filtered points file already exists
            filtered_file=os.path.dirname(input_file)+'/'+ Name_parts[0] + '_ATLAS08_WaterSurface_points.'+extension
            if os.path.isfile(filtered_file) and overwrite==False:
                print('The water surface points shapefile already exists. Change overwrite option to True')
                return(None)
            #Filter the ICESat1 GLAS file using the defined flags
            if extension=='shp' or extension=='json':
                atlas_file=gpd.read_file(input_file)
            else:
                atlas_file=pd.read_csv(input_file)
            filtered_data=atlas_file[atlas_file.sat_flag==0] #Filter acording to the saturation flag
            filtered_data=filtered_data[filtered_data.snr>1.1]#Filter according to signal to noise ratio
            filtered_data=filtered_data[filtered_data.layer_flag==0]
            filtered_data=filtered_data[filtered_data.msw_flag==0]
            filtered_data=filtered_data[filtered_data.msw_flag==0]
            #Remove possible outliers above the defined elevation threshold.
            if elevation_threshold!= None:
                filtered_data=filtered_data[filtered_data.h_best_fit<elevation_threshold]
            column_subset=['FILE','group','Lat','Lon','time','h_best_fit','dem_h','night_flag','snr','urban_flag','seg_snowcv','geometry']
        if aditional_columns is not None:
            aditional_columns = [x for x in aditional_columns if x not in column_subset] #Make sure the column is not repeated
            aditional_columns = [x for x in aditional_columns if x in list(filtered_data.columns.values)] #Make sure the column is in the data
            [column_subset.append(y) for y in aditional_columns]
        filtered_data=filtered_data.loc[:,column_subset]    
        #Save the filtered data to a temporary file
        temp_file=os.path.dirname(input_file)+'/'+'temporary_file.'+extension
        if extension=='shp' or extension=='json':
            filtered_data.to_file(temp_file)
        else:
            filtered_data.to_csv(temp_file)
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
        #remove temporary files
        wp_files=[os.path.join(output_directory, file) for file in os.listdir(output_directory)]
        remove=[file for file in wp_files if ('2018_percent' in file or '2019_percent' in file)]
        temp_file_names=[temp_file,re.sub(extension, 'dbf',temp_file), re.sub(extension, 'cpg',temp_file),re.sub(extension, 'prj',temp_file),re.sub(extension, 'shx',temp_file)]
        remove.append(temp_file_names)
        remove= [item for sublist in remove for item in sublist]
        [os.remove(file) for file in remove if os.path.isfile(file)]   
    if extension=='shp' or extension=='json':
        filtered_data.to_file(filtered_file)
    else:
        filtered_data.to_csv(filtered_file)

