#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 07:56:36 2021

@author: aparrar1
Code for GEDI, ICESat GLAS and ICESat2 ATLAS data processing, for the Global Delta project
This code defines functions for creating a workflow for querying, downloading, and subsetting GEDI and ICESat data
This code was based on publicly available code and resources.
Original code:
    ICESat query and download: Download script from the NSIDC webpage. Example: https://nsidc.org/data/GLAH06/versions/34 in 'download data' click in 'Download script'.
    GEDI query: GEDI finder web service https://lpdaacsvc.cr.usgs.gov/services/gedifinder
    GEDI subsetting: GEDI_Subsetter.py code by Cole Krehbiel https://git.earthdata.nasa.gov/projects/LPDUR/repos/gedi-subsetter/browse/GEDI_Subsetter.py 
    Earth Data login setup code EarthdataLoginSetup.py created by Cole Krehbiel https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_python/browse/EarthdataLoginSetup.py
#-----------------------------------GEDI L2 data--------------------------------------------------------------------#
GEDI L2 Elevation and Height Metrics Data short name is GEDI02_A and GEDI02_B. there are 2 versions available 001 and 002
The data_group for GEDI are any of the BEAM names ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
GEDI L2A files have 765 science dataset layers per beam. The default layer list only includes 26 layers. 
GEDI L2B files have 251 science dataset layers per beam. The default layer list only includes 38 layers. 
The approximate pixel size for this GEDI data is 25 m. 
#-----------------------------------ICESat GLAS data--------------------------------------------------------------------#
GLAS/ICEsat L2 Global Land Surface Altimetry Data short name is GLAH14, the latest version for this product is version 034
The data_group for GLAS are the HZ groups ['Data_1HZ', 'Data_40HZ']
The approximate pixel size for this ICESat data is 65 m.
#-----------------------------------ICESat2 ATLAS L8 data--------------------------------------------------------------------#
ICESat2 ATLAS L3A Land and Vegetation Height short name is ATL08, the latest version for this product is version 004.
The data_group for ATLAS are any of the ground tracks names ['gt1l', 'gt1r', 'gt2l','gt2r','gt3l','gt3r']
The approximate pixel size for this ICESat2 data is 100 m

"""

import os
import sys
import requests
import itertools
from shapely.geometry import Polygon
import geopandas as gp
from netrc import netrc
import shutil
import h5py
import pandas as pd
import numpy as np
from datetime import timedelta, datetime, timezone
import re
from multiprocessing import Pool, Queue
from io import StringIO
import gc

"""
Functions for file query and data extraction
"""

def ROI_to_BBOX(ROI_file):
#function for getting the bounding box from a shapefile or a list of upper left, lower right coordinates.
#Get the coordinates information from the ROI BBOX
#Make sure the CRS is WGS84, if not, reproject to WGS84
    if "/" in ROI_file:    
        if ROI_file.endswith('json') or ROI_file.endswith('.shp'):
            try:
                ROI = gp.GeoDataFrame.from_file(ROI_file)
                if ROI.crs != 'EPSG:4326':
                    ROI = ROI.to_crs("EPSG:4326")
                if len(ROI) > 1:
                    print('Multi-feature polygon detected. Only the first feature will be used to subset the data.')
                    ROI = ROI.iloc[[0]]
            except:
                print('error: unable to read input geojson file or the file was not found')
                sys.exit(2)
    else:
        if type(ROI_file) == list:
            try:
                ROI = Polygon([(ROI_file[1], ROI_file[0]), (ROI_file[3], ROI_file[0]), (ROI_file[3], ROI_file[2]), (ROI_file[1], ROI_file[2])]) 
                ROI = gp.GeoDataFrame(index=[0], geometry=[ROI]) 
                ROI.crs = 'EPSG:4326'
            except:
                print('error: unable to read input bounding box coordinates, the required format is: ul_lat,ul_lon,lr_lat,lr_lon')
                sys.exit(2)
    return(ROI)

def GEDI_query_wrapper(ROI_file, output_directory, product='GEDI02_A', save_output=True, file_name='GEDI_query'):
#The GEDI query wrapper uses the GEDI finder web service and only works with version 001 products. For querying GEDI v002 products use the CMR_query_wrapper
#Get the coordenates of a ROI, and construct a GEDI Finder web service URL 
#Get the names of the files that cover the ROI

#Get the coordinates information from the ROI BBOX
#Make sure the CRS is WGS84, if not, reproject to WGS84
    BBOX=ROI_to_BBOX(ROI_file)
    lon=[BBOX[0],BBOX[2]] 
    lat=[BBOX[1],BBOX[3]]       
    BBOX=str([max(lat), min(lon), min(lat), max(lon)])  

#Construct the URL link
    url_base='https://lpdaacsvc.cr.usgs.gov/services/gedifinder?product='
    url_name=url_base+product+'&version=001&bbox='+BBOX+'&output=json'

#Get the available files that overlap the ROI and save or return the list
    DAAC_response= requests.get(url_name).text
    response_dict= json.loads(DAAC_response)
    available_files=list(response_dict["data"])
    if save_output==True:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        save_file=output_directory+'/'+file_name
        GEDI_files=open(save_file, "w") 
        GEDI_files.write(str(available_files))
        GEDI_files.close()
    print("A total of "+str(len(available_files))+' GEDI files overlap the ROI')
    return(available_files)

def cmr_query_wrapper(ROI_file, output_directory, short_name='', version='', time_start='', time_end='',save_output=True, file_name='file_query'):
#Code for querying the Common Metadata Repository (CMR). This metadata system catalogs the Earth Science data and associated service metadata records.
#-----------------------------Step 1------------------------------
#Make sure the parameters for short_name and version are included
    if short_name=='' or version=='':
        print('error: short_name and/or version not defined')
        sys.exit(2)
#-----------------------------Step 2------------------------------
#Get the coordinates information from the ROI BBOX
#Make sure the CRS is WGS84, if not, reproject to WGS84
    BBOX=ROI_to_BBOX(ROI_file)
    BBOX=BBOX.total_bounds
    lon=[BBOX[0],BBOX[2]] 
    lat=[BBOX[1],BBOX[3]]
#The bounding box should be formated west south east north       
    BBOX=[min(lon), min(lat), max(lon), max(lat)]
    final_bbox=str(BBOX[0])+","+str(BBOX[1])+","+str(BBOX[2])+","+str(BBOX[3])
#---------------------------------Step 3---------------------------------
#Construct the URL link
    CMR_URL = 'https://cmr.earthdata.nasa.gov'
    CMR_PAGE_SIZE = 2000
    if 'GEDI' in short_name:
        provider='LPDAAC_ECS'
    else:
        provider='NSIDC_ECS'
    CMR_FILE_URL = ('{0}/search/granules.json?provider={1}'
                '&sort_key[]=start_date&sort_key[]=producer_granule_id'
                '&scroll=true&page_size={2}'.format(CMR_URL,provider, CMR_PAGE_SIZE))
    params = '&short_name={0}'.format(short_name)
    params += '&version={0}'.format(version)
    params += '&temporal[]={0},{1}'.format(time_start, time_end)
    params += '&bounding_box={0}'.format(final_bbox)

    cmr_query_url =CMR_FILE_URL + params
    response= requests.get(cmr_query_url)
    cmr_info=response.json()
    cmr_info=cmr_info['feed']
    results_list=cmr_info['entry']
#---------------------------------Step 4---------------------------------
#Get the available files that overlap the ROI and save or return the list
    entries = [e['links']
               for e in results_list
               if 'links' in e]
    # Flatten "entries" to a simple list of links
    links = list(itertools.chain(*entries))

    urls = []
    unique_filenames = set()
   
    for link in links:
        if 'href' not in link:
            # Exclude links with nothing to download
            continue
        if 'inherited' in link and link['inherited'] is True:
            # This refers to links that indicate higher level pages, not the links to the files
            continue
        if 'rel' in link and 'data#' not in link['rel']:
            # Exclude links which are not classified by CMR as "data" or "metadata"
            continue
        if 'title' in link and 'opendap' in link['title'].lower():
            # Exclude OPeNDAP links--they are responsible for many duplicates
            # This is a hack; when the metadata is updated to properly identify
            # non-datapool links, we should be able to do this in a non-hack way
            continue
        filename = link['href'].split('/')[-1]
        if filename in unique_filenames:
            # Exclude links with duplicate filenames (they would overwrite)
            continue
        unique_filenames.add(filename)
        urls.append(link['href'])

    if save_output==True:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        save_file=output_directory+'/'+file_name
        avail_files=open(save_file, "w") 
        avail_files.write(str(urls))
        avail_files.close()
    url_h5=[e for e in urls if 'h5' in e.lower() and not "MET" in e]
    print("A total of "+str(len(url_h5))+' files overlap the ROI')
    return(urls)

def earthdata_download(f, temporary_directory):
#Code for downloading the available files for the ROI
#A netrc_file needs to be created to be able to download the files. 
#Follow the instructions on the Earth Data login setup code EarthdataLoginSetup.py created by Cole Krehbiel to construct the netrc_file:
#https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_python/browse/EarthdataLoginSetup.py
    urs = 'urs.earthdata.nasa.gov' 
    netrc_file=os.path.expanduser("~/.netrc")
    if not os.path.exists(temporary_directory):
        os.makedirs(temporary_directory)
    saveName = os.path.join(temporary_directory, f.split('/')[-1].strip())   
    # Create and submit request and download file
    with requests.get(f.strip(), verify=False, stream=True, auth=(netrc(netrc_file).authenticators(urs)[0], netrc(netrc_file).authenticators(urs)[2])) as response:
        if response.status_code != 200:
            print("{} not downloaded. Verify that your username and password are correct in {}".format(f.split('/')[-1].strip(), netrc_file))
        else:
            response.raw.decode_content = True
            content = response.raw
            with open(saveName, 'wb') as d:
                shutil.copyfileobj(content, d)
            print('Downloaded file: {}'.format(saveName))

def GEDI_group_process(file1Name,b,file1, file1SDS, sdsSubset,ROI):
#Function for extracting the data layers from GEDI HDF5 files
    groupSDS = [s for s in file1SDS if b in s]   
    # Search for delta time, latitude, longitude, and shot number SDS
    d_time=[l for l in groupSDS if 'delta_time' in l][0]
    lat=[l for l in groupSDS if 'lat_lowestmode' in l][0] 
    lon = [l for l in groupSDS if 'lon_lowestmode' in l][0]
    shot = [l for l in groupSDS if 'shot_number' in l][0]        
        
    # Open latitude, longitude, and shot number SDS
    times = file1[d_time][()]
    shots = file1[shot][()]
    lats = file1[lat][()]
    lons = file1[lon][()]
    #the delta time is defined as seconds since 2018-01- 01. transform to complete date and time
    initial_date=datetime(2018, 1, 1,12,0,0, tzinfo=timezone.utc)
    delta = [str(initial_date+timedelta(seconds=l)) for l in times]
    delta=np.array(delta)
       
    # Append BEAM, shot number, latitude, longitude and an index to the GEDI dataframe
    geoDF = pd.DataFrame({'FILE': len(shots) * [file1Name], 'BEAM': len(shots) * [b], 'time':delta, 'shot_no': shots, 
                              'Lat':lats, 'Lon':lons })
    del times, shots, lats, lons, delta
    #OPEN SDS AND APPEND TO GEODATAFRAME
    groupSDS = [s for s in groupSDS if b in s and not any(s.endswith(d) for d in ['delta_time','lat_lowestmode','lon_lowestmode','shot_number'])]
    for s in groupSDS:
        sName = s.split('/')[-1].replace('/', '_')
        #Change the science dataset layer's names to shorter versions
        sName = re.sub('digital_elevation_model', 'DEM', sName)
        sName = re.sub('elev_highestreturn', 'elevHighRe', sName)
        sName = re.sub('elev_lowestmode', 'elevLowMod', sName)        
        sName = re.sub('elevation_bin0', 'elevBin0', sName)
        sName = re.sub('elevation_bin0_error', 'elevBn0Er', sName)
        sName = re.sub('elevation_lastbin', 'elevLstBin', sName)
        sName = re.sub('elevation_lastbin_error', 'elevLBnEr', sName)
        sName = re.sub('height_bin0', 'heightBin0', sName)
        sName = re.sub('height_lastbin', 'heightLBin', sName)
        sName = re.sub('leaf_off_flag', 'leafOffFlg', sName)
        sName = re.sub('modis_nonvegetated', 'MODnonveg', sName)
        sName = re.sub('urban_proportion', 'urbanProp', sName)
        sName = re.sub('num_detectedmodes', 'noDetctMod', sName)
        sName = re.sub('_treecover', 'TC', sName)
        sName = re.sub('water_persistence', 'WP', sName)
        sName = re.sub('algorithmrun_flag', 'algoRunFlg', sName)
        sName = re.sub('selected_rg_algorithm', 'selecRGalg', sName)
        sName = re.sub('_quality_flag', '_QF', sName)
        sName = re.sub('sensitivity', 'sensitivit', sName)
        sName = re.sub('stale_return_flag', 'staleRetFl', sName)
        sName = re.sub('surface_flag', 'surfFl', sName)
        sName = re.sub('degrade_flag', 'degradeFlg', sName)
        sName = re.sub('solar_elevation', 'solarElev', sName)
        sName = re.sub('mean_sea_surface', 'meanSeaSur', sName)
        sName = re.sub('selected_mode', 'selectMode', sName)
        sName = re.sub('quality_flag', 'QF', sName)
        sName = re.sub('elevation_bias_flag', 'elevBiasFl', sName)
        sName = re.sub('selected_algorithm', 'selectAlgo', sName)
        
        # Datasets with consistent structure as shots
        if file1[s].shape == file1[shot].shape:
            temp = pd.DataFrame()
            temp[sName] = file1[s][()] 
            geoDF=pd.concat((geoDF, temp), axis=1)
                        
        # Datasets with a length of one 
        elif len(file1[s][()]) == 1:
            temp = pd.DataFrame()
            temp[sName] = [file1[s][()][0]] * len(geoDF) # create array of same single value
            geoDF=pd.concat((geoDF, temp), axis=1)
                
        # Multidimensional datasets
        elif len(file1[s].shape) == 2 and 'surface_type' not in s: 
            allData = file1[s][()][()]
                    
        # For each additional dimension, create a new output column to store those data
            for i in range(file1[s].shape[1]):
                temp = pd.DataFrame()
                temp[f"{sName}_{i}"] = allData[:,i]
                geoDF=pd.concat((geoDF, temp), axis=1)
                
        # Waveforms
        elif s.endswith('waveform') or s.endswith('pgap_theta_z'):
            waveform = []
                    
            if s.endswith('waveform'):
                # Use sample_count and sample_start_index to identify the location of each waveform
                start = file1[f'{b}/{s.split("/")[-1][:2]}_sample_start_index'][()]
                count = file1[f'{b}/{s.split("/")[-1][:2]}_sample_count'][()]
                    
            # for pgap_theta_z, use rx sample start index and count to subset
            else:
                # Use sample_count and sample_start_index to identify the location of each waveform
                start = file1[f'{b}/rx_sample_start_index'][()]
                count = file1[f'{b}/rx_sample_count'][()]
            wave = file1[s][()]
                    
            # in the dataframe, each waveform will be stored as a list of values
            for k in range(len(start)):
                singleWF = wave[int(start[k] - 1): int(start[k] - 1 + count[k])]
                waveform.append(','.join([str(q) for q in singleWF]))
            geoDF[sName] = waveform
                
        # Surface type 
        elif s.endswith('surface_type'):
            surfaces = ['land', 'ocean', 'sea_ice', 'land_ice', 'inland_water']
            allData = file1[s][()]
            for i in range(file1[s].shape[0]):
                geoDF[f'{surfaces[i]}'] = allData[i][()]
        else:
            print(f"SDS: {s} not found")
        
    # Convert lat/lon coordinates to shapely points and append to geodataframe
    geoDF2 = geoDF.copy()
    geoDF = gp.GeoDataFrame(geoDF2, geometry=gp.points_from_xy(geoDF2.Lon, geoDF2.Lat))
    geoDF.crs='EPSG:4326'
    # Clip to only include points within the user-defined bounding box
    geoDF = geoDF[geoDF['geometry'].within(ROI.geometry[0].envelope)]
    return(geoDF)

def GLAS_group_process(file1Name, b,file1, file1SDS, sdsSubset,ROI):
#Function for extracting the data layers from ICESat GLAS HDF5 files
    groupSDS = [s for s in file1SDS if b in s]
    if b=='Data_1HZ':
        date = '/Data_1HZ/DS_UTCTime_1'
    else:
        date = '/Data_40HZ/DS_UTCTime_40' 
   
    # Search for latitude, longitude, and shot number SDS
    lat = [l for l in groupSDS if 'lat' in l][0]  
    lon = [l for l in groupSDS if 'lon' in l][0]
    index = [l for l in groupSDS if 'i_rec_ndx' in l][0]
    shot = [l for l in groupSDS if 'i_shot_count' in l][0]         
    
    # Open the date, index number, shot number, latitude, and longitude
    times = file1[date][()]
    indices = file1[index][()]
    shots = file1[shot][()]
    lats = file1[lat][()]
    lons = file1[lon][()]  
    initial_date=datetime(2000, 1, 1,12,0,0, tzinfo=timezone.utc)
    delta = [str(initial_date+timedelta(seconds=l)) for l in times]
    delta=np.array(delta)  
    
    #Check if the longitud data is in the right coordinate format
    #The GLAS datasets sometimes have the longitud values from 0 to 360, instead of using negative values
    if max(lons[lons!=max(lons)])>180:
        lons[lons!=1.7976931348623157e+308]=((lons[lons!=1.7976931348623157e+308]+180)%360)-180 #1.7976931348623157e+308 is the no data value
    
    #Append data group, date, index number, shot number, latitude, and longitude to the ICESat dataframe
    geoDF = pd.DataFrame({'FILE': len(shots) * [file1Name],'group': len(shots) * [b],'time':delta, 'id': indices+shots, 
                              'Lat':lats, 'Lon':lons})
    del times, indices, shots, lats, lons, delta
    groupSDS=[s for s in groupSDS if not any(s.endswith(d) for d in ['lat','lon','i_rec_ndx','i_shot_count'])]
         
    for s in groupSDS:
        sName = s.split('/')[-1].replace('/', '_')
        
        #Change the science dataset layer's names to shorter versions
        sName = re.sub('d_SigBegOff', 'SigBegOff', sName) 
        sName = re.sub('d_SigEndOff', 'SigEndOff', sName) 
        sName = re.sub('d_gpCntRngOff', 'gpCntRng', sName) 
        sName = re.sub('d_deltaEllip', 'deltaEllip', sName) 
        sName = re.sub('d_satElevCorr', 'satElevCor', sName)
        sName = re.sub('d_ElevBiasCorr', 'ElevBiasCo', sName)
        sName = re.sub('sat_corr_flg', 'satCorrFlg', sName)
        sName = re.sub('elev_use_flg', 'elevUseFlg', sName)   
        sName = re.sub('sigma_att_flg', 'sigmaAttFl', sName)   
         
        # Datasets with consistent structure as shots
        if file1[s].shape == file1[shot].shape:
            temp = pd.DataFrame()
            temp[sName] = file1[s][()] 
            geoDF=pd.concat((geoDF, temp), axis=1)
            
        # Datasets with a length of one 
        elif len(file1[s][()]) == 1:
            temp = pd.DataFrame()
            temp[sName] = [file1[s][()][0]] * len(geoDF) # create array of same single value
            geoDF=pd.concat((geoDF, temp), axis=1)
           
        # Multidimensional datasets
        elif len(file1[s].shape) == 2: 
            allData = file1[s][()][()]
                
            # For each additional dimension, create a new output column to store those data
            for i in range(file1[s].shape[1]):
                temp = pd.DataFrame()
                temp[f"{sName}_{i}"] = allData[:,i]
                geoDF=pd.concat((geoDF, temp), axis=1)
            
        else:
            print(f"SDS: {s} not found")
    # Convert lat/lon coordinates to shapely points and append to geodataframe
    geoDF2 = geoDF.copy()
    geoDF = gp.GeoDataFrame(geoDF2, geometry=gp.points_from_xy(geoDF2.Lon, geoDF2.Lat))
    geoDF.crs='EPSG:4326'
    #remove points with no data values for Lat and Lon
    geoDF=geoDF[geoDF['Lat']!=1.7976931348623157e+308]
    # Clip to only include points within the user-defined bounding box
    geoDF = geoDF[geoDF['geometry'].within(ROI.geometry[0].envelope)]
    return(geoDF)

def ATLAS08_group_process(file1Name, b,file1, file1SDS, sdsSubset,ROI):
#Function for extracting the data layers from ICESat2 ATLAS HDF5 files
    groupSDS = [s for s in file1SDS if b in s]
    # Search for latitude, longitude, and shot number SDS
    lat = [l for l in groupSDS if 'latitude' in l][0]  
    lon = [l for l in groupSDS if 'longitude' in l][0]
    date= [l for l in groupSDS if 'delta_time' in l][0]   
    
    # Open the date, index number, shot number, latitude, and longitude
    times = file1[date][()]
    lats = file1[lat][()]
    lons = file1[lon][()]  
    initial_date=datetime(2018, 1, 1,0,0,0, tzinfo=timezone.utc)
    delta = [str(initial_date+timedelta(seconds=l)) for l in times]
    delta=np.array(delta)  
     
    #Append data group, date, index number, shot number, latitude, and longitude to the ICESat dataframe
    geoDF = pd.DataFrame({'FILE': len(times) * [file1Name],'group': len(times) * [b],'time':delta, 
                              'Lat':lats, 'Lon':lons})
    del times, lats, lons, delta
    groupSDS=[s for s in groupSDS if not any(s.endswith(d) for d in ['latitude','longitude','delta_time'])] 
         
    for s in groupSDS:
        sName = s.split('/')[-1].replace('/', '_')
        
        #Change the science dataset layer's names to shorter versions
        sName = re.sub('h_te_best_fit', 'h_best_fit', sName) 
        sName = re.sub('h_te_interp', 'h_interp', sName) 
        sName = re.sub('h_te_median', 'h_median', sName)     
        sName = re.sub('h_te_uncertainty', 'h_uncerta', sName)      
        sName = re.sub('subset_te_flag', 'subset_flg', sName)     
        sName = re.sub('terrain_slope', 'slope', sName)
        sName = re.sub('dem_removal_flag', 'dem_remv_f', sName)
        sName = re.sub('ph_removal_flag', 'ph_remv_fl', sName)      
        sName = re.sub('segment_landcover', 'seg_landcv', sName)   
        sName = re.sub('segment_snowcover', 'seg_snowcv', sName)   
        sName = re.sub('segment_watermask', 'watermask', sName)   
        sName = re.sub('canopy_flag', 'canopy_flg', sName)   
        sName = re.sub('canopy_openness', 'canopyOpen', sName)   
        sName = re.sub('h_canopy_abs', 'h_canopy_a', sName)  
        sName = re.sub('canopy_openness', 'canopyOpen', sName)  
        sName = re.sub('h_canopy_uncertainty', 'hCanopyUnc', sName)  
        sName = re.sub('h_max_canopy', 'hMaxCanopy', sName) 
        sName = re.sub('h_canopy_uncertainty', 'hCanopyUnc', sName)     
        sName = re.sub('h_mean_canopy', 'hMeanCanop', sName) 
        sName = re.sub('h_min_canopy', 'hMinCanopy', sName) 
        sName = re.sub('landsat_flag', 'landsatFlg', sName) 
        sName = re.sub('landsat_perc', 'landsatPer', sName)
        sName = re.sub('subset_can_flag', 'subsCanFlg', sName) 
        sName = re.sub('toc_roughness', 'toc_roughn', sName)          
         
        # Datasets with consistent structure as shots
        if file1[s].shape == file1[lat].shape:
            temp = pd.DataFrame()
            temp[sName] = file1[s][()] 
            geoDF=pd.concat((geoDF, temp), axis=1)
            
        # Datasets with a length of one 
        elif len(file1[s][()]) == 1:
            temp = pd.DataFrame()
            temp[sName] = [file1[s][()][0]] * len(geoDF) # create array of same single value
            geoDF=pd.concat((geoDF, temp), axis=1)
           
        # Multidimensional datasets
        elif len(file1[s].shape) == 2: 
            allData = file1[s][()][()]
                
            # For each additional dimension, create a new output column to store those data
            for i in range(file1[s].shape[1]):
                temp = pd.DataFrame()
                temp[f"{sName}_{i}"] = allData[:,i]
                geoDF=pd.concat((geoDF, temp), axis=1)
            
        else:
            print(f"SDS: {s} not found")
      
    # Convert lat/lon coordinates to shapely points and append to geodataframe
    geoDF2 = geoDF.copy()
    geoDF = gp.GeoDataFrame(geoDF2, geometry=gp.points_from_xy(geoDF2.Lon, geoDF2.Lat))
    geoDF.crs='EPSG:4326'
    # Clip to only include points within the user-defined bounding box
    geoDF = geoDF[geoDF['geometry'].within(ROI.geometry[0].envelope)]
    return(geoDF)

def ATLAS03_group_process(file1Name, b,file1, file1SDS, sdsSubset,ROI):
#Function for extracting the data layers from ICESat2 ATLAS HDF5 files
    groupSDS = [s for s in file1SDS if b in s]
    # Search for latitude, longitude, and photon id number SDS
    #Get the information at the photon level  
    lat = [l for l in groupSDS if 'lat_ph' in l][0]  
    lon = [l for l in groupSDS if 'lon_ph' in l][0]
    date= [l for l in groupSDS if 'delta_time' in l][0]   
    channel = [l for l in groupSDS if 'ph_id_channel' in l][0]  
    count = [l for l in groupSDS if 'ph_id_count' in l][0]  
    pulse = [l for l in groupSDS if 'ph_id_pulse' in l][0] 
    
    # Open the date, index number, id number, latitude, and longitude
    times = file1[date][()]
    lats = file1[lat][()]
    lons = file1[lon][()]  
    initial_date=datetime(2018, 1, 1,0,0,0, tzinfo=timezone.utc)
    delta = [str(initial_date+timedelta(seconds=l)) for l in times]
    delta=np.array(delta)
    channels = file1[channel][()]
    counts = file1[count][()]
    pulses = file1[pulse][()]
    #ph_id=channels.astype(str)+counts.astype(str)+pulses.astype(str)
    #Append data group, date, index number, shot number, latitude, and longitude to the ICESat dataframe
    photonDF = pd.DataFrame({'FILE': len(times) * [file1Name],'group': len(times) * [b],'time':delta,'Lat':lats, 'Lon':lons, 'chanel': channels, 
                             'count': counts,'pulse': pulses,'ph_ndx': np.arange(1, len(lats)+1)})
    del times, lats, lons, delta, channels, counts, pulses
    groupSDS=[s for s in groupSDS if not any(s.endswith(d) for d in ['lat_ph','lon_ph','delta_time','ph_id_channel','ph_id_count','ph_id_pulse'])]   
    #Separate the remaining parameters in 20m segment paramters and individual photon parameters 
    part1 =[s for s in groupSDS if 'heights' in s]
    part2 =[s for s in groupSDS if not 'heights' in s] 
    #Get the information at the photon level  
    for s in part1:
        sName = s.split('/')[-1].replace('/', '_')     
        #Change the science dataset layer's names to shorter versions
        sName = re.sub('signal_conf_ph', 'sgnlConf', sName)       
        # Datasets with consistent structure as shots
        if file1[s].shape == file1[lat].shape:
            temp = pd.DataFrame()
            temp[sName] = file1[s][()] 
            photonDF=pd.concat((photonDF, temp), axis=1)
            del temp
        # Datasets with a length of one 
        elif len(file1[s][()]) == 1:
            temp = pd.DataFrame()
            temp[sName] = [file1[s][()][0]] * len(photonDF) # create array of same single value
            photonDF=pd.concat((photonDF, temp), axis=1)
            del temp
        # Multidimensional datasets
        elif len(file1[s].shape) == 2: 
            allData = file1[s][()][()]    
            # For each additional dimension, create a new output column to store those data
            for i in range(file1[s].shape[1]):
                temp = pd.DataFrame()
                temp[f"{sName}_{i}"] = allData[:,i]
                photonDF=pd.concat((photonDF, temp), axis=1)   
                del temp
        else:
            print(f"SDS: {s} not found")
    
    #Get the information for the 20m segments
    # Search for the photon index and counts per segment
    ph_ndx = [l for l in part2 if 'ph_index_beg' in l][0]  
    ph_count = [l for l in part2 if 'segment_ph_cnt' in l][0]
    rgt = 'orbit_info/rgt'
    seg_id= [l for l in part2 if 'segment_id' in l][0]
    # Open the photon index begining and the photon count per segment
    ph_index = file1[ph_ndx][()]
    counts = file1[ph_count][()]
    rgts= file1[rgt][()]
    seg_ids = file1[seg_id][()]    

    #create a separate data frame with the reference ground track, the segment id, the photon index beginning and end
    geoDF = pd.DataFrame({'rgt':[rgts[0]]*len(ph_index),'segm_id':seg_ids,'index_beg':ph_index, 'index_end':ph_index+counts-1, 'ph_count':counts})
    del ph_index, counts, rgts, seg_ids
    part2=[s for s in part2 if not any(s.endswith(d) for d in ['ph_index_beg','segment_ph_cnt','segment_id'])]
     
    for s in part2:
        sName = s.split('/')[-1].replace('/', '_')     
        sName = re.sub('full_sat_fract', 'fullSatFrc', sName) 
        sName = re.sub('near_sat_fract', 'nearSatFrc', sName) 
        sName = re.sub('podppd_flag', 'podppdFlg', sName) 
        sName = re.sub('reference_photon_index', 'refPhIndx', sName) 
        sName = re.sub('solar_elevation', 'solarElev', sName) 
        sName = re.sub('geoid_free2mean', 'geoFr2mea', sName) 
        sName = re.sub('tide_earth_free2mean', 'tideEFr2me', sName) 
        sName = re.sub('tide_equilibrium', 'tideEquili', sName) 
        sName = re.sub('tide_oc_pole', 'tideOCpole', sName) 
        sName = re.sub('surf_type', 'surfType', sName) 
        
        # Datasets with consistent structure as shots
        if file1[s].shape == file1[ph_ndx].shape:
            temp = pd.DataFrame()
            temp[sName] = file1[s][()] 
            geoDF =pd.concat((geoDF, temp), axis=1)
            del temp
        # Datasets with a length of one 
        elif len(file1[s][()]) == 1:
            temp = pd.DataFrame()
            temp[sName] = [file1[s][()][0]] * len(geoDF) # create array of same single value
            geoDF =pd.concat((geoDF, temp), axis=1)
            del temp
        # Multidimensional datasets
        elif len(file1[s].shape) == 2: 
            allData = file1[s][()][()]   
            # For each additional dimension, create a new output column to store those data
            for i in range(file1[s].shape[1]):
                temp = pd.DataFrame()
                temp[f"{sName}_{i}"] = allData[:,i]
                geoDF =pd.concat((geoDF, temp), axis=1)
                del temp
        else:
            print(f"SDS: {s} not found")
    
    #Extend the geoDF dataframe so that the number of rows matches the number of rows of the photonDF
    #Extend the geoDF using the ph_count, and add a new photon index column
    geoDF2=geoDF.loc[geoDF.index.repeat(geoDF.ph_count)]
    geoDF2['ph_ndx']=np.arange(1, len(geoDF2)+1)   
    #Merge the 2 dataframes using the 'ph_ndx' column
    df_list = [photonDF, geoDF2]
    for df in df_list:
        df.set_index(['ph_ndx'], inplace=True)
    finalDF = pd.concat(df_list, axis=1) 
    finalDF.reset_index(inplace=True)
    #Remove the points from the finalDF that do not overlap the ROI BBOX
    BBOX=ROI.total_bounds
    lon=[BBOX[0],BBOX[2]] 
    lat=[BBOX[1],BBOX[3]]
    #Remove the points that fall outside the ROI bounding box
    finalDF=finalDF[finalDF.Lat>=min(lat)]
    finalDF=finalDF[finalDF.Lat<=max(lat)]
    finalDF=finalDF[finalDF.Lon>=min(lon)]
    finalDF=finalDF[finalDF.Lon<=max(lon)]
    # Convert lat/lon coordinates to shapely points
    finalDF = gp.GeoDataFrame(finalDF, geometry=gp.points_from_xy(finalDF.Lon, finalDF.Lat))
    finalDF.crs='EPSG:4326'
    #Further filter the data according to quality flags quality_ph and the the signal_conf_ph
    #Only use values with quality_ph==0 (nominal) and Events associated with the land surface type high confidence
    finalDF=finalDF[finalDF.quality_ph==0]
    finalDF=finalDF[finalDF.surfType_0==1]
    finalDF=finalDF[finalDF.sgnlConf_0==4]
    del geoDF2
    del photonDF
    return(finalDF)
     
def HDF5_data_subset(ROI_file,input_file,output_directory,data_group ='',layerSubset=None, file_ext='shp'):
#Wrapper function for data extraction for GEDI, ICESat, and ICESat2 HDF5 files    
#-----------------------------Step 1------------------------------
#Make sure the parameter for data_group is included
    if data_group=='':
        print('error: the data_group parameter is not defined')
        sys.exit(2)
#---------------------------------Step 2-----------------------------------------------------------------------#
#Get the extent information from the ROI shapefile or coordinates input
    ROI=ROI_to_BBOX(ROI_file) 
# Keep the exact input geometry for the final clip to ROI
    finalClip = gp.GeoDataFrame([1], geometry=[ROI.geometry[0]], crs='EPSG:4326')    
#--------------------------------Step 3----------------------------------------------------------------------------#
#Make sure that the input file and the output directory exist
    if not os.path.isfile(input_file) and not input_file.lower().endswith('.h5') and not any('GLA' or 'ATL' or 'GEDI' in input_file):    
        print('error: input file provided does not exist or was not found, or is not a GEDI or ICESat GLAS/ATLAS file')
        sys.exit(2)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
# --------------------------------Step 4---------------------------------------------------------------------------#
#Define the default Science Dataset layers to be subset and exported, see the different product data dictionaries for information on available file layers

#GEDI L2A product
    l2aSubset = ['/lat_lowestmode', '/lon_lowestmode', '/channel', '/shot_number', '/delta_time', 
                 '/digital_elevation_model','/digital_elevation_model_srtm', '/elev_lowestmode', 'elev_highestreturn','mean_sea_surface','rh100'
                 '/sensitivity','selected_mode', '/quality_flag', '/degrade_flag', '/elevation_bias_flag', '/surface_flag',  
                  '/num_detectedmodes',  '/selected_algorithm',  '/solar_elevation','/land_cover_data/landsat_treecover',
                  '/land_cover_data/landsat_water_persistence', '/land_cover_data/leaf_off_flag','/land_cover_data/modis_nonvegetated',
                  '/land_cover_data/modis_treecover','/land_cover_data/pft_class', '/land_cover_data/urban_proportion']   
#GEDI L2B product    
    l2bSubset = ['/geolocation/lat_lowestmode', '/geolocation/lon_lowestmode', '/channel', '/geolocation/shot_number','/geolocation/delta_time',
                 '/geolocation/digital_elevation_model', '/geolocation/elev_highestreturn','/geolocation/elev_lowestmode','/geolocation/elevation_bin0',
                 '/geolocation/elevation_bin0_error','/geolocation/elevation_lastbin', '/geolocation/elevation_lastbin_error','/geolocation/height_bin0','/geolocation/height_lastbin', 
                 '/cover', '/fhd_normal','/omega', '/pai', '/rhov',  '/rhog', '/rh100',
                 '/land_cover_data/landsat_treecover', '/land_cover_data/landsat_water_persistence', '/land_cover_data/leaf_off_flag', '/land_cover_data/modis_nonvegetated', 
                 '/land_cover_data/modis_treecover', '/land_cover_data/pft_class', '/land_cover_data/urban_proportion',
                 '/algorithmrun_flag','/num_detectedmodes','/selected_rg_algorithm','/l2a_quality_flag', '/l2b_quality_flag',  '/sensitivity',  
                 '/stale_return_flag', '/surface_flag', '/geolocation/degrade_flag',  '/geolocation/solar_elevation']
#ICESat GLAS products
    GLASSubset = ['/Geolocation/d_lat', '/Geolocation/d_lon', '/Time/i_rec_ndx','/Time/i_shot_count',
                  '/Elevation_Surfaces/d_elev','/Elevation_Surfaces/d_refRng',
                  '/Elevation_Offsets/d_ldRngOff','/Elevation_Offsets/d_SigBegOff','/Elevation_Offsets/d_SigEndOff','/Elevation_Offsets/d_gpCntRngOff',
                  '/Geophysical/d_DEM_elv','/Geophysical/d_deltaEllip','/Geophysical/d_erElv','/Geophysical/d_gdHt','/Geophysical/d_ldElv','/Geophysical/d_ocElv','/Geophysical/d_poTide',    
                  '/Elevation_Corrections/d_satElevCorr','/Elevation_Corrections/d_ElevBiasCorr','/Elevation_Corrections/d_GmC',
                  '/Quality/sat_corr_flg', '/Quality/elev_use_flg', '/Waveform/i_numPk','/Quality/sigma_att_flg','/Reflectivity/d_reflctUC']    
#ICESat2 ATLAS08 product
    ATLAS8Subset = ['/land_segments/latitude', '/land_segments/longitude', '/land_segments/delta_time',
                   '/land_segments/terrain/h_te_best_fit','/land_segments/terrain/h_te_interp','/land_segments/terrain/h_te_max','/land_segments/terrain/h_te_mean',
                   '/land_segments/terrain/h_te_median','/land_segments/terrain/h_te_min','/land_segments/terrain/h_te_mode','/land_segments/terrain/h_te_uncertainty',
                   '/land_segments/terrain/subset_te_flag','/land_segments/terrain/terrain_slope',
                   '/land_segments/dem_flag','/land_segments/dem_h','/land_segments/dem_removal_flag','/land_segments/h_dif_ref','/land_segments/layer_flag',
                   '/land_segments/msw_flag','/land_segments/n_seg_ph','/land_segments/night_flag','/land_segments/ph_removal_flag',
                   '/land_segments/sat_flag','/land_segments/segment_landcover','/land_segments/segment_snowcover','/land_segments/segment_watermask',
                   '/land_segments/snr','/land_segments/surf_type','/land_segments/terrain_flg','/land_segments/urban_flag',
                    '/land_segments/canopy/canopy_flag','/land_segments/canopy/canopy_openness','/land_segments/canopy/h_canopy',
                    '/land_segments/canopy/h_canopy_abs','/land_segments/canopy/h_canopy_uncertainty','/land_segments/canopy/h_max_canopy',
                    '/land_segments/canopy/h_mean_canopy','/land_segments/canopy/h_min_canopy','/land_segments/canopy/landsat_flag',
                    '/land_segments/canopy/landsat_perc','/land_segments/canopy/subset_can_flag','/land_segments/canopy/toc_roughness']
#ICESat2 ATLAS03 product
    ATLAS3Subset = ['/heights/lat_ph', '/heights/lon_ph', '/heights/delta_time',
                   '/heights/h_ph','/heights/ph_id_channel','/heights/ph_id_count','/heights/ph_id_pulse',
                   '/heights/quality_ph','/heights/signal_conf_ph',
                   '/geolocation/full_sat_fract','/geolocation/near_sat_fract','/geolocation/ph_index_beg','/geolocation/podppd_flag',
                   '/geolocation/reference_photon_index','/geolocation/segment_id','/geolocation/segment_ph_cnt',
                   '/geolocation/solar_elevation','/geolocation/surf_type',
                   '/geophys_corr/dac', '/geophys_corr/dem_flag','/geophys_corr/dem_h','/geophys_corr/geoid',
                   '/geophys_corr/geoid_free2mean','/geophys_corr/tide_earth','/geophys_corr/tide_earth_free2mean','/geophys_corr/tide_equilibrium',
                   '/geophys_corr/tide_load','/geophys_corr/tide_oc_pole','/geophys_corr/tide_ocean','/geophys_corr/tide_pole']
    
# ----------------------------------Step 5------------------------------------------------------------------------- #   
#Get the science datasets from each file and subset the available points to the ROI
    print(f"Processing file: {input_file}")
    file1 = h5py.File(input_file, 'r')      # Open file
    file1Name = re.split(".h5", input_file.split("/")[-1], flags=re.IGNORECASE)[0] # Keep original filename 
    file1_objs = []            
    file1.visit(file1_objs.append)  # Retrieve list of datasets  

    # Search for relevant Science Datasets (SDS) inside data file
    file1SDS = [str(o) for o in file1_objs if isinstance(file1[o], h5py.Dataset)] #Make sure the available names are layers and not group headers
    
    # Define subset of layers based on product
    if 'GEDI02_A' in input_file:
        sdsSubset = l2aSubset 
    elif 'GEDI02_B' in input_file:
        sdsSubset = l2bSubset   
    elif 'GLAH'in input_file:
        sdsSubset = GLASSubset
    elif 'ATL08' in input_file:
        sdsSubset = ATLAS8Subset
    elif 'ATL03' in input_file:
        sdsSubset = ATLAS3Subset
    else:
        print(f"{input_file} is not a GEDI or ICESat file, or it is not the right product version.")
        sys.exit(2)
    
    # Append additional datasets to extract besides the default ones if provided
    if layerSubset is not None:
        [sdsSubset.append(y) for y in layerSubset]
    
    # Subset to the selected datasets
    file1SDS = [c for c in file1SDS if any(c.endswith(d) for d in sdsSubset)]
        
    # Get unique list of data sets and subset to user-defined subset
    groups = []
    for h in file1SDS:
        group = h.split('/', 1)[0]
        if group not in groups and group in data_group:
            groups.append(group)

    file1DF = pd.DataFrame()  # Create empty dataframe to store the datasets    
    
    # Loop through each data group and create a geodataframe with lat/lon for each shot, then clip to ROI
    if 'GEDI02_A' in input_file or 'GEDI02_B' in input_file:
        for b in groups:
            geoDF=GEDI_group_process(file1Name,b,file1, file1SDS, sdsSubset,ROI)
            file1DF = file1DF.append(geoDF)
            del geoDF
    elif 'GLAH'in input_file:
        for b in groups:
            geoDF=GLAS_group_process(file1Name,b,file1, file1SDS, sdsSubset,ROI)
            file1DF = file1DF.append(geoDF)
            del geoDF
    elif 'ATL08' in input_file:
        for b in groups:
            geoDF=ATLAS08_group_process(file1Name,b,file1, file1SDS, sdsSubset,ROI)
            file1DF = file1DF.append(geoDF)
            del geoDF
    elif 'ATL03' in input_file:
        for b in groups:
            geoDF=ATLAS03_group_process(file1Name,b,file1, file1SDS, sdsSubset,ROI)
            file1DF = file1DF.append(geoDF)
            del geoDF
    else:
        print(f"{input_file} is not a GEDI or ICESat file, or it is not the right product version.")
        sys.exit(2)
    #Check if the resulting data frame is empty.
    if file1DF.shape[0] == 0:
        print(f"{input_file} does not have points intersecting the bounding box of the input ROI.")
        return None
    # Convert to geodataframe and add crs
    file1DF = gp.GeoDataFrame(file1DF)
    file1DF.crs = 'EPSG:4326'
    
    # Subset the output DF to the actual boundary of the input ROI
    outDF = gp.overlay(file1DF, finalClip)
    del outDF[0]            
# --------------------------------Step 6---------------------------------------------- #
#Export the resulting points to a shapefile or geoJSON file. 
# Check for empty output dataframe
    if outDF.shape[0] == 0:
        print(f"{input_file} intersects the bounding box of the input ROI, but no shots intersect final clipped ROI.")
        return None
    try:
        if file_ext=='json':
            output_file=output_directory+'/'+file1Name+'.json'
            outDF.to_file(output_file, driver='GeoJSON')
        else:
            output_file=output_directory+'/'+file1Name+'.shp'
            outDF.to_file(output_file)
    except ValueError:
        print("Make sure the file extension is shp or json")

"""
Worflow functions for parallel process
"""

def download_workflow(q):
    while not q.empty():
        try:
            f, temporary_directory = q.get()
            earthdata_download(f, temporary_directory)
        except ValueError as val_error:
            print(val_error)
        except Exception as error:
            print(error)

def HDF5_subset_workflow(q):
    while not q.empty():
        try:
            ROI_file,input_file,output_directory,data_group = q.get()
            HDF5_data_subset(ROI_file,input_file,output_directory, data_group)
        except ValueError as val_error:
            print(val_error)
        except Exception as error:
            print(error)
        except SystemExit:
            print("the file does not containt points on the ROI")

def per_file_workflow(q):
    while not q.empty():
        try:
            f, temporary_directory, output_directory, ROI_file, data_group = q.get()
            fileName = os.path.join(temporary_directory, f.split('/')[-1].strip())
            #Check if the file has already been downloaded
            if not os.path.isfile(fileName):
                earthdata_download(f, temporary_directory)
            #Check if the file has already been processed
            SHPfileName = os.path.join(output_directory, re.sub('h5', 'shp', f.split('/')[-1].strip(),flags=re.IGNORECASE))
            if not os.path.isfile(SHPfileName): 
                old_stdout=sys.stdout
                message=StringIO()
                sys.stdout=message
                HDF5_data_subset(ROI_file,fileName,output_directory, data_group)
                sys.stdout=old_stdout
                if re.search("no shots intersect", message.getvalue()) or re.search("does not have points intersecting", message.getvalue()):
                    os.remove(fileName)
            else:
                continue
            gc.collect()
        except ValueError as val_error:
            print(val_error)
        except Exception as error:
            print(error)
        except SystemExit:
            print("System Exit occurred")

"""
Functions for parallel workflow execution
"""

def parallel_download(available_files,temporary_directory2, num_workers):            
        q = Queue()
        for f in available_files:
            if not os.path.isfile(os.path.join(temporary_directory2, f.split("/")[-1])):
                q.put((f, temporary_directory2))
        workers = Pool(num_workers, download_workflow,(q,))
        workers.close()
        workers.join()

def parallel_subset(input_files,ROI_file,output_directory2,data_group,num_workers):            
        q = Queue()
        for f in input_files:
            if not os.path.isfile(os.path.join(output_directory2, re.split(".h5", f.split("/")[-1],flags=re.IGNORECASE)[0])+'.shp'):
                q.put((ROI_file,f,output_directory2,data_group))
        workers = Pool(num_workers, HDF5_subset_workflow,(q,))
        workers.close()
        workers.join()
        
def parallel_workflow_per_file(available_files, temporary_directory2, ROI_file, output_directory2, data_group, num_workers):
    q= Queue()
    for f in available_files:
        q.put((f, temporary_directory2, output_directory2, ROI_file, data_group)) 
    workers = Pool(num_workers, per_file_workflow,(q,), maxtasksperchild=5)
    workers.close()
    workers.join()
 