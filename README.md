# Global deltas #
Code for GEDI, ICESat GLAS and ICESat2 ATLAS data processing, for the Global Delta project.

This code defines functions for creating a workflow for querying, downloading, and subsetting GEDI and ICESat data
This code was based on publicly available code and resources.

## Original code sources:
- ICESat query and download: Download script from the NSIDC webpage. Example: https://nsidc.org/data/GLAH06/versions/34 in 'download data' click in 'Download script'.
- GEDI query: GEDI finder web service https://lpdaacsvc.cr.usgs.gov/services/gedifinder
- GEDI subsetting: GEDI_Subsetter.py code by Cole Krehbiel https://git.earthdata.nasa.gov/projects/LPDUR/repos/gedi-subsetter/browse/GEDI_Subsetter.py 
- Earth Data login setup code EarthdataLoginSetup.py created by Cole Krehbiel https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_python/browse/EarthdataLoginSetup.py

## Datasets information
### GEDI L2 data
GEDI L2 Elevation and Height Metrics Data short name is GEDI02_A and GEDI02_B. there are 2 versions available 001 and 002
- The data_group for GEDI are any of the BEAM names ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
- GEDI L2A files have 765 science dataset layers per beam. The default layer list only includes 26 layers. 
- GEDI L2B files have 251 science dataset layers per beam. The default layer list only includes 38 layers. 
- The approximate pixel size for this GEDI data is 25 m. 

### ICESat GLAS data
GLAS/ICEsat L2 Global Land Surface Altimetry Data short name is GLAH14, the latest version for this product is version 034
- The data_group for GLAS are the HZ groups ['Data_1HZ', 'Data_40HZ']
- The approximate pixel size for this ICESat data is 65 m.


### ICESat2 ATLAS L8 data
ICESat2 ATLAS L3A Land and Vegetation Height short name is ATL08, the latest version for this product is version 004.
- The data_group for ATLAS are any of the ground tracks names ['gt1l', 'gt1r', 'gt2l','gt2r','gt3l','gt3r']
- The approximate pixel size for this ICESat2 data is 100 m

