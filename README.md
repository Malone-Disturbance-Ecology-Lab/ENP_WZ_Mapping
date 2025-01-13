# ENP_WZ_Mapping
# Author: Boya Zhang, Sparkle Malone

**Project Overview**: 
This project uses Landsat optical remote sensing observations to map the land-cover changes of the South Florida Everglades coastal wetlands since 1984. Land covers were classified into "Vegetation Zone" (VZ, mostly mangroves), "White Zone" (WZ, dead aboveground biomass), and "Transition Zone" (TZ, with optical signature between the VZ and WZ). About 1500 point samples were manually classified for eight years out of the 40 years time span (with ~6 years interval) using available high-resolution Google Earth imagery. We used an Machine Learning (ML) approach to spatially and temporally extrapolate the classification to the regional scale for the past 40 years. We conducted spatial modelling for the classification time series products to track the trend of white zone expansion and mangrove migration.   

**Input data**:
1. Manually classified point samples of eight years (1986, 1994, 1999, 2003, 2009, 2013, 2017, and 2023) in the format of shapefile
2. Optical observations were obtained from Google Earth Engine (GEE), and the script was available in this repo


**File Structure**:
This repository includes both Google Earth Engine (GEE) and Python jupyter notebook code.
GEE code collects time series of optical observations. The Python jupyter notebook code compiles the optical data and manually labelled vegetation classification data into a machine learning framework to predict regional vegetation classification since 1984. 

To reproduce the results, one need to follow the four-step procedure. The order of steps were described as below.


GEE code (Step 2):
1. LS_timeseris. Extract the time series optical data given a lat and lon values
2. point_loop. Extract time series optical data given mulitple pairs of lat and lon values

Python Jupyter Notebook:

1. Data_Prep.ipynb (Step 1) processes kml files and converting into shapefiles that will be used for Google earth engine (GEE).
2. ML_Model.ipynb (Step 3). This script needs two type of data inputs: (1) the data file generated from GEE processing including optical spectral values of all the bands plus NDVI, NDMI, and SAVI; (2) point shapefile with classification information. The script conducts temporal and spatial transfer learning to map the classification map time series.
3. Spatial modelling (Step 4). This script needs the inputs of vegetation classification time series to model the changes of spatial pattern of white zone over the past 40 years. 
