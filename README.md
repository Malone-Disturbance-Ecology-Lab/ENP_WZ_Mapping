# ENP_WZ_Mapping
# Author: Boya Zhang, Sparkle Malone

# Part 1 Remote Sensing 
**Project Overview**: 
This project uses Landsat optical remote sensing observations to map the land-cover changes of the Southeast Saline Everglades (SESE) coastal wetlands between 1984 to 2024. Land covers were classified into "Vegetation Zone" (VZ, mangroves), "White Zone" (WZ, dead aboveground biomass with exposed marl), and "Transition Zone" (TZ, with optical signature between the VZ and WZ). 

About 1500 point samples were manually classified for each of the eight years out of the 41 years time span (with ~6 years interval) using available high-resolution Google Earth basemap imagery. We used an Machine Learning (ML) approach to spatially and temporally for classification for the SESE for the past 40 years. 

The results land-cover maps will be used for spatial modelling (next part) to track the trend of white zone expansion and mangrove encroachment.   

Note: The scripts pull data using absolute directory on the server. When a script calls another script/function, it uses relative directory. This way allows users to directly run the code.

**Input data**:
1. Manually classified point samples of eight years (1986, 1994, 1999, 2003, 2009, 2013, 2017, and 2023) in the format of shapefile
2. Optical observations were obtained from Google Earth Engine (GEE), and the script was available in this repo

These data are stored on the server: /Volumes/Malonelab/Research/Whitezone/Data/KML/


**Script Structure**:
This repository includes both Google Earth Engine (GEE) and Python jupyter notebook code.
GEE code collects time series of optical observations. The Python jupyter notebook code integrates the optical data and label data into a machine learning framework to predict regional vegetation classification since 1984 at an annual basis. 

To reproduce the results, one need to follow the six-step procedure. The order of steps were described as below.


Step 1 Generating shapefile for GEE to extract LS values
Step 2 GEE code to extract annual median values of LS data
Step 3 Sensitivity analysis for selecting the best VZ threshold
Step 4 Using the GE years to test spatial and temporal transfer learning
Step 5 Emulation process to generate land cover maps for the past 41 years
Step 6 Evaluate products for selecting optimal change detection threshold & Plot time series maps

The final output of time series raster (Change detection threshold 0f 0.6 with Closest GE years) is output for spatial modelling.

# Part 2 Spatial Modelling

Part 2 uses R scripts to define transect line, and get land cover info by intersecting with the 41 years of classification maps. Correlation between WZ proportion and salinity/wl was calculated at the transect level to infer driver mechanisms. 

The entire work flow (Spatial_Work_Flow.R) was divided into three sub work flows (Spatial_Analysis_step*.R):
Step 1 (Spatial_Analysis_step1.R) makes buffer of transects to intersect with center points of pixels, generates the changes of leading and trailing edges between 1984 and 2024 to calculate WZ expansion/retraction and rates. 

Step 2 (Spatial_Analysis_step2.R) loads the proportion of each land-cover class in time for all the transects (resulting from Step 1), then make figure plots, including time series of total area changes (Figure 7 in the Manuscript). Step 2 reads all the salinity and wl files to compile them as annual average and then calculate transect-based correlation between WZ proportion and salinity/wl respectively, and make interactive plots in R studio. 

Step 3 (Spatial_Analysis_step3.R) uses the correlation results from step 2 to model driver mechanism of WZ expansion.

**Input data**:
Raster files of land cover classification time series produced by Step 1
Boundary lines, including coast line and inland line (i.e., trailing and leading edges). 

**File Structure**:
Flow_Transect_Analysis.R is the work flow that calls all the functions
Functions:
  1. function_MakeVegetationZoneLayers.R polygonize the rasters
  2. function_Make_TransectLayers.R  make transects
  3. Spatial_Work_Flow.R analyze the transects and to produce plots (Spatial_Work_Flow is the entire flow, which is splitted into three sub scripts
      including Spatial_Analysis_step1.R, Spatial_Analysis_step2.R, Spatial_Analysis_step3.R
  4. function_called_by_main_script.R includes the functions called by Spatial_Work_Flow.R for spatial analysis or generating plots

