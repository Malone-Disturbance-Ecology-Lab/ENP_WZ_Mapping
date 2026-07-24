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

**Main Workflow**
Step 1 (`01_Generating_pixels_for_GEE.ipynb`): Generating shapefile for GEE to extract LS values.
Step 2 (`02_GEE_multi_year_median_extraction.rtf`, `02_GEE_single_year_median_extraction.rtf`): GEE code to extract annual median values of LS data.
Step 3 (`03_Sensitivity_analysis_for_VZ_threshold.ipynb`): Sensitivity analysis for selecting the best VZ threshold.
Step 4 (`04_Spatial_Transfer_for_GE_year_maps.ipynb`, `04_Temporal_Transfer_for_GE_years_and_Verification.ipynb`): Using the GE years to test spatial and temporal transfer learning.
Step 5 (`05_Roy_Temporal&Spatial_Transferfor_non-GE_years_L5_1984-2003.ipynb`, `05_Roy_Temporal&Spatial_Transferfor_non-GE_years_L7&8_2004-2023.ipynb`): Emulation process to generate land cover maps for the past 41 years.
Step 6 (`06_Evaluation_with_Independent_Point_Samples.ipynb`, `06_Illustration_of_Time_Series.ipynb`): Evaluate products for selecting optimal change detection threshold & Plot time series maps.

**Producing Figures and Tables **
Figure 1 (`Fig1_Generating_time_line_for_methodology.ipynb`, `Figure1_Conceputal_Framework_plots.ipynb`)
Figure 3 (`Fig3_Spatial_analysis_sa_subregions.ipynb`, `Fig3_Spatial_analysis_wl_subregions.ipynb`, `Fig3_time_trend_of_env_variables.ipynb`)
Figure 8 (`Fig8_Projection_with_rate.ipynb`)
Figure S1 (`FigS1_CERP_analysis.ipynb`)
Figure S3 (`FigS3_Validation_with_sawgrass_proof_of_concept.ipynb`)
Figure S5 (`FigS5_Validation_MODIS_LAI.ipynb`)
Figure S6 (`FigS6_Landsat_signature_violin.ipynb`)
Note S4 (`Response_Directional_Analysis_of_VZ.ipynb`)
Table S3 (`TableS3_cross_validation&temporal_transfer_accuracy.ipynb`)

The final output of time series rasters are output for spatial modelling. The rasters have been published on Environmental Data Initiative (EDI). https://doi.org/10.6073/pasta/ab0c9bb8d372251785aeb3e55a62377e.

# Part 2 Spatial Modelling

Part 2 uses R scripts to define transect line, and get land cover info by intersecting with the 41 years of classification maps. Correlation between WZ proportion and salinity/wl was calculated at the transect level to infer driver mechanisms. 

**Input data**:
Raster files of land cover classification time series produced by Step 1
Boundary lines, including coast line and inland line (i.e., trailing and leading edges). 

**Function Scripts**:

Functions:
  1. `function_MakeVegetationZoneLayers.R`
`polygonize_SpatRast()`: polygonize the rasters 
`Trans_Intersection()`: calculates percentage of intersection over transect

  2. `function_Make_TransectLayers.R`
`makeTransect()`: produce transects that are perpendicular to two reference lines
 
  3. `function_called_by_main_script.R`
includes the functions called by Spatial_Work_Flow_cleanup.R for spatial analysis or generating plots

**Main Workflow**:

The entire work flow can be run by `Spatial_Work_Flow_cleanup.R`. This master script:
Makes buffer of transects to intersect with center points of pixels, generates the changes of leading and trailing edges between 1984 and 2024 to calculate WZ expansion/retraction and rates. 
Loads the proportion of each land-cover class in time for all the transects (resulting from Step 1), then make figure plots, including time series of total area changes (Figure 7 in the Manuscript). Step 2 reads all the salinity and wl files to compile them as annual average and then calculate transect-based correlation between WZ proportion and salinity/wl respectively, and make interactive plots in R studio. 
Uses the correlation results from Step 2 to model driver mechanism of WZ expansion.

