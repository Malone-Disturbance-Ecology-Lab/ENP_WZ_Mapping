# ENP_WZ_Mapping
# Author: Boya Zhang, Sparkle Malone

All files are jupyter notebook Python code
Data_Prep.ipynb is the script processing kml files and converting into shapefiles that will be used for Google earth engine (GEE).
ML_Model.ipynb processes the data file as a result of GEE processing and the file includes the optical spectral values of all the bands plus NDVI, NDMI, and SAVI. It also read in the label informatio from the shapefile for model training and testing.
