This repository contains the basic Google Earth Engine and Matlab scripts for the automated mapping of avalanches in Sentinel-1 images.

The preprocessing of the Sentinel-1 RGB images is done in Google Earth Engine (S1_preprocessing_GEE.txt) and the mapping of avalanches is done with the AvalancheMapping_Hispar.m script which gives an example for the Hispar region. This script calls the function geotiffcrop_shp.m.

The shifting of the avalanche outlines caused by glacier elevation shange is done with the matlab function Shift_outlines_Hispar_dh.m
