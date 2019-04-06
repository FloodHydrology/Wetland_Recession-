#######################################################################
#Title: Spatial Data Preperation
#Coder: Nate Jones
#Date: 4/5/2019
#Description: Prep spatial data on gages accross Delmarva
#######################################################################

#######################################################################
#Setup workspace-------------------------------------------------------
#######################################################################
#Clear workspace
remove(list=ls())

#load appropriate packages
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(raster)
library(sf)

#Set data dir
data_dir<-'/nfs/njones-data/Research Projects/Delmarva_Hysteresis/gage_analysis/'

#download, mask, and combine fac adn fdr rasters
fac_a<-raster(paste0(data_dir, "NHDPlus02/FAC_FDR_Unit_a/fac"))
fac_b<-raster(paste0(data_dir, "NHDPlus02/FAC_FDR_Unit_b/fac"))
fdr_a<-raster(paste0(data_dir, "NHDPlus02/FAC_FDR_Unit_a/fac"))
fdr_b<-raster(paste0(data_dir, "NHDPlus02/FAC_FDR_Unit_a/fac"))
mask<-st_read(paste0(data_dir,"mask.shp"))
fac_a<-crop(fac_a, mask)
fac_b<-crop(fac_b, mask)
fdr_a<-crop(fdr_a, mask)
fdr_b<-crop(fdr_b, mask)
fac<-mosaic(fac_a,fac_b, fun=min)
fdr<-mosaic(fdr_a,fdr_b, fun=min)
remove(list=c("fac_a","fac_b", "fdr_a","fdr_b"))

#download gage shapefile
gages<-st_read(paste0(data_dir, "NHDPlus02/StreamGageEvent.shp"))
gages<-st_transform(gages, crs=crs(mask))
gages<-st_zm(gages)
gages<-gages[mask,]

#######################################################################
#Delineate Watershed---------------------------------------------------
#######################################################################



