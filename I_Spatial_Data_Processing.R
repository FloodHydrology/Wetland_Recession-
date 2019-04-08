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
library(rgdal)

#Set data dir
data_dir<-'//storage.research.sesync.org/njones-data/Research Projects/Delmarva_Hysteresis/gage_analysis/'
wbt_dir<-     "C:\\WBT/whitebox_tools"
scratch_dir<- "C:\\ScratchWorkspace/"

#download, mask, and combine fac adn fdr rasters
mask<-readOGR(paste0(data_dir,"mask.shp"))
dem_a<-raster(paste0(data_dir, "NHDPlus02/Elev_Unit_a/elev_cm"))
dem_b<-raster(paste0(data_dir, "NHDPlus02/Elev_Unit_b/elev_cm"))
dem_a<-crop(dem_a, extent(mask))
dem_b<-crop(dem_b, extent(mask))
dem<-mosaic(dem_a,dem_b, fun=min)
remove(list=c("dem_a","dem_b"))

#download gage shapefile
mask<-st_as_sf(mask)
gages<-st_read(paste0(data_dir, "NHDPlus02/StreamGageEvent.shp"))
gages<-st_transform(gages, crs=crs(mask))
gages<-st_zm(gages)
gages<-gages[mask,]

#download wetland
fp_wetlands<-st_read(paste0(data_dir, "fp_wetlands.shp"))
nfp_wetlands<-st_read(paste0(data_dir, "nfp_wetlands.shp"))
wetlands<-rbind(fp_wetlands, nfp_wetlands)
remove(fp_wetlands, nfp_wetlands)

#######################################################################
#Delineate Watershed---------------------------------------------------
#######################################################################
#Preprocess DEM--------------------------------------------------------
#Export DEM and stream layer to local working directory
writeRaster(dem, 
            paste0(scratch_dir,"dem.tif"), 
            overwrite=T)

#Gaussian Filter
system(paste(paste(wbt_dir), 
             "-r=GaussianFilter", 
             paste0("--wd=",scratch_dir),
             "-i='dem.tif'", 
             "-o='dem_filter.tif'",
             "--sigma=3"))

#Fill "single cell" depressions
system(paste(paste(wbt_dir),
             "-r=FillSingleCellPits",
             paste0("--wd=",scratch_dir),
             "--dem='dem_filter.tif'",
             "-o='dem_breach_minor.tif'"))

#Breach larger depressions
system(paste(paste(wbt_dir), 
             "-r=BreachDepressions", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_minor.tif'", 
             "-o='dem_breach_major.tif'"))

#Define stream network------------------------------------------------
#Create Flow Accumulation Raster
system(paste(paste(wbt_dir), 
             "-r=D8FlowAccumulation", 
             "--out_type='cells'",
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_major.tif'", 
             "-o='fac.tif'"))

#Define stream network based on fac threshold
flowgrid<-raster(paste0(scratch_dir,"fac.tif"))
flowgrid[flowgrid<1e3]<-NA
flowgrid<-flowgrid*0+1
flowgrid@crs<-dem@crs
writeRaster(flowgrid,paste0(scratch_dir,"flowgrid.tiff"), overwrite=T)

#Delineate wateshed-------------------------------------------------
#Export pour point to scratch directory 
st_write(gages, paste0(scratch_dir,"pnts.shp"), delete_layer = T)

#Run flow direction [note we can skip breaching and/or filling sinks b/c we are using NHD data
system(paste(paste(wbt_dir), 
             "-r=D8Pointer", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_major.tif'", 
             "-o='fdr.tif'",x
             "--out_type=sca"))

#Delineate watershed
system(paste(paste(wbt_dir),
             "-r=Watershed", 
             paste0("--wd=",scratch_dir),
             "--d8_pntr='fdr.tif'", 
             "--pour_pts='pnts.shp'",
             "-o='watershed.tif"))

#Download watershed
ws_grd<-raster(paste0(scratch_dir,"watershed.tif"))
ws_shp<-rasterToPolygons(ws_grd, dissolve = T)
