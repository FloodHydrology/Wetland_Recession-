#######################################################################
#Title: Watershed Delineation
#Coder: Nate Jones
#Date: 4/5/2019
#Description: Delineate watershed for each gage and estimate wetland area
#######################################################################

#######################################################################
#Setup workspace-------------------------------------------------------
#######################################################################
#Clear workspace
remove(list=ls())

#load appropriate packages
library(tidyverse)
library(raster)
library(sf)
library(rgdal)
library(stars)
library(parallel)
library(rslurm)   

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
wetlands<-st_transform(wetlands, crs=paste(dem@crs))
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

#Read fac and define values
fac<-raster(paste0(scratch_dir,"fac.tif"))
fac<-c(cellStats(fac, min), cellStats(fac, max))

#Define stream network based on fac threshold
system(paste(paste(wbt_dir), 
             "-r=Reclass", 
             paste0("--wd=",scratch_dir),
             "-i='fac.tif'",
             "-o='flowgrid.tif",
             paste0("--reclass_vals='0;",fac[1],";1000;1;1000;",fac[2]+1)
))

#Create Pour Points--------------------------------------------------
#Create UID 
gages$id<-seq(0, nrow(gages)-1)

#Export pour point to scratch directory 
st_write(gages, paste0(scratch_dir,"pnts.shp"), delete_layer = T)

#Create pour pnt raster
system(paste(paste(wbt_dir), 
             "-r=VectorPointsToRaster", 
             paste0("--wd=",scratch_dir),
             "-i='pnts.shp'", 
             "--field=UID",
             "-o=pp.tif",
             "--assign=min",
             "--nodata",
             "--base=dem.tif"))

#Jenson Snap Pour point
system(paste(paste(wbt_dir),
             "-r=JensonSnapPourPoints", 
             paste0("--wd=",scratch_dir),
             "--pour_pts='pp.tif'", 
             "--streams='flowgrid.tif'",
             "-o='pp_snap.tif",
             "--snap_dist=1000"))

#Convert back to point file
snapgrid<-raster(paste0(scratch_dir,"pp_snap.tif"))
snapvector<-getValues(snapgrid)
snapvalues<-na.omit(snapvector)
snappnts<- tibble(snap_length=which(snapvector %in% snapvalues)) %>% 
  mutate(x = (snap_length %% ncol(snapgrid))*res(snapgrid)[1]+extent(snapgrid)[1]-res(snapgrid)[1]/2, 
         y = extent(snapgrid)[4]-((ceiling(snap_length/ncol(snapgrid)))*res(snapgrid)[2])+res(snapgrid)[2]/2,
         id= snapvalues)

#Create sf file and export to workspace
snappnts<-st_as_sf(snappnts, 
                   coords=c("x","y"), 
                   crs=paste(dem@crs))
st_write(snappnts, paste0(scratch_dir,"snap.shp"), delete_layer = T)

#Delineate wateshed-------------------------------------------------
#Run flow direction 
system(paste(paste(wbt_dir), 
             "-r=D8Pointer", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_major.tif'", 
             "-o='fdr.tif'",
             "--out_type=sca"))
#Delineate watershed
system(paste(paste(wbt_dir),
             "-r=unnest_basins", 
             paste0("--wd=",scratch_dir),
             "--d8_pntr='fdr.tif'", 
             "--pour_pts='snap.shp'",
             "-o='watershed.tif"))

#Clean up workspace
remove(list=ls()[ls()!='gages'    &
                 ls()!='wetlands' &
                 ls()!='data_dir' &
                 ls()!='wbt_dir'  &
                 ls()!='scratch_dir'])

#######################################################################
#Estimate metrics for each watershed-----------------------------------
#######################################################################
#Download PP file 
pp<-raster(paste0(scratch_dir,"pp_snap.tif"))

#Locate files with watershed shapes
files<-list.files(scratch_dir)
files<-files[substr(files, 1,9)=="watershed"]
index_fun<-function(i){
  temp_watershed<-raster(paste0(scratch_dir, files[i]))
  data.frame(wuid = unique(temp_watershed), 
             file = files[i])
}
watershed_index<-lapply(X = seq(1, length(files)), FUN=index_fun)
watershed_index<-do.call(rbind, watershed_index)

#Create function to estimate wetland and watershed area
attributes_fun<-function(j){

  #Download watershed file of interest
  watershed<-raster(paste0(scratch_dir, watershed_index$file[j]))
  
  #Identify watershed [use watershed 1 for now]
  w_grd<-watershed %in% watershed_index$wuid[j]
  w_grd[w_grd==0]<-NA
  
  #crop raters to reasonable extent
  w_pnts <- tibble(w_length=which(w_grd@data@values)) %>% 
    mutate(x = (w_length %% ncol(w_grd))*res(w_grd)[1]+extent(w_grd)[1]-res(w_grd)[1]/2, 
           y = extent(w_grd)[4]-((ceiling(w_length/ncol(w_grd)))*res(w_grd)[2])+res(w_grd)[2]/2)
  w_grd<-crop(w_grd, c(min(w_pnts$x,na.rm = T)-res(w_grd)[1]*10,
                       max(w_pnts$x,na.rm = T)+res(w_grd)[1]*10, 
                       min(w_pnts$y,na.rm = F)-res(w_grd)[2]*10, 
                       max(w_pnts$y,na.rm = F)+res(w_grd)[2]*10))
  
  #convert to polygon
  w_shp<- w_grd %>% st_as_stars() %>% st_as_sf(., as_points = FALSE, merge = TRUE) #rasterToPolygons(w_grd, dissolve = T)
  w_shp<-st_as_sf(w_shp)
  st_crs(w_shp)<-st_crs(wetlands)
  
  #Estimate number of wetlands
  wet_shp<-wetlands[w_shp,]
  
  #Identify pour point of interest 
  pp_grd<-crop(pp, c(min(w_pnts$x,na.rm = T)-res(pp)[1]*10,
                     max(w_pnts$x,na.rm = T)+res(pp)[1]*10, 
                     min(w_pnts$y,na.rm = F)-res(pp)[2]*10, 
                     max(w_pnts$y,na.rm = F)+res(pp)[2]*10))
  pp_shp<-data.frame(rasterToPoints(pp_grd))
  pp_shp<-st_as_sf(pp_shp, 
                   coords=c("x","y"), 
                   crs=st_crs(wetlands))
  w_line<-st_cast(w_shp, "LINESTRING")
  pp_shp$dist<-st_distance(pp_shp, w_line, by_element = T)
  pp_shp<-pp_shp[pp_shp$dist==min(pp_shp$dist, na.rm=T),]
  
  #Export watershed and wetland areas (ha)
  c(pp_shp$pp_snap, as.numeric(st_area(w_shp)/10000), sum(wet_shp$hectares))
}

#Execute function
t0<-Sys.time()
n.cores<-detectCores() #detect number of cores
cl <- makePSOCKcluster(n.cores) #Create Clusters
clusterEvalQ(cl, library(raster))  #Send clusters the libraries used
clusterEvalQ(cl, library(sf))  #Send clusters the libraries used
clusterEvalQ(cl, library(tidyverse))  #Send clusters the libraries used
clusterEvalQ(cl, library(stars))  #Send clusters the libraries used
clusterExport(cl, c('wetlands', 'pp','watershed_index', 'scratch_dir'), env=environment())  #Send Clusters function with the execute function
x<-parLapply(cl, seq(1,nrow(watershed_index)), attributes_fun) #Run execute Function
stopCluster(cl)  #Turn clusters off
tf<-Sys.time()
tf-t0

#Unlist
output<-data.frame(do.call(rbind,x))
colnames(output)<-c("PP_ID", "Watershed_Area", "Wetland_Area")

#Join with gages
gages$PP_ID<-seq(1,nrow(gages))
gages<-left_join(gages, output, by='PP_ID')
gages<-gages[,c("SOURCE_FEA", "Watershed_Area", "Wetland_Area")]

#Export data
write.csv(gages,"output/wetland_area.csv")
