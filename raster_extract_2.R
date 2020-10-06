
install.packages("AzureStor")
install.packages("raster")
install.packages("sf")
install.packages("rgdal")
install.packages("rgeos")
# load shit 

library(tidyverse)
library(data.table)
library(raster)
library(sf)
library(AzureStor)
library(rgdal)
library(rgeos)

setwd("/home/Team8_Vera/")



#------ 
#container <- "https://hlssa.blob.core.windows.net/hls"
# substring(matches_df$date,1,4)
# lat = 46.6101 
# lon = -120.2015 # yakima, WA 10GFQ 
#
# year    = '2019'
# daynum  = '001'   # 1-indexed day-of-year
# folder  = 'S309'  # 'S309' for Sentinel, 'L309' for Landsat
# product = 'S30'   # 'S30' for Sentinel, 'L30' for Landsat
# band    = '01'
# 
# # we need a test of band existance - error here for yakima band
# tile_id = '10GFQ' # See hls.gsfc.nasa.gov/wp-content/uploads/2016/10/S2_TilingSystem2-1.txt
# # this list does not have dates. We need a function that pulls every image of this region. I did not rewrite from python. 
# version = 'v1.4'  # Currently always v1.4
# 
# blob_name = paste0(folder, '/HLS.',product,'.T',tile_id,'.',
#                    year,daynum,'.',version,'_',band,'.tif')
# 
# test_it = paste0(container,"/",folder, '/HLS.',product,'.T',"10GFQ")
# 
# 
# url <- paste0(container,"/",blob_name)
# 
# url_band2 <- gsub(pattern = "*01.tif", "02.tif", url)
# url_band4 <- gsub(pattern = "*01.tif", "04.tif", url)
# url_band9 <- gsub(pattern = "*01.tif", "09.tif", url)
# # 
# 
# # download needed files bands 2 4 8 (others )
# band2_name <-paste0(product,tile_id,".",year, daynum,"_02.tif")
# band4_name <-paste0(product,tile_id,".",year, daynum,"_04.tif")
# band9_name <-paste0(product,tile_id,".",year, daynum,"_09.tif")

#### loading rasters from Azure, calculating derivatives, and building a stack ---------------- 


# work out matches list --------------------------
matches <- fread("Matches.csv")[,2]  ############### this path will need to change
matches <- matches[2:nrow(matches),]
matches_df <- matches %>% 
  as.character() %>%
  strsplit(.,split="[.|,]")%>%
  unlist()

matches_df <- cbind.data.frame(split(matches_df, rep(1:7, times=length(matches_df)/7)))
names(matches_df) <- c("folder", "product","tile_id", "datestring", "version","band","file")

tile_ids <- unique(matches_df$tile_id)

matches_df <- matches_df %>%
  mutate(year = substring(datestring,1,4), 
         doy = substring(datestring,5,8),
         date = format(as.Date(datestring, format='%Y%j'),format ='%m/%d/%Y'), 
         fullnames = matches)

matches_04 <- matches_df$fullnames[month(as.Date(matches_df$date,format ='%m/%d/%Y'))==4]
matches_07 <- matches_df$fullnames[month(as.Date(matches_df$date,format ='%m/%d/%Y'))==7]
matches_exact <- matches_df$fullnames[matches_df$date == "07/02/2017"]

# pull in data ---------------------------
container <- "https://hlssa.blob.core.windows.net/hls"
i=12
for (i in seq(1:nrow(matches_exact))){
  download_from_url(paste0(container,"/",matches_exact[i]), 
                     paste0("/home/Team8_Vera/",matches_exact[i]), overwrite=T)  ## change
}


# build a raster stack ----------------------------------
mypath <- getwd()  ###changethis 
raster_files <- list.files(mypath,full.names = T, pattern = "*.tif") #use pattern = '.tif$' or something else if you have multiple files in this folder

r_name <- list.files(mypath,full.names = F, pattern = "*.tif")

rList <- list() # to save raster values

L <- raster(raster_files[1])
r <- stack(L[[1]])
for(i in 2:length(raster_files)){
  temp <- raster(raster_files[i])
    r <- addLayer(r, temp)
}

#### derivatives function ----------------------------

dervies <- function(rasterstack){## add the indexing (PLEASE someone review)
  #test number of bands
  EVI <- (2*rasterstack[[9]] - rasterstack[[4]])/(rasterstack[[9]] + 6*rasterstack[[4]] -7.5*rasterstack[[2]]+1) 
  NVDI <- (rasterstack[[9]]-rasterstack[[4]])/(rasterstack[[9]]+rasterstack[[4]])
  Albedo <- 0.356*rasterstack[[02]]+0.130*rasterstack[[04]]+0.373*rasterstack[[08]]+0.085*rasterstack[[10]]+0.072*rasterstack[[11]]-0.0018
  GNDVI <- (rasterstack[[09]]-rasterstack[[03]])/(rasterstack[[09]]+rasterstack[[03]])
  MSR <- (rasterstack[[09]] / rasterstack[[04]] - 1) / (sqrt(rasterstack[[09]]/rasterstack[[04]])+1)
  OSAVI <- (rasterstack[[09]]**2-rasterstack[[04]])*1.5*rasterstack[[09]]**2+rasterstack[[04]]+0.5
  PVR <- (rasterstack[[03]]-rasterstack[[04]])/(rasterstack[[03]]+rasterstack[[04]])
  VARIG <- (rasterstack[[03]] - rasterstack[[04]])/(rasterstack[[03]]+rasterstack[[04]]-rasterstack[[02]])
  SIPI <- (rasterstack[[08]]-rasterstack[[02]])/(rasterstack[[08]]+rasterstack[[02]])
  GARI <- (rasterstack[[08]]-rasterstack[[03]]+rasterstack[[02]]-rasterstack[[04]])/(rasterstack[[08]]+rasterstack[[03]]-rasterstack[[02]]+rasterstack[[04]])
  mNVDI <- (rasterstack[[08]]-rasterstack[[04]])/(rasterstack[[08]]+rasterstack[[04]]-2*rasterstack[[02]])
  EPIChlb <- (0.0337*(rasterstack[[04]]/rasterstack[[03]])**1.8695)
  NDII <- (rasterstack[[08]] - rasterstack[[10]]) / (rasterstack[[08]] + rasterstack[[10]])
  restack <- stack(EVI, NVDI,Albedo,GNDVI,MSR,OSAVI,PVR,VARIG, SIPI,GARI,mNVDI, EPIChlb, NDII)
  return(restack)
}

#### build the derivatives step for each tile ----------------------------

r2 <- dervies(r)
names(r2) <- c("EVI", "NVDI","Albedo","GNDVI","MSR","OSAVI","PVR",
               "VARIG", "SIPI","GARI","mNVDI", "EPIChlb", "NDII")

## above here is run once. 

#load shape files and data and reproject to UTM 

shapes <- shapefile("tGrant_Yakima.shp")
extent <- extent(r2)
shapes_utm_all <- spTransform(shapes, "+proj=utm +zone=10 +ellps=WGS84 +units=m +no_defs")

test <- raster::intersect(shapes_utm_all, extent)
shapes_utm_all_crop <- crop(shapes_utm_all, r2[[1]])
shape_df <- shapes@data %>%
  mutate(date = as.Date(LstSrvD),
         year = format(as.Date(LstSrvD), '%Y'))


### duplicate this and extract
shape_df_sub <- subset(shape_df,year=="2017")
shapes_sub <- shapes[shapes$ID %in% shape_df_sub$ID,]
#check this it
shapes_utm <- spTransform(shapes_sub, "+proj=utm +zone=10 +ellps=WGS84 +units=m +no_defs")

# shape_df_sub2 <- subset(shape_df,year=="2017")[51:250,]
# shapes_sub2 <- shapes[shapes$ID %in% shape_df_sub2$ID,]
# #check this it
# shapes_utm2 <- spTransform(shapes_sub2, "+proj=utm +zone=10 +ellps=WGS84 +units=m +no_defs")

# shape_df_sub3 <- subset(shape_df,year=="2017")[251:2500,]
# shapes_utm3<- shapes_utm_all[shapes_utm_all$ID %in% shape_df_sub3$ID,]
# 
# plot(r2[[1]])
# lines(shapes_utm_all)
# 
# plot(shapes_utm_all)
# lines(extent) 
# 
# plot(shapes_utm3)
# lines(extent)

### extraction bit---------------------------

pathout <- "points_data.csv"    ### Path out

bands_out_func <- function (band_stack, field_shapes, pathout){
    bands_out <-list()
    for(i in seq(1:dim(band_stack)[3])){
      bands_out[[i]] <- raster::extract(band_stack[[i]], field_shapes)
      }
    bands_out <- na.exclude(bands_out)
    bands_out_sum <- lapply(bands_out, function(x) sapply(x, mean))
    bands_out_sum <- do.call(rbind, bands_out_sum) %>% t() %>% as.data.frame()
    names(bands_out_sum) <- names(band_stack)
    shape_df_out <- cbind(field_shapes@data, bands_out_sum)
    fwrite(shape_df_out, pathout)
}
