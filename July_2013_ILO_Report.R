################################################################################
# Vegetation Index / Rainfall Estimate Validation Exercises 
# Created by: Anthony Louis D'Agostino (ald2187@columbia.edu)
# Date Created: June 18, 2013  
################################################################################


# Edits: 

# ======== PURPOSE =============================================================
# 
# to emulate results from sniiddemo (around line 1020)

# variables of interest:
  # vi.mat - creates matrix of agreement percentage with ARC data by VI product 
  # arc.ranks - used for ensuring results match earlier work - gives ECDF ranks across early and late phases 
  # 



# ======= Requirements ========

###  sitetable.csv - copied from Dry2012Satellite's common.data folder under the sniidharita repository, current as of June 18, 2013.  Will identify any discrepancies between this script's output and the earlier ILO reporting data   
###  rainfall files from sniidharita script 

# ================================

# ======== NOTES ===============================================================
# Datasets may be applicable only for some regions, check beforehand
# Need to determine which rainfall products we are comparing against before pushing forward with design
# ==============================================================================

# ======== TASKS ===============================================================
# * Variable lag length for VI data against rainfall (creating a vector with multiple lags enables us to do immediate comparisons - this will be in terms of the 16-day composites)  Also consider taking the average of multiple lags (e.g., 1st 16-day composite lag behind, and 2nd)
# * Masking out pixels with spatial correlation values below threshold
# * Changing aggregate pixel size (may be based on some criteria) - currently all pixels within grid are evenly weighted
# *               
# ==============================================================================

# Next Steps: 
#   * Determine how best to save the pixel-level data in a meaningful way 
#Apply some weighting criteria to grid based on distance to village pixel 
# Determine appropriate graphing facilities for exporting data  
# User specifies which products interested in, use "if" statements to match product with URL assembly

#???? Uses some of the core existing functions written for sniidharita 



# ---- Initialization and GPS Coordinate Acquisition -----

rm(list = ls())
library(zoo)
library(ggmap)  # for visualizations
library(ggplot2)
library(reshape2)
# require(abind)   # for reshaping matrices into arrays 
 
# =============================================================================#
# Which folder is the script in and any generated output files? 
    base.path   <-  "~/Documents/Index_Insurance/ILO/"

      # load auxiliary functions like "vi.compare" 
      source(paste0(base.path,"vi_functions.R"))
# =============================================================================#





updaterainfall <- FALSE 


#########################################
# ===== READ IN PRECIPITATION DATA ======
#########################################


#range of years covered in analysis, specified in Phase One report 
  years <- as.character(seq(2003,2011))  
  
  # Choose the scenario to work with (draws upon existing files in that folder)
    setwd("~/Documents/Index_Insurance/sniidharita_fist/WIIET_Data/Dry2013Initial")

#precip.path <- "~/Documents/Index_Insurance/sniidharita_fist/WIIET_Data/Dry2013Initial" 
  site.data <- read.csv(paste0("common.data/sitetable.csv"), header = TRUE, row.names = 1)
  dekad.cal <- read.csv(paste0("common.data/CalendarYearDayDecad.csv"), header = F, skip = 2, row.names = 1)
  scen <- "scenarios/"
  rbf.fn <- "/EthRFEadj/rainbyphase.csv"

 



# establish parameters for comparison table 
  prods <- c("EVI", "NDVI", "NDWI")  
  satellite <- "Mod"  
  phase <- "Early" 
  badyear.thres <- 1/6

# create matrices for early, late comparison of arc estimates
  arc.ranks <- array(data = NA, dim = c(nrow = dim(site.data)[1], length(years), 2))
  dimnames(arc.ranks) <- list(rownames(site.data), years, c("Early", "Late"))
  
# create a matrix where comparison data across products will be kept 
  vi.mat <- matrix(ncol = length(prods), nrow = length(rownames(site.data)))
  rownames(vi.mat) <- rownames(site.data)
  colnames(vi.mat) <- prods

# create an array with rankings for each product, each site, across all years
  vi.rankings <- array(data = NA, dim = c(length(rownames(site.data)), length(years), length(prods)))
  dimnames(vi.rankings) <- list(rownames(site.data), years, prods)

for (site in rownames(site.data)){
  
  print(paste0("Now calculating comparisons for ",site))
  
    rfe <- read.csv(paste0(scen,site,rbf.fn), header = T, row.names = 1)
    rfe.sub <- rfe[years,]  #only years we're interested in 
    rfe.p1 <- rfe.sub[years,"X1"] # first phase 
    rfe.p2 <- rfe.sub[years,"X3"] # second phase 
  # How to use ecdf on a matrix, rather than splitting into vectors?  
  
    Fn1 <- ecdf(rfe.p1)
    Fn2 <- ecdf(rfe.p2) 
    r1 <- round(Fn1(rfe.p1), digits = 2)   # ranks for 1st phase
    r2 <- round(Fn2(rfe.p2), digits = 2)   # ranks for 2nd phase
  
  # identify the years in the bottom "badyear.thres" - problem with ties (what to do?)
    rfe.years.early <- rownames(rfe.sub)[which(rfe.p1 <= quantile(rfe.p1,probs = badyear.thres))]   
    rfe.years.late <- rownames(rfe.sub)[which(rfe.p2 <= quantile(rfe.p2,probs = badyear.thres))]
  
  # define sowing windows to calculate ARC rainfall by using dekad calendar 
  # days given for comparison against earlier reports 
    early.start <- min(which(dekad.cal == site.data[site,"EarlyFirst"]))
    early.end <- max(which(dekad.cal == site.data[site,"EarlyLast"])) 
    late.start <- min(which(dekad.cal == site.data[site,"LateFirst"]))
    late.end <- max(which(dekad.cal == site.data[site,"LateLast"]))
    early.dates <- c(rownames(dekad.cal)[early.start], rownames(dekad.cal)[early.end])
    late.dates <- c(rownames(dekad.cal)[late.start], rownames(dekad.cal)[late.end])
    
  # read ARC precip file in day-year format
    arc <- read.csv(paste0(scen,site,"/met.data/precip.daily.csv"), header = T, row.names = 1)
    colnames(arc) <- substr(colnames(arc),2,5)   # remove X from character string preceding colnames
  
  # create subset of only arc data over years in "years", take sum over window timings 
    arc.sub <- arc[,years]
    earlies <- colSums(arc.sub[early.start:(early.end),])
    lates <- colSums(arc.sub[late.start:(late.end),]) 
  
  # identify which years satisfy the quantile requirement 
    arc.years.early <- colnames(arc.sub)[which(earlies <= quantile(earlies, probs = badyear.thres))]
    arc.years.late <- colnames(arc.sub)[which(lates <= quantile(lates, probs = badyear.thres))]
  
  # write-in the ecdf rankings of arc estimates
    arc.Fn <- ecdf(earlies)
    arc.ranks[site,,phase] <- round(arc.Fn(earlies), digits = 2) 
  
    arc.Fn <- ecdf(lates)
    arc.ranks[site,,"Late"] <- round(arc.Fn(lates), digits = 2)  
  # returns vector of rankings overlap results
    print(sapply(prods,vi.compare))
    vi.mat[site,] <- sapply(prods, vi.compare) 
  
  # converting matrix to array for incorporation into master array
  # write ecdf rankings values to master array
    ranks.mat <- sapply(prods, vi.ecdf)
    dim(ranks.mat) <- c(1, length(years), length(prods)) #individual sites
    dimnames(ranks.mat) <- list(site,years,prods)
    vi.rankings[site,,]  <-  ranks.mat
  
  # add non-vi items to matrix 
  # non.vi.mat
} # end of across sites loop 

  

  # write product-level ecdf rankings to csv
    f.out <- paste0(base.path,prods,years[1],"-",years[length(years)],"rankings.csv")
  
  # how could i approach this using a non-for loop approach?
    for (i in 1:length(prods)){
      write.csv(vi.rankings[,,i], f.out[i])
      }
    # comparison of these against existing tables looks good  


  # Create plot lattice with ecdf rankings by year and gps coords
  # hopefully make this work as an array rather than separate early/late matrices - how to do?
  
  # need to add lat/lon & site names as separate entry since previously only used as rownames 
    arc.early <-cbind(arc.ranks[,,"Early"],site.data[,c("Latitude","Longitude")],rownames(site.data))
      colnames(arc.early)[12] <- "site"
      rownames(arc.early) <- NULL
    arc.late <- cbind(arc.ranks[,,"Late"],site.data[,c("Latitude","Longitude")], rownames(site.data))
      colnames(arc.late)[12] <- "site"
      rownames(arc.early) <- NULL
  
  # convert to appropriate form for ggplot2/ggmap use  
    gg.arc.early <- melt(arc.early, id.vars = c("site","Latitude", "Longitude"))
    colnames(gg.arc.early)[5] <- "Ranking"
    cm.key <- "3ba6f5c05bc142209d423981fcbacb4a"  # from Cloud Made API


  # problems when points are outside the bounding box grabbed from OSM - ensure they are inside 
    theme_set(theme_bw(16))

l <- (min(gg.arc.early$Longitude)-0.5)
r <- (max(gg.arc.early$Longitude)+0.5)
t <- (max(gg.arc.early$Latitude)+0.5)
b <- (min(gg.arc.early$Latitude)-0.5)



    Eth <- get_cloudmademap(bbox = c(left = l, bottom = b, right = r, top = t), api_key = cm.key)
    #Eth <- get_map(location = "Ethiopia", zoom = 6, source = c("osm"), api_key = cm.key)
    #ggmap(Eth) + geom_point(data = gg.arc.early, aes(x = Longitude, y = Latitude, color = Ranking)) + facet_wrap(~ variable)
    ggmap(Eth) + stat_density2d(aes(x=Longitude, y = Latitude, fill = Ranking, alpha = ..level..), size = 2, bins = 4, data = gg.arc.early, geom = "polygon") + facet_wrap(~ variable)




  # make sure we have the same order of sites across objects
    if (identical(rownames(site.data), rownames(vi.mat))){
      
      # combine latitude and longitude coordinates for each site  
      # prepare for ggplot2 plotting
      vi.mat <- data.frame(site.data[,c("Latitude", "Longitude")],vi.mat)
      vi.mat <- data.frame(site = rownames(vi.mat), vi.mat)
      rownames(vi.mat) <- NULL
      vi.2 <- melt(vi.mat, id.vars = c("site","Latitude", "Longitude"))
    
      } #end if identical clause
  

  # save table to file
    f.out <- "vi_comparison_table.csv"
    write.csv(vi.mat, file = paste0(base.path, f.out))

  # save sorted version
    vi.sort <- vi.mat[with(vi.mat, order(site)),]  
    f.out <- "vi_comparison_table_sorted.csv"
    write.csv(vi.sort, file = paste0(base.path, f.out)) 
  




# === GRAPHICAL OUTPUT =========================================================
#
#     Here we use a combinations of existing map packages with ggplot2 
#
#===============================================================================

#take the coordinate range to generate a bounding box which is the base layer of the map outputs
lat.r <- c(min(site.data[,"Latitude"]), max(site.data[,"Latitude"])) 
lon.r <- c(min(site.data[,"Longitude"]), max(site.data[,"Longitude"])) 


map <- get_stamenmap(bbox = c(left = lon.r[1], bottom = (lat.r[1]+3), top = (lat.r[2]-1), right = lon.r[2]), zoom = 12, maptype = "toner")
ggmap(map) + geom_point(aes(x= Longitude, y = Latitude, size = RVar), data = vi.mat, colour = 'red') 

eth <- get_map(location = "Ethiopia", zoom = 12, color = "bw", source = "osm")
ethMap <- ggmap(eth, base_layer = ggplot(aes(x = Longitude, y = Latitude), data = vi.2))
ethMap + 
  
  stat_density2d(aes(x = Longitude, y = Latitude, fill = value),bins = 5, geom = "polygon", data = vi.2) + scale_fill_gradient(low = "blue", high = "red") + facet_wrap(~ variable)











######## Variable Grid Size Rank Correlation ######### 

#Enables user to specify a preferred grid size to determine how closely it
#matches to alternative index (e.g., ARC2, RFE, etc.)


#########################################################################
##### User specifies: How many kilometers wide do you want the aggregate gridding box?  

    diameter.range <- 5  #(km) Should be <= 10 km ? To use only the given site pixel, use a value of 0

##########################################################################


    range.number <- 10    # number of values to be computed across the diameter
                          # must be >= 2 

  ## Product pixel sizes (square)
    MODIS.pixel.size  <-  0.00221704 
    SPOT.pixel.size   <-  0.008928572 


long.eq <- 111.32 # number of kilometers between longitude degrees at equator, if we decide to make adjustments to account for shrinking pixel sizes 








# loop through all the sites in sitetable.csv, using GPS coordinates and parameters above to map out an adjacency space 
for (site in rownames(site.data)){
  lat <- site.data[site,"Latitude"]
  lon <- site.data[site,"Longitude"]
#  coord.list <- make.coords(lat,lon,diameter.range,range.number)  #complete list of unique lat/lon pairs given parameter values above - we may not need this since DL simply requires edge values 
#  for (pair in coord.list){  #May not be necessary if only taking edge values
#    # Execute function to get product 
#  }#end of pair for loop
  coords <- coord.range(lat,lon,diameter.range)
  
  
  url.1         <-  "http://iridl.ldeo.columbia.edu/expert/"
  url.2         <-  "SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X"  #standard first part of   
                    #URL query (presumes EVI - have to make more generic)
  url.lonlat    <-  paste0(coords[2,1],"/",coords[2,2],"/RANGEEDGES/Y/",coords[1,1],"/",coords[1,2],"                       
                           /RANGEEDGES[X/Y]average")  #fills in the Lat/Lon specifics
  postxyaddress <-  "T+exch+table-+text+text+-table++.csv"
  address       <-  paste(url.front,url.lonlat,postxyaddress, sep = "/")  
  print(paste("Downloading for", site))
  d.frame       <-  read.csv(address, header = TRUE)

    
   
}# end of site for loop


# SPATIAL CORRELATION ########
# Explanation: Mask out values below a specified threshold
# Theoretically could return a distance value which would inform the earlier function
###############################

# Starting URL http://iridl.ldeo.columbia.edu/expert/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X/37.85/39.85/RANGEEDGES/Y/13.107/15.107/RANGEEDGES/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X/38.85/VALUE/Y/14.107/VALUE/X/removeGRID/Y/removeGRID/exch%5BT%5Dcorrelate/0.9/maskle/figviewer.html?map.url=X+Y+fig-+colors+-fig&my.help=more+options


url.1 <- "http://iridl.ldeo.columbia.edu/expert/"
url.2 <- "SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI"
url.3 <- paste0("X/", lon1,"/",lon2,"/RANGEEDGES/Y/",lat1,lat2,"/RANGEEDGES/")




if (updaterainfall){
  # get last date available in net.cdf file:
  # OLD: download.file("http://iridl.ldeo.columbia.edu/expert/%28/beluga/data/deo/iri_html/iridldata/ethiopiaARCprcp.cdf%29readCDF/.est_prcp/T/last/VALUE//pointwidth/1/def/T/%28%25Y%20%25m%20%25d%29strftimes/data.ch","tempdatefile.csv")
  download.file("http://iridl.ldeo.columbia.edu/expert/home/.deo/.harita/.arc/.precip/T/last/VALUE/pointwidth/1/def/T/%28%25Y%20%25m%20%25d%29strftimes/data.ch","tempdatefile.csv")
  download.file("http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/est_prcp/T/last/VALUE/pointwidth/1/def/T/%28%25Y%20%25m%20%25d%29strftimes/data.ch","tempdatefile.csv")
  lastdate<-read.csv("tempdatefile.csv") # I have to create a stupid tempfile
  # xxx delete tempfile or figure out how to directly load the variable in a future version
  
  print(c("Latest date available in rainfall dataset",lastdate,". All future dates are calculated using average historical rainfall for that day of the year."))
} # end update data

