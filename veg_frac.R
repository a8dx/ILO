
veg.frac.weight <- function(Lat, Lon, Month, BoundingSize){
  environment()  
  # This function uses vegetation fraction values from Marshall Rogers-Martinez'
  # analysis in ENVI.  This is a preliminary view.  
  
  # produces a time-series of monthly, X-Y averaged 1-mo lag of EVI
  # for a specified 3-letter month 
  #
  # Borrows heavily from evi.corr.regrid() 
  # 
  # Lat : user-specified latitude value
  # Lon : user-specified longitude value
  # BoundingSize : value in degrees to create the bounding box: this is added and 
  #                 and subtracted from Lat, Lon coords to create rectangle
  # Month : 3-letter string specifying which months to focus on in identifying worst years 
  
  # create IRI-DL url in stages given input values 
  # keeping as small segments to improve flexibility for future functioning 
  ad1 <- "http://iridl.ldeo.columbia.edu/expert/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.NDVI/T/3837.0/3837.0/RANGE/"
  ad2 <- paste0("X/",Lon-BoundingSize,"/", Lon+BoundingSize, "/RANGEEDGES/")
  ad3 <- paste0("Y/", Lat-BoundingSize, "/", Lat+BoundingSize, "/RANGEEDGES/")
  ad4 <- "dup/home/.ald2187/.ald2187_lvegwgs1984/.the_geom/rasterize/home/.ald2187/.ald2187_lvegwgs1984/.class_id/mul/[gid]/sum/"
  ad5 <- ad2
  ad6 <- ad3
  ad7 <- "T/removeGRID/"
  ad8 <- "SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/"
  ad9 <- ad2
  ad10 <- ad3
  ad11 <- "T/0.9/monthlyAverage/"
  ad12 <- paste0("T/1/1/1/shiftdatashort/")  # stock lag value 
  ad13 <- "mul/"
  ad14 <- paste0("T/(",Month, ")RANGE/")
  ad15 <- "[X/Y]average/"
  ad16 <-  "T+exch+table-+text+text+-table++.csv"
  

  # need to ensure get function calls variables from proper environment
  myfunc <- function(x){
    return(get(x,envir = as.environment(-1)))
  }    
  
  # pieces together segments to generate a unified, usable URL 
  vec <- paste0("ad",1:16) 
  url.name <- capture.output(cat(unlist(lapply(vec,myfunc)), sep = "", collapse = ""))
  
  # specify output .csv filename and location
  fout <- paste0(out.path,"evicorr_temp.csv")
  download.file(url.name, fout, cacheOK = FALSE)
  return(read.csv(fout, header = T))
    
} # end veg.frac.weight
