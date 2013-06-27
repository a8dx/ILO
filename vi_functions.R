# Created:    June 26, 2013
# Purpose:    Auxiliary Functions Used in VI Comparisons 
# Author:     Anthony D'Agostino (ald2187@columbia.edu)
# Last Edit:  June 26, 2013 




vi.common <- function(product){
  # product is "EVI", "NDVI", etc.
  # phase and satellite must already be defined: satellite is "Mod" or "Spo", phase is "early" or "late"  
  fn <- paste0(scen,site,"/met.data/",phase,product,satellite,"R.csv")
  
  # need to create this if loop for sites which don't have complete satellite data 
  if (file.exists(fn)){
    vi <- read.csv(fn, header = T, row.names = 1)  # read appropriate VI file 
    # tryCatch(stop("you threw an error"), error = function(e) print(e$message))
    rownames(vi) <- substr(rownames(vi), 5, 9)  # removing month from rownames 
    return(vi[years,])  # take only the relevant years 
  }

}




vi.compare <- function(product){ 
  
  # This function reads in existing VI data from a scenario folder, identifies
  # which years are the worst (according to a pre-defined "badyears.thres"), and
  # then returns the percentage of years in which these worst years correspond
  # to ARC data
  
  # product is "EVI", "NDVI", etc.
  # phase and satellite must already be defined: satellite is "Mod" or "Spo", phase is "early" or "late"  
  
  fn <- paste0(scen,site,"/met.data/",phase,product,satellite,"R.csv")
  
  # need to create an if loop for sites which don't have complete satellite data 
  if (file.exists(fn)){
    vi <- read.csv(fn, header = T, row.names = 1)  # read appropriate VI file 
    # tryCatch(stop("you threw an error"), error = function(e) print(e$message))
    rownames(vi) <- substr(rownames(vi), 5, 9)  # removing month from rownames 
    vi <- vi[years,] # take only the relevant years 
    
    # Identify which years are bad according to pre-defined threshold 
      vi.bad.years <- years[which(vi <= quantile(vi, probs = badyear.thres))]
      if (phase == "Early"){
        # calculate overlap 
        ol <- intersect(arc.years.early,vi.bad.years)
        pct <- length(ol) / max(length(arc.years.early), length(vi.bad.years))
      } 
      if (phase == "Late"){
        ol <- intersect(arc.years.late,vi.bad.years)
        pct <- length(ol) / max(length(arc.years.late), length(vi.bad.years))   
      }
  }
    
    # if VI file does not exist 
  else {pct = NA}
  
return(pct)
}  #end of vi.compare function 


vi.ecdf <- function(product){
  # Reads in VI data for single phase, single sat products 
  # Calculates ranks using 'ecdf' 
  # Output table is generated in main script 

  # product is "EVI", "NDVI", etc.
  # phase and satellite must already be defined: satellite is "Mod" or "Spo", phase is "early" or "late"  
  
  # read in existing vi data file 
    fn <- paste0(scen,site,"/met.data/",phase,product,satellite,"R.csv")
  
  # need to create an if loop for sites which don't have complete satellite data 
    if (file.exists(fn)){
      vi <- read.csv(fn, header = T, row.names = 1)  # read appropriate VI file 
      # tryCatch(stop("you threw an error"), error = function(e) print(e$message))
      rownames(vi) <- substr(rownames(vi), 5, 9)  # removing month from rownames 
      vi <- vi[years,] # take only the relevant years 
    
      Fn1 <- ecdf(vi)
      res <- round(Fn1(vi), digits = 2) 
    }

  
  # if VI file does not exist 
    else {res <- rep(NA,length(years))}
    
  return(res)  
}
  
  
  
  long.cor <- function(longitude,latitude){  
    
    # Calculates distance between longitude degrees 
    # Need to determine whether these corrections get included in final version 
    
    return(cos(latitude/(180/pi))*long.eq)
  }
  
  
  
  make.coords <- function(Lat,Lon,Range,Points){ 
    #make.coords creates a matrix of all evenly-spaced lat/lon pairs, to create a grid of VI data values 
    coords.mat <- matrix(data = NA, ncol = 2)
    lat.range <- seq(Lat - Range / (2*long.eq), Lat + Range / (2*long.eq), by = (2*Range/long.eq)/(Points-1))    #extent of lat range
    lon.range <- seq(Lon - Range / (2* long.eq), Lon + Range / (2*long.eq), by = (2*Range/long.eq)/(Points-1))    #extent of lon range
    coords.mat <- expand.grid(lat.range,lon.range)
    colnames(coords.mat) <- c("Latitude", "Longitude")
    return(coords.mat)
  }
  
  coord.range <- function(Lat,Lon,Range){ 
    #Provides the edge values of a lat-lon grid given user-specified distance value from site pixel
    coords.out <- matrix(data = NA, nrow = 2, ncol = 2)
    lat.range <- c(Lat - Range / (2*long.eq), Lat + Range / (2*long.eq))
    lon.range <- c(Lon - Range / (2*long.eq), Lon + Range / (2*long.eq))
    coords.out <- rbind(lat.range,lon.range)
    rownames(coords.out) <- c("Latitude","Longitude")
    return(coords.out)
  }