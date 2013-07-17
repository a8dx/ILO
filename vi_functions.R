# Created:    June 26, 2013
# Purpose:    Auxiliary Functions Used in VI Comparisons 
# Author:     Anthony D'Agostino (ald2187@columbia.edu)
# Last Edit:  July 1, 2013 



# EDIT HISTORY
#
# July 15 (ALD) - modified URL in evi.corr.regrid to better reflect the calculation
#                 of masked EVI values.  

library(gtools)

vi.common <- function(product = c("EVI")){
  # Product has to be a single product, default is EVI unless otherwise specified
  # this code [will be] common to several functions in the vi_functions file and 
  # read in VI data from existing .csv's  
  #
  # product is "EVI", "NDVI", etc.
  # phase and satellite must already be defined: satellite is "Mod" or "Spo", phase is "early" or "late"  
  fn <- paste0(base.path, scen,site,"/met.data/",phase,product,satellite,"R.csv")
  
  # need to create this if loop for sites which don't have complete satellite data 
  if (file.exists(fn)){
    vi <- read.csv(fn, header = T, row.names = 1)  # read appropriate VI file 
    # tryCatch(stop("you threw an error"), error = function(e) print(e$message))
    rownames(vi) <- substr(rownames(vi), 5, 9)  # removing month from rownames 
    return(vi[years,])  # take only the relevant years 
  }
  else {return(NA)}

}




vi.compare <- function(product){ 
  
  # This function reads in existing VI data from a scenario folder, identifies
  # which years are the worst (according to a pre-defined "badyears.thres"), and
  # then returns the percentage of years in which these worst years correspond
  # to ARC data
  
  # product is "EVI", "NDVI", etc.
  # phase and satellite must already be defined: satellite is "Mod" or "Spo", phase is "early" or "late"  
  
    if(any(!is.na(vi.common(product)))){ #returns NA if file does not exist 
      
      vi <- vi.common(product)
    
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
  # `phase' and `satellite' must already be defined: `satellite' is "Mod" or "Spo", `phase' is "early" or "late"  
  
  if(any(!is.na(vi.common(product)))){ #returns NA if file does not exist 
      vi <- vi.common(product)
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
  
  
  
  make.coords2 <- function(Lat,Lon,Range,Points){ 
    #make.coords creates a matrix of all evenly-spaced lat/lon pairs, to create a grid of VI data values 
    coords.mat <- matrix(data = NA, ncol = 2)
    lat.range <- seq(Lat - Range / (2*long.eq), Lat + Range / (2*long.eq), by = (2*Range/long.eq)/(Points-1))    #extent of lat range
    lon.range <- seq(Lon - Range / (2* long.eq), Lon + Range / (2*long.eq), by = (2*Range/long.eq)/(Points-1))    #extent of lon range
    coords.mat <- expand.grid(lat.range,lon.range)
    colnames(coords.mat) <- c("Latitude", "Longitude")
    return(coords.mat)
  }
  


make.coords <- function(df,Buffer){
  # Reads in the max/min lat/lon of a df with Latitude and Longitude columns 
  # and generates a list of bounding coordinates to be read into a ggmap setting, 
  # given a user-specified buffer, e.g., (0.1 degrees)
  
  return(list(l = min(df$Longitude)-Buffer, b = min(df$Latitude)-Buffer, r = max(df$Longitude)+Buffer, t = max(df$Latitude)+Buffer))

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



evi.corr.regrid <- function(Lat,Lon,Size,CorrThreshold,Month,RegridSize=0.05,Lags = FALSE){
  environment()  
  # produces a time-series of monthly, X-Y averaged 1-mo lag of EVI
    # for a specified 3-letter month 
    # Lat : user-specified latitude value
    # Lon : user-specified longitude value
    # Size : value in degrees to create the bounding box
    # RegridSize : in degrees, what scale is he underlying MODIS data regridded to
    # CorrThreshold : threshold value for which pixel-level correlations between EVI and ARC under this value are masked out as NA 
    # Month : 3-letter string specifying which months to focus on in identifying worst years 
    # Lags : enables capture of XY average correlation values across pre-specified lags
      
    # create IRI-DL url in stages given input values 
    # keeping as small segments to improve flexibility for future functioning 
  lagval <- paste0("T/",lag.start,"/1/",lag.end,"/shiftdatashort")
      ad1 <- "http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/.est_prcp"
      ad2 <- paste0("/X/",Lon,"/VALUE/Y/",Lat,"/VALUE/")
      ad3 <- "X/removeGRID/Y/removeGRID/T/1.0/monthlyAverage"
      ad4 <- "/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/"
      ad5 <- paste0("X/",Lon-Size,"/", Lon+Size, "/RANGEEDGES/")
      ad6 <- paste0("Y/", Lat-Size, "/", Lat+Size, "/RANGEEDGES/")
      ad7 <- paste0("X/", RegridSize, "/", "0.9/boxAverage/Y/", RegridSize, "/", "0.9/boxAverage/")
      ad8 <- "T/0.9/monthlyAverage/"
      ad9 <- "dup/"   # have to duplicate before taking the lag 
      ad10 <- paste0("T/1/1/1/shiftdatashort")  # stock lag value 
      ad11 <- "/3/-1/roll"
      ad12 <- "/[T]/correlate/"
      ad13 <- paste0(CorrThreshold, "/maskle/") 
      ad14 <- "dataflag/0/maskle/mul/"
      ad15 <- "[X/Y]average/"
      ad16 <- paste0("T/(",Month, ")RANGE/")
      ad17 <- "T+exch+table-+text+text+-table++.csv"
      
      # need to ensure get function calls variables from proper environment
  myfunc <- function(x){
    return(get(x,envir = as.environment(-1)))
  }    
  
  # modifying url to reflect capture only of correlation values 
  if (Lags == TRUE){
    ad7 <- ""
    ad9 <- lagval
    ad10 <- "" 
    ad11 <- ""
    ad12 <- ad12
    ad13 <- ""
    ad14 <- ""
    ad15 <- ad15
    ad16 <- ""
    ad17 <- "T_lag+exch+table-+text+text+-table++.csv"
}
  
  
  
  # lag URL will look something like this: "http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/.est_prcp/X/38.726/VALUE/Y/7.845/VALUE/X/removeGRID/Y/removeGRID/T/1.0/monthlyAverage/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X/38.526/38.926/RANGEEDGES/Y/7.645/8.045/RANGEEDGES/T/1.0/monthlyAverage/T/0/1/0/shiftdatashort/[T]/correlate/[X/Y]average/  T_lag"
  
  
  # pieces together segments to generate a unified, usable URL 
    vec <- paste0("ad",1:17) 
    url.name <- capture.output(cat(unlist(lapply(vec,myfunc)), sep = "", collapse = ""))

  # specify output .csv filename and location
    fout <- paste0(out.path,"evicorr_temp.csv")
      download.file(url.name, fout)
      return(read.csv(fout, header = T))
  } # end evi.corr.regrid function 



# more work needs to be done on this.  
evi.regrid <- function(GridScale){
  # gridscale is a value in degrees 
 
   out <- eval(evi.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Month = month, CorrThreshold = 0, Size = GridScale, RegridSize = MODIS.pixel.size))
  return(out)
} # end evi.regrid 






df.melt <- function(dframe.pre){
  # dframe = a standard dataframe which has not yet been melted for gg use 
  # because of scoping rules (i believe), have to melt the dataframe
  # outside of arc.vi.vis to avoid an environment problem 
  #
  return(melt(dframe.pre, id.vars = c("site","Latitude", "Longitude")))
}



arc.vi.vis <- function(x.df){
  #     
  # dataframe must have Latitude and Longitude columns, appropriately named
  # imports buffer size (b.s) variable from main script 
    #
  cm.key <- "3ba6f5c05bc142209d423981fcbacb4a"  # from Cloud Made API
  
  
  
  l <- min(x.df$Longitude) - b.s
  b <- min(x.df$Latitude) - b.s
  r <- max(x.df$Longitude) + b.s
  t <- max(x.df$Latitude) + b.s

  
  if(override){
    l <- override.box[1]; b <- override.box[2]; r <- override.box[3]; t <- override.box[4]
    
  }
  
  # convert to appropriate form for ggplot2/ggmap use  
    colnames(x.df)[5] <- "Ranking"
  
  Latitude <- x.df$Latitude
  Longitude <- x.df$Longitude
  
  theme_set(theme_bw(16))
  outmap <- get_cloudmademap(bbox = c(left = l, bottom = b, right = r, top = t), api_key = cm.key)  
  
  # modify this portion if want to change gradient colors 
    ggmap(outmap) + geom_point(aes_string(x = Longitude, y = Latitude, color = Ranking), data = x.df) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable)
          
} # end arc.vi.vis 

find.worst <- function(obj,colname){
  # Returns the years of an obj with column "Year" falling below passed in 'badyear.thres' value
 # # Requires obj to also have "correlation" column 
  # Could generalize and include an input parameter corresponding to column to be ranked
  Fn <- ecdf(obj[,colname]) 
  ranks <- Fn(obj[,colname])
  return(obj$Year[which(ranks <= badyear.thres)])
}


bench.corr.compare <- function(corr.value){
    # corr.value corresponds to an earlier run which denotes a unique .csv file
    # with output from the ARC-VI spatial correlation section of the ILO report
    # script 
  
  df.both <- read.csv(paste0(out.path,"VIspatialcorrelation_",corr.value,"_.csv"), header = T)
  df.early <- df.both[which(df.both[,"Window"] == "Early"),]      
  df.late <- df.both[which(df.both[,"Window"] == "Late"),]
  
  if (identical(df.early[,"site"],benchmark.early[,"site"])){
    print("Sites match, Early comparison success!")
    out.early <- data.frame(df.early, fin = (df.early[,"Agreement"]-benchmark.early[,"value"]))
  }else{print("Mismatch in site list - early window comparison failed")}
 
      
    if (identical(df.late[,"site"],benchmark.late[,"site"])){print("Sites match, Late comparison success!")
        out.late <- data.frame(df.late, fin = (df.late[,"Agreement"]-benchmark.late[,"value"] ))
    }else{print("Mismatch in site list - late window comparison failed")}

  
        out <- rbind(out.early,out.late)
      # use common naming system for different outputs   
        fn <- paste0(out.path,"CorrThreshPerfDiff_",corr.value,"_.")
        
        ti.c <- paste0("Agreement Performance Gain [ ", corr.value, "Correlation Threshold - Benchmark]")
    # save as .png and .csv
        ggmap(outmap) + geom_point(aes(x = Longitude, y = Latitude, color = fin), data = out) + scale_color_gradient(low = "red", high = "green") + labs(title = ti.c)
        ggsave(file = paste0(fn,"png"))
        write.csv(out,paste0(fn,"csv"))
        return(out)
  }  # end of bench.corr.compare function 

# read.csv
# subtract from benchmark, both early and late
# save output in ggmap format as well as .csv

