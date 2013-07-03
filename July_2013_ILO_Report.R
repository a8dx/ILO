# ******************************************************************************
#
#       Vegetation Index / Rainfall Estimate Validation Procedures 
#       Created by: Anthony Louis D'Agostino (ald2187@columbia.edu)
#       Date Created: June 18, 2013  
#
# ******************************************************************************




# =========================  Edit History  ===================================== 
#
# 7/2 (ALD) - improved plotting using ggmap, streamlined documentation 
# 7/1 (ALD) - created evi.corr.regrid function in vi_functions.R, auxiliary to
# this file

# ******************************************************************************



# ===========================  Purpose  ======================================== 
# 
# This script performs 3 discrete tasks 1. Calls existing rainfall estimate and
# VI data from a sniidharita scenario to  emulate results produced in earlier
# ILO reports. 2. Creates a rankings comparion against ARC2 for a
# variable-sized, upscaled VI pixel (i.e., user can specify whether or not to
# spatially average over a 2km x 2km grid to determine correlations and payout
# agreement with rainfall). 3. Creates a rankings comparison against ARC2 by
# masking out pixels with VI spatial correlation coefficients below some
# user-specified level.  The remaining pixels are spatially averaged to produce
# a single VI time-series per site.
# 
# ******************************************************************************




# ============================  Requirements ===================================
# 
# * sitetable.csv - copied from Dry2012Satellite's common.data folder under the
# sniidharita repository, current as of June 18, 2013.  Will identify any
# discrepancies between this script's output and the earlier ILO reporting data 
#
# * ARC2[/RFE] rainfall files from sniidharita script
# 
# ******************************************************************************



# ============================= Notes ========================================== 
#
# Datasets may be applicable only for some regions, check beforehand
# Need to determine which rainfall products we are comparing against before pushing forward with design
#
# * Currently only focusing on ARC2 data, but have also included RFE rankings 
#
# VARIABLES OF INTEREST
# * vi.mat = matrix of agreement percentage with ARC data by VI product
# * arc.ranks = for validating earlier ILO results: gives ECDF ranks across 
#     early and late phases 
# ******************************************************************************


# ============================== Tasks ========================================= 
#
# * Variable lag length for VI data against rainfall (creating a vector with multiple lags enables us to do immediate comparisons - this will be in terms of the 16-day composites)  Also consider taking the average of multiple lags (e.g., 1st 16-day composite lag behind, and 2nd)
# * Masking out pixels with spatial correlation values below threshold
# * Changing aggregate pixel size (may be based on some criteria) - currently all pixels within grid are evenly weighted
# * Identify which pixels fall below a threshold correlation level with ARC2 and then determine whether those pixels have a better correspondence with the worst years                
#   * Determine how best to save the pixel-level data in a meaningful way 
#Apply some weighting criteria to grid based on distance to village pixel 
# Determine appropriate graphing facilities for exporting data  
# User specifies which products interested in, use "if" statements to match product with URL assembly

#???? Uses some of the core existing functions written for sniidharita 
#  
# ******************************************************************************



# ---- Initialization --------------------

    rm(list = ls())
    library(zoo)
    library(ggmap)  # for visualizations
    library(ggplot2)
    library(reshape2) 
 
#=============================================================================#
  # Specify folder containing this script and related files 
    setwd("~/Documents/Index_Insurance/ILO/ILO_github/")

  # load auxiliary functions like "vi.compare" 
    source("vi_functions.R")

  # Names and creates (in first run) folder for outputted rankings files 
    out.path <- "Output"
    if (!file.exists(out.path)){dir.create(out.path)}
#=============================================================================#


#======================= Validation Comparison ===============================#
#
# Comparing worst years of ARC2 against select VI products.  Results should 
# match output in earlier ILO reports.  
#
#=============================================================================#



# ---- Read in Precipitation Data -----

  # range of years covered in analysis, specified in Phase One report 
    years <- as.character(seq(2003,2011))  
  
  # choose the scenario to work with (draws upon existing file/folder structure)
    base.path <- "~/Documents/Index_Insurance/sniidharita_fist/WIIET_Data/Dry2013Initial/"
  # read in sitetable and common dekad files 
    site.data <- read.csv(paste0(base.path,"common.data/sitetable.csv"), header = TRUE, row.names = 1)
    dekad.cal <- read.csv(paste0(base.path,"common.data/CalendarYearDayDecad.csv"), header = F, skip = 2, row.names = 1)
    dekadmonth <- read.csv(paste0(base.path, "common.data/dekadmonth.csv"))
  
  # folder name under which site-level data is located
    scen <- "scenarios/"
  # folder structure for rbf files 
    rbf.fn <- "/EthRFEadj/rainbyphase.csv"

 



  # establish parameters for comparison table 
    prods <- c("EVI", "NDVI", "NDWI")  
    satellite <- "Mod"  
    phase <- "Early" 
    badyear.thres <- 1/6

  # create array for early and late arc ecdf rankings 
    arc.ranks <- array(data = NA, dim = c(nrow = dim(site.data)[1], length(years), 2))
    dimnames(arc.ranks) <- list(rownames(site.data), years, c("Early", "Late"))
  

  # create a matrix where comparison data for a single product at a time will be kept 
    vi.mat <- matrix(ncol = length(prods), nrow = length(rownames(site.data)))
    rownames(vi.mat) <- rownames(site.data)
    colnames(vi.mat) <- prods

  # create an array with rankings for each product, each site, across all years
    vi.rankings <- array(data = NA, dim = c(length(rownames(site.data)), length(years), length(prods)))
    dimnames(vi.rankings) <- list(rownames(site.data), years, prods)

  # begin comparison process on a site-by-site basis 
    for (site in rownames(site.data)){
      
      # for troubleshooting purposes 
       print(paste0("Now calculating comparisons for ",site))
  
      # read in rainfall by phase data  
        rfe <- read.csv(paste0(base.path,scen,site,rbf.fn), header = T, row.names = 1)
        rfe.sub <- rfe[years,]  #only years we're interested in 
        rfe.p1 <- rfe.sub[years,"X1"] # first phase 
        rfe.p2 <- rfe.sub[years,"X3"] # second phase 
       
      # generate rankings 
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
        arc <- read.csv(paste0(base.path,scen,site,"/met.data/precip.daily.csv"), header = T, row.names = 1)
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
       
      # returns vector of rankings agreement (overlap) results
      # uses "vi.compare" from vi_functions script  
        print(sapply(prods,vi.compare))
        vi.mat[site,] <- sapply(prods, vi.compare) 
      
      # converting matrix to array for incorporation into master array
      # write ecdf rankings values to master array
        ranks.mat <- sapply(prods, vi.ecdf)
        dim(ranks.mat) <- c(1, length(years), length(prods)) #individual sites
        dimnames(ranks.mat) <- list(site,years,prods)
        vi.rankings[site,,]  <-  ranks.mat
    
  } # end of across sites loop 

  # write product-level ecdf rankings to csv
  # how could i approach this using a non-for loop approach?
    f.out <- paste0(out.path,prods,years[1],"-",years[length(years)],"rankings.csv")
    for (i in 1:length(prods)){
      write.csv(vi.rankings[,,i], f.out[i])
      }
    # comparison of these against existing tables looks good  
  
  # create all-site rankings table with lat/lon coordinates 
  # need to add lat/lon & site names as separate entry since previously only used as rownames 
    arc.early <-cbind(arc.ranks[,,"Early"],site.data[,c("Latitude","Longitude")],rownames(site.data))
      colnames(arc.early)[12] <- "site"
      rownames(arc.early) <- NULL
    arc.late <- cbind(arc.ranks[,,"Late"],site.data[,c("Latitude","Longitude")], rownames(site.data))
      colnames(arc.late)[12] <- "site"
      rownames(arc.early) <- NULL
  
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
    write.csv(vi.mat, file = paste0(out.path, f.out))
    
  # save sorted version
    vi.sort <- vi.mat[with(vi.mat, order(site)),]  
    f.out <- "vi_comparison_table_sorted.csv"
    write.csv(vi.sort, file = paste0(out.path, f.out)) 
    



# ---- Visualization --------------------

# set bounding box coordinates for map background. 
# May consider modifying to exclude Michael Debir for crisper presentation. 
# Problems arise when df points are outside the bounding box grabbed from OSM, 
# therefore ensure they are inside.


  # API key from Cloud Made Maps
    cm.key <- "3ba6f5c05bc142209d423981fcbacb4a"  


  b.s <- 0.2  # buffer size in degrees, for how large to create the bounding box around the selected bounding box edge pixels  

override.box <- c("L","B","R","T")    # active only when override = TRUE, vector sequence is left, bottom, right, and top pixels (L,R: longitude, T,B: latitude) - must replace with values in degrees 
    
    override <- FALSE   # if TRUE, overrides default calculation of bounding box 
                        #dimensions in the viz function 

 
    if(override){
      l <- override.box[1]; b <- override.box[2]; r <- override.box[3]; t <- override.box[4]
    }
else
  {    # default bounding box coordinates, including the earlier specified b.s
  l <- min(arc.early$Longitude) - b.s
  b <- min(arc.early$Latitude) - b.s
  r <- max(arc.early$Longitude) + b.s
  t <- max(arc.early$Latitude) + b.s
  }


  # convert input dataframe to appropriate form for ggplot2/ggmap use  
    # arc early rankings 
    a_e <- melt(arc.early, id.vars = c("site","Latitude", "Longitude"))
    colnames(a_e)[5] <- "Ranking"
    # arc late rankings 
    a_l <- melt(arc.late, id.vars = c("site","Latitude", "Longitude"))
    colnames(a_l)[5] <- "Ranking"


  # set standard map base layers   
    theme_set(theme_bw(16))
    outmap <- get_cloudmademap(bbox = c(left = l, bottom = b, right = r, top = t), api_key = cm.key)  

  # generate early window ARC rankings 
    print(ggmap(outmap, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = a_e)) + geom_point(aes(color = Ranking)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable) + labs(title = "ARC Early Window Rankings"))    

  # generate late ARC ranks 
    print(ggmap(outmap, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = a_l)) + geom_point(aes(color = Ranking)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable) + labs(title = "ARC Late Window Rankings"))


#Comparing worst years of ARC2 against select VI products.  Results should 
# match output in earlier ILO reports.


# ====  These functions should work, but encountering problems =====
  # Have opted to pull out the commands outside the function - can improve for later iteration 
#  arc.early <- df.melt(arc.early)  # creates a melted version for gg
#  arc.late <- df.melt(arc.late)  # creates a melted version for gg
#  arc.vi.vis(arc.early)
#  arc.vi.vis(arc.late)



# ======== Comparing VI and ARC ranks ======= 
# 
# Also used for validating against previous ILO reports, this time comparing 
# individual VI products against ARC and with the badyear.thres value defining
# which are the worst years 
#
# ******************************************************************************

  # ARC-VI Agreement on Worst Years
    # this can serve as the benchmark against which other versions can be compared to evaluate any performance gains from regridding, spatial correlation, etc.  
    print(ggmap(outmap, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = vi.2)) + geom_point(aes(color = value)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable) + labs(title = "ARC-VI Worst Year Agreement %"))









# working version of point-based output 
# print(ggmap(Eth, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = gg.arc.early)) + geom_point(aes(color = factor(Ranking), size = 3)) + facet_wrap(~ variable))  

# working version using scale gradient 
#print(ggmap(Eth, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = gg.arc.early)) + geom_point(aes(color = Ranking)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable))  




# VI-ARC agreement rankings 


    
    



#http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/.est_prcp/X/38.85/VALUE/Y/14.1509/VALUE/X/removeGRID/Y/removeGRID/T/1.0/monthlyAverage/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X/33/0.1/48/GRID/Y/3.4/0.1/14.9/GRID/T/1.0/monthlyAverage/T/0/1/3/shiftdatashort%5BT%5Dcorrelate/0.7/maskle/figviewer.html?my.help=more+options&map.T_lag.plotvalue=0.0&map.Y.units=degree_north&map.Y.plotlast=14.95N&map.url=X+Y+fig-+colors+-fig&map.domain=+%7B+/T_lag+0.0+plotvalue+%7D&map.domainparam=+/plotaxislength+432+psdef+/plotborder+72+psdef&map.zoom=Zoom&map.Y.plotfirst=3.35N&map.X.plotfirst=32.95E&map.X.units=degree_east&map.X.modulus=360&map.X.plotlast=48.05E&map.correlation.plotfirst=-0.8637556&map.correlation.units=unitless&map.correlation.plotlast=0.8637556&map.newurl.grid0=X&map.newurl.grid1=Y&map.newurl.land=draw+...&map.newurl.plot=colors&map.plotaxislength=432&map.plotborder=72&map.fnt=Helvetica&map.fntsze=12&map.XOVY=auto&map.color_smoothing=1&map.framelbl=framelabelstart&map.framelabeltext=&map.iftime=25&map.mftime=25&map.fftime=200




# ======================= VI Scaling Exercises ===============================#
#
#   
#
# How many km wide do you want the aggregate gridding box?  [0 < x <= 10 km]
#
    diameter.range <- 5  
#
#=============================================================================#






# ---- Graphical Output -----
#
# Here we use a combinations of existing map packages with ggplot2 
#








#======================= ARC - VI SPATIAL CORRELATION =========================#
#
# Uses a correlation threshold 'r' such that all pixels whose ARC-VI correlation 
# coefficients exceeding 'r' inside the specified bounding box are averaged 
# out, generating a single VI value for each time period.  
#

  # specify the correlation coefficient threshold for masking out pixels 
        corr.value <- 0.7

  # specify years of interest - uncomment 1st to reuse years from beginning of analysis 
      # years <- years   
        years <- as.character(seq(2000,2012))  
  
#=============================================================================#



  # create template df from which early/late windows will be drawn - final agreement % values will be stored here as we move through the loop 
    agree.df <- data.frame(Latitude = site.data[,"Latitude"], Longitude = site.data[,"Longitude"], Agreement = NA)
    rownames(spatial.agree.df) <- rownames(site.data)
    
    agree.early <- agree.df 
    agree.late <- agree.df 

# PSEUDO - much of this needs to be turned into a function since repeated for 
# both early and late windows 

for (site in rownames(site.data)){
  # month of VI values determined in addVEG*.R script : reads in dekad month 
  # 
  
  # have to read in site-specific contract parameters to determine which dekads to pull VI products for 
  contract.pth<-paste0(base.path,scen,site,"/payout.data/contract.R")
  source(contract.pth)
  midearly<-as.integer(t$swfirst+(t$phases[2,1]+t$phases[1,1])/2)
  midlate<-as.integer(t$swfirst+(t$phases[2,3]+t$phases[1,3])/2)
  
  # 4 dekads later was the [adjustable] standard delay 
  # delay already inherent in IRI Data Library code in evi.corr.regrid 
  dekdelay <- 0 
  earlymonth <- as.character(dekadmonth[(midearly+dekdelay)%%36,"Month"])
  latemonth <- as.character(dekadmonth[(midlate+dekdelay)%%36,"Month"])
  
  # generates lagged EVI values for all years available   
  evi.early <- data.frame(evi.corr.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Size = b.s, CorrThreshold = corr.value, Month = earlymonth), Window = "Early") 
  evi.late <- data.frame(evi.corr.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Size = b.s, CorrThreshold = corr.value, Month = latemonth), Window = "Late") 
  
  # now subset to the years of interest and identify which years agree with ARC worst
  # early.worst / late.worst produce worst years, while *.ranks give total ranks vectors
    evi.early <- data.frame(evi.early, Month = substr(evi.early$Time,1,3)[1], Year = as.numeric(substr(evi.early$Time, 5, 8)))
    evi.late  <- data.frame(evi.late,  Month = substr(evi.late$Time,1,3)[1], Year = as.numeric(substr(evi.late$Time, 5, 8)))
  
  # drop the original "Time" column 
    evi.early <- subset(evi.early, select = -Time)
    evi.late  <- subset(evi.late, Year == years, select = -Time)
  
  # PSEUDO - create code to ensure that all years are providing data for the same month

    evi.early <- evi.early[is.element(evi.early$Year,years),]
    evi.late  <- evi.late[is.element(evi.late$Year,years),]
  
    Fn <- ecdf(evi.early$flag) # since under the function, IRIDL generates "flag" 
    early.ranks <- round(Fn(evi.early$flag),digits = 2)
    early.worst <- evi.early$Year[which(early.ranks <= badyear.thres)]
  
    Fn <- ecdf(evi.late$flag) # since under the function, IRIDL generates "flag" 
    late.ranks <- round(Fn(evi.late$flag),digits = 2)
    late.worst <- evi.late$Year[which(late.ranks <= badyear.thres)]
    
  print(paste0("In ", site, " the worst ", round(badyear.thres, digits = 2), "early window years for the ", years[1], years[length(years)], " period are ", early.worst))
  
  print(paste0("In ", site, " the worst ", round(badyear.thres, digits = 2), "late window years for the ", years[1], years[length(years)], " period are ", late.worst)) 
  
  # PSEUDO - now run an agree function to find out what percent of years are matching     
  
    # looking at arc.years.early and arc.years.late 
  
  arc.early <- as.numeric(arc.years.early)
  arc.late <- as.numeric(arc.years.late)
  
    agree.early[site,"Agreement"] <- length(intersect(arc.early,early.worst)) / max(length(arc.years.early), length(early.worst))
  
  agree.late[site,"Agreement"] <- length(intersect(arc.late,late.worst)) / max(length(arc.years.late), length(late.worst))
  } # end of spatial correlation for loop across sites 

#PSEUDO - want to compare these results against the baseline to identify any performance improvements 

#PSEUDO - insert visualization function here, when it's working 
arc.evi <- rbind(cbind(agree.early, Window = "Early"), cbind(agree.late, Window = "Late"))

ti <- paste0("ARC-EVI Agreement for Pixels with R > ", corr.value)

print(ggmap(outmap, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = arc.evi)) + geom_point(aes(color = Agreement)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ Window) + labs(title = ti))    


























######## Variable Grid Size Rank Correlation ######### 

#Enables user to specify a preferred grid size to determine how closely it
#matches to alternative index (e.g., ARC2, RFE, etc.)







    range.number <- 10    # number of values to be computed across the diameter
                          # must be >= 2 

  ## Product pixel sizes (square)
    MODIS.pixel.size  <-  0.00221704 
    SPOT.pixel.size   <-  0.008928572 


long.eq <- 111.32 # number of kilometers between longitude degrees at equator, if we decide to make adjustments to account for shrinking pixel sizes 











