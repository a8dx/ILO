# ******************************************************************************
#
#       Vegetation Index / Rainfall Estimate Validation Procedures 
#       Created by: Anthony Louis D'Agostino (ald2187@columbia.edu)
#       Date Created: June 18, 2013  
#
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
# Some plots will be sent straight to PDF in the Output folder 
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



# =========================  Edit History  ===================================== 
# 7/9 (ALD) - added both woreda and site-level map plots so that in write-up
#             we can refer to specific areas when making comments about results        
# 7/3 (ALD) - further documentation work, some pseudo-coding for necessary work
#             problems accessing ggplot2 functions inside function, requires 
#             some repetition of procedures - need to discuss with Helen 
# 7/2 (ALD) - improved plotting using ggmap, streamlined documentation 
# 7/1 (ALD) - created evi.corr.regrid function in vi_functions.R, auxiliary to
# this file

# ******************************************************************************




# ============================= Notes ========================================== 
#
# Datasets may be applicable only for some regions, given how IRI-DL data is 
# divided into tiles - check beforehand.  
# Need to determine which rainfall products we are comparing against before 
# pushing forward with design -> should move to creating generic functions 
#
# * Currently only focusing on ARC2 data, but have also included RFE rankings 
# * May encounter download problems with IRI-DL, solution is to try again.  
#
# VARIABLES OF INTEREST
# * vi.mat = matrix of agreement percentage with ARC data by VI product
# * arc.ranks = for validating earlier ILO results: gives ECDF ranks across 
#     early and late phases 
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
    out.path <- "Output/"
    if (!file.exists(out.path)){dir.create(out.path)}
#=============================================================================#


#===================== ARC Worst Years Assessment =============================
#
# Comparing worst years of ARC2 against select VI products.  Results should 
# match output in earlier ILO reports.  
#
#=============================================================================#



# ---- Read in Precipitation Data -----

  # range of years covered in analysis, specified in Phase One report 
    years <- as.character(seq(2000,2011))  
  
  # choose the scenario to work with (draws upon existing file/folder structure)
    base.path <- "~/Documents/Index_Insurance/sniidharita_fist/WIIET_Data/Dry2012satellite/"
  # read in sitetable and common dekad files 
    site.data <- read.csv(paste0(base.path,"common.data/sitetable.csv"), header = TRUE, row.names = 1)
    dekad.cal <- read.csv(paste0(base.path,"common.data/CalendarYearDayDecad.csv"), header = F, skip = 2, row.names = 1)
    dekadmonth <- read.csv(paste0(base.path, "common.data/dekadmonth.csv"))
  
  # folder name under which site-level data is located
    scen <- "scenarios/"
  # folder structure for rbf files 
    rbf.fn <- "/EthRFEadj/rainbyphase.csv"

 
# PSEUDO - 'phase' is passed to vi.common (in vi_functions), yet we are really 
# interested in assessing all phases.  Here I take the short-cut of focusing 
# on a single phase since this is about reproducing an earlier report, not
# developing a new core function.  Should be fixed in the future.  

  # establish parameters for comparison table 
    prods <- c("EVI", "NDVI", "NDWI")  
    satellite <- "Mod"  
    phases <- c("Early","Late")
    phase <- "Early"
  # compare the worst x % of years for agreement between products 
    badyear.thres <- 1/6    
  # create array for early and late arc ecdf rankings 
    arc.ranks <- array(data = NA, dim = c(nrow = dim(site.data)[1], length(years), length(phases)))
    dimnames(arc.ranks) <- list(rownames(site.data), years, phases)
  

  # matrix with agreement % values for all products, all phases 
    vi.mat <- matrix(ncol = length(phases)*length(prods), nrow = length(rownames(site.data)))
    rownames(vi.mat) <- rownames(site.data)

  # not happy with this, but need to access subset of columns later on - this should be fixed 
    phaseprod1 <- paste0(phases[1],prods)
    phaseprod2 <- paste0(phases[2],prods)
    colnames(vi.mat) <- c(phaseprod1, phaseprod2)

  # create an array with ecdf rankings for each product, each site, across all years
    vi.rankings <- array(data = NA, dim = c(length(rownames(site.data)), length(years), length(prods)))
    dimnames(vi.rankings) <- list(rownames(site.data), years, prods)

  # create single matrix to hold the worst years for all ARC windows  
    event.years <- 2*badyear.thres*length(years)
    arc.worst <- matrix(data = NA, ncol = event.years, nrow = dim(site.data)[1])
    rownames(arc.worst) <- rownames(site.data)
    colnames(arc.worst) <- rep(phases, each = event.years/2) # gives "Early Early..."



  # begin comparison process on a site-by-site basis 
    for (site in rownames(site.data)){
      
      # for troubleshooting purposes 
       print(paste0("Now calculating comparisons for ",site))
  
      # read in rainfall by phase data  
        rfe <- read.csv(paste0(base.path,scen,site,rbf.fn), header = T, row.names = 1)
        rfe.sub <- rfe[years,]  #only years we're interested in 
        rfe.p1 <- rfe.sub[years,"X1"] # first phase 
        rfe.p2 <- rfe.sub[years,"X3"] # second phase 

       
      ##### Have not worried about RFE for this analysis, but base work is here #  
       
      # generate rankings 
#        Fn1 <- ecdf(rfe.p1)
#        Fn2 <- ecdf(rfe.p2) 
#        r1 <- round(Fn1(rfe.p1), digits = 2)   # ranks for 1st phase
#        r2 <- round(Fn2(rfe.p2), digits = 2)   # ranks for 2nd phase
      
      # identify the years in the bottom "badyear.thres" - problem with ties (what to do?)
#        rfe.years.early <- rownames(rfe.sub)[which(rfe.p1 <= quantile(rfe.p1,probs = badyear.thres))]   
#        rfe.years.late <- rownames(rfe.sub)[which(rfe.p2 <= quantile(rfe.p2,probs = badyear.thres))]
      
      # define sowing windows to calculate ARC rainfall by using dekad calendar 
      # days given for comparison against earlier reports 
       
       # following variables are defined on sitetable values, not individual contracts
      #  early.start <- min(which(dekad.cal == site.data[site,"EarlyFirst"]))
      #  early.end   <- max(which(dekad.cal == site.data[site,"EarlyLast"])) 
      #  late.start  <- min(which(dekad.cal == site.data[site,"LateFirst"]))
      #  late.end    <- max(which(dekad.cal == site.data[site,"LateLast"]))
       
      # these variables defined based on site contracts, following format of sniiddemo code
       
      # read in window data from site-specific contract 
       
      contract <- paste0(base.path, scen, site, "/payout.data/contract.R") 
      source(contract)
       
      early.start <- min(which(dekad.cal == (t$phases[1,1]+t$swfirst-1)))
      early.end <- max(which(dekad.cal == (t$phases[2,1]+t$swfirst-1)))
      late.start <- min(which(dekad.cal == (t$phases[1,3]+t$swfirst-1)))
      late.end <-  max(which(dekad.cal == (t$phases[2,3]+t$swfirst-1)))
       
      # give the actual MMM-DD days for comparison purposes  
       early.dates <- c(rownames(dekad.cal)[early.start], rownames(dekad.cal)[early.end])
        late.dates  <- c(rownames(dekad.cal)[late.start], rownames(dekad.cal)[late.end])
        
      # read ARC precip file in day-year format
        arc <- read.csv(paste0(base.path,scen,site,"/met.data/precip.daily.csv"), header = T, row.names = 1)
        colnames(arc) <- substr(colnames(arc),2,5)   # remove X from character string preceding colnames
      
      # create subset of only arc data over years in "years", take sum over window timings 
        arc.sub   <- arc[,years]
        earlies   <- colSums(arc.sub[early.start:(early.end),])
        lates     <- colSums(arc.sub[late.start:(late.end),]) 
      
      # identify which years satisfy the quantile requirement, NA as a flag if  
        arc.years.early <- colnames(arc.sub)[which(earlies <= quantile(earlies, probs = badyear.thres))]
       
      # since I modify the rowname, need an alternative way to call the rowname value 
       rowname.num <- which(rownames(arc.worst) == site)
       if (length(arc.years.early) != badyear.thres*length(years)){
         # coerce into being the right vector length, raise flag by changing site name
         rownames(arc.worst)[rowname.num] <- paste0(rownames(arc.worst)[rowname.num],"E",length(arc.years.early))
         arc.years.early <- arc.years.early[1:(badyear.thres*length(years))] 
       }
       
        arc.years.late  <- colnames(arc.sub)[which(lates <= quantile(lates, probs = badyear.thres))]
       if (length(arc.years.late) != badyear.thres*length(years)){
         # coerce into being the right vector length, raise flag by changing site name
         rownames(arc.worst)[rowname.num] <- paste0(rownames(arc.worst)[rowname.num],"L",length(arc.years.late))
         arc.years.late <- arc.years.late[1:(badyear.thres*length(years))] 
       }     
       
       
       
      # write-in the ecdf rankings of arc estimates
        arc.Fn <- ecdf(earlies)
        arc.ranks[site,,"Early"]  <- round(arc.Fn(earlies), digits = 2) 
      
        arc.Fn <- ecdf(lates)
        arc.ranks[site,,"Late"]   <- round(arc.Fn(lates), digits = 2)  
       
      # returns vector of rankings agreement (overlap) results for early window
      # uses "vi.compare" from vi_functions script  
      # vi.mat includes both early and late rankings, whereas other objs don't   
        phase <- "Early"
        print(sapply(prods,vi.compare))
        vi.mat[site,phaseprod1] <- sapply(prods, vi.compare) 
       
        phase <- "Late"
        print(sapply(prods,vi.compare))
        vi.mat[site,phaseprod2] <- sapply(prods, vi.compare)     
      
      # converting matrix to array for incorporation into master array
      # write ecdf rankings values to master array
        ranks.mat <- sapply(prods, vi.ecdf)
        dim(ranks.mat) <- c(1, length(years), length(prods)) #individual sites
        dimnames(ranks.mat) <- list(site,years,prods)
        vi.rankings[site,,]  <-  ranks.mat
    
        arc.worst[rowname.num,] <- as.numeric(c(arc.years.early,arc.years.late))         
       
  } # end of across sites loop 

  # write worst ARC years by window to file 
    f.out <- paste0(out.path,"ARC_worst_years",years[1],"-", years[length(years)],".csv")
    write.csv(arc.worst, f.out) 

  # write single-phase, product-level ecdf rankings to csv
  # how could i approach this using a non-for loop approach?
    f.out <- paste0(out.path,prods,years[1],"-",years[length(years)],"rankings.csv")
    for (i in 1:length(prods)){
      write.csv(vi.rankings[,,i], f.out[i])
      }
    # comparison of these against existing tables looks good  
  
  # create all-site rankings table with lat/lon coordinates 
  # need to add lat/lon & site names as separate entry since previously only used as rownames 
    arc.early <-cbind(arc.ranks[,,"Early"],site.data[,c("Latitude","Longitude")],rownames(site.data))
      colnames(arc.early)[(length(years)+3)] <- "site"
      rownames(arc.early) <- NULL
    arc.late <- cbind(arc.ranks[,,"Late"],site.data[,c("Latitude","Longitude")], rownames(site.data))
      colnames(arc.late)[(length(years)+3)] <- "site"
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

  # API key from Cloud Made Maps
    cm.key <- "3ba6f5c05bc142209d423981fcbacb4a"  

  # create map scenes that have site label text atop map layer, split into N and S 
  # first divide into two scenes, since much of the land in the larger bounding box 
  # does not include any sites 
    s.sites <- subset(site.data, Latitude < 9, select = c(Latitude, Longitude, Woreda))
    n.sites <- subset(site.data, Latitude > 12, select = c(Latitude, Longitude, Woreda))
  
  # create bounding box coordinates, specify buffer size, function in vi_functions.R
    s.coords <- make.coords(s.sites, 0.1)
    n.coords <- make.coords(n.sites, 0.1)

  # read in cloudmade maps given bounding information from make.coords functions
    s.map <- get_cloudmademap(bbox = c(left = s.coords$l, bottom = s.coords$b, right = s.coords$r, top = s.coords$t), api_key = cm.key)  
    n.map <- get_cloudmademap(bbox = c(left = n.coords$l, bottom = n.coords$b, right = n.coords$r, top = n.coords$t), api_key = cm.key)  

  # maps of site labels, saved to output folder 
    ggmap(s.map) + geom_text(aes(x=Longitude, y=Latitude, label=rownames(s.sites)), data = s.sites, size = 3) + labs(title = "Ethiopia Harita Sites (S)")
    ggsave(filename = paste0(out.path,"EthSouthSitesMap.pdf"))

    ggmap(n.map) + geom_text(aes(x=Longitude, y=Latitude, label=rownames(n.sites)), data = n.sites, alpha = 1.0, size = 3) + labs(title = "Ethiopia Harita Sites (N)")
    ggsave(filename = paste0(out.path,"EthNorthSitesMap.pdf"))

  # maps of woreda labels, saved to output folder 
    ggmap(s.map) + geom_text(aes(x=Longitude, y=Latitude, label=Woreda), data = s.sites, alpha = 1.0, size = 5) + labs(title = "Ethiopia Harita Woredas (S)")
    ggsave(filename = paste0(out.path,"EthSouthWoredasMap.pdf"))

    ggmap(n.map) + geom_text(aes(x=Longitude, y=Latitude, label=Woreda), data = n.sites, alpha = 1.0, size = 3) + labs(title = "Ethiopia Harita Woredas (N)")
    ggsave(filename = paste0(out.path,"EthNorthWoredasMap.pdf"))



# set bounding box coordinates for map background. 
# May consider modifying to exclude Michael Debir for crisper presentation. 
# Problems arise when df points are outside the bounding box grabbed from OSM, 
# therefore ensure they are inside.





  b.s <- 0.2  # buffer size in degrees, for how large to create the bounding box around the selected bounding box edge pixels  

override.box <- c("L","B","R","T")    # active only when override = TRUE, vector sequence is left, bottom, right, and top pixels (L,R: longitude, T,B: latitude) - must replace with values in degrees 
    
    override <- FALSE   # if TRUE, overrides default calculation of bounding box 
                        #dimensions in the viz function 

 
    if(override){
      l <- override.box[1]; b <- override.box[2]; r <- override.box[3]; t <- override.box[4]
    } else{
      # default bounding box coordinates, including the earlier specified b.s
eth.coords <- make.coords(arc.early, b.s)
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
    outmap <- get_cloudmademap(bbox = c(left = eth.coords$l, bottom = eth.coords$b, right = eth.coords$r, top = eth.coords$t), api_key = cm.key)  

  # generate early window ARC rankings, can modify number of rows output appears in  
    ggmap(outmap) + geom_point(aes(x= Longitude, y = Latitude, color = Ranking), data = a_e) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable, nrow = 1) + labs(title = "ARC Early Window Rankings")
    ggsave(file = paste0(out.path, "ARCearlyranks.pdf"), scale = 1.5)

  # generate late ARC ranks 
  # wrap in print() to display plot on-screen  
    ggmap(outmap) + geom_point(aes(x=Longitude, y=Latitude, color = Ranking), data = a_l) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable, nrow = 1) + labs(title = "ARC Late Window Rankings")
    ggsave(file = paste0(out.path,"ARClateranks.pdf"), scale = 1.5)
    

#Comparing worst years of ARC2 against select VI products.  Results should 
# match output in earlier ILO reports.


# ====  These functions should work, but encountering problems =====#
  # Have opted to pull out the commands outside the function - can improve for later iteration 
#  arc.early <- df.melt(arc.early)  # creates a melted version for gg
#  arc.late <- df.melt(arc.late)  # creates a melted version for gg
#  arc.vi.vis(arc.early)
#  arc.vi.vis(arc.late)



# ================= Comparing VI and ARC ranks =================================
# 
# Also used for validating against previous ILO reports, this time comparing 
# individual VI products against ARC and with the badyear.thres value defining
# which are the worst years 
#
# *****************************************************************************#

  # ARC-VI Agreement on Worst Years
    # this can serve as the benchmark against which other versions can be compared to evaluate any performance gains from regridding, spatial correlation, etc.  
    ggmap(outmap) + geom_point(aes(x=Longitude, y=Latitude, color = value), data = vi.2) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable, nrow = 1 ) + labs(title = "ARC-VI Worst Year Agreement %")
    ggsave(file = paste0(out.path, "ARCEVIagree.pdf"), scale = 1.5)









# working version of point-based output 
# print(ggmap(Eth, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = gg.arc.early)) + geom_point(aes(color = factor(Ranking), size = 3)) + facet_wrap(~ variable))  

# working version using scale gradient 
#print(ggmap(Eth, base_layer = ggplot(aes(x=Longitude, y=Latitude), data = gg.arc.early)) + geom_point(aes(color = Ranking)) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ variable))  


# ======================= VI Grid Size Scaling ================================= 
#
# An alternative to comparing VI results scaled at the same pixel size as ARC2 
# is to 1.) compare across a smaller aggregate pixel size (under 10k x 10k), 
# and 2.) within that smaller box, to only compare pixels whose correlation 
# values 
#
# Since much of the code is common to the correlation work below, outputs are
# in that chunk of code.  
# 
# Will create a .csv output file titled with diameter value for comparison  
#
# How many km wide do you want the aggregate gridding box?  [0 < x <= 10 km]
#
  #  diameter.range <- 5  
    diameter.range <- 1:5
    range.number <- 10    # do i want this?  number of values to be computed across the diameter
#
#==============================================================================#

  # default pixel size in degrees, from IRI-DL
    MODIS.pixel.size  <-  0.00221704 
    SPOT.pixel.size   <-  0.008928572 

    long.eq <- 111.32 # number of kilometers between longitude degrees at equator, 
                      # if we decide to make adjustments to account for 
                      # shrinking pixel sizes when moving towards poles 

    diameter.range <- diameter.range / long.eq  # returns input argument in deg
  









#======================= ARC - VI Spatial Correlation ==========================

#
# Uses a correlation threshold 'r' such that all pixels whose ARC-VI correlation 
# coefficients exceeding 'r' inside the specified bounding box are averaged 
# out, generating a single VI value for each time period.  
#

  # specify the correlation coefficient threshold for masking out pixels 
        corr.value <- 0.7

  # specify years of interest - uncomment 1st to reuse years from beginning of analysis 
      # years <- years   
        years <- as.character(seq(2002,2012))  
  
#=============================================================================#



  # create template df from which early/late windows will be drawn - final agreement % values will be stored here as we move through the loop 
    agree.df <- data.frame(Latitude = site.data[,"Latitude"], Longitude = site.data[,"Longitude"], Agreement = NA)
    rownames(agree.df) <- rownames(site.data)
    
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
  # delay already inherent in IRI Data Library code in evi.corr.regrid since using
  # shiftdatashort with a 1 month lag 
  dekdelay <- 0 
  earlymonth <- as.character(dekadmonth[(midearly+dekdelay)%%36,"Month"])
  latemonth <- as.character(dekadmonth[(midlate+dekdelay)%%36,"Month"])
  
  
   
  
  
  # PSEUDO - this is where the rescaled outputs will be generated 
  month <- earlymonth
 # evi.scale.early <- evi.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Month = earlymonth, Size = )  

  
  # generates lagged EVI values for all years available   
  evi.early <- data.frame(evi.corr.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Size = b.s, CorrThreshold = corr.value, Month = earlymonth), Window = "Early") 
  
  evi.late <- data.frame(evi.corr.regrid(Lat = site.data[site,"Latitude"], Lon = site.data[site,"Longitude"], Size = b.s, CorrThreshold = corr.value, Month = latemonth), Window = "Late") 
  
  # now subset to the years of interest and identify which years agree with ARC worst
  # early.worst / late.worst produce worst years, while *.ranks give total ranks vectors
    evi.early <- data.frame(evi.early, Month = substr(evi.early$Time,1,3)[1], Year = as.numeric(substr(evi.early$Time, 5, 8)))
    evi.late  <- data.frame(evi.late,  Month = substr(evi.late$Time,1,3)[1], Year = as.numeric(substr(evi.late$Time, 5, 8)))
  
  # drop the original "Time" column 
    evi.early <- subset(evi.early, select = -Time)
   # evi.late  <- subset(evi.late, Year == years, select = -Time)
    evi.late  <- subset(evi.late, select = -Time)
  
  # PSEUDO - create code to ensure that all years are providing data for the same month

    evi.early <- evi.early[is.element(evi.early$Year,years),]
    evi.late  <- evi.late[is.element(evi.late$Year,years),]
  
    Fn <- ecdf(evi.early$correlation) # since under the function, IRIDL generates "correlation" 
    early.ranks <- round(Fn(evi.early$correlation),digits = 2)
    early.worst <- evi.early$Year[which(early.ranks <= badyear.thres)]
  
    Fn <- ecdf(evi.late$correlation) # since under the function, IRIDL generates "flag" 
    late.ranks <- round(Fn(evi.late$correlation),digits = 2)
    late.worst <- evi.late$Year[which(late.ranks <= badyear.thres)]
    
  print(paste0("In ", site, " the worst ", round(badyear.thres, digits = 2), " early window years for the [", years[1], ",", years[length(years)], "]", " period are ", early.worst[1], ",", early.worst[2]))
  
  print(paste0("In ", site, " the worst ", round(badyear.thres, digits = 2), " late window years for the [", years[1], ",", years[length(years)], "]", " period are ", late.worst[1], ",", late.worst[2])) 
  
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

  # plot title 
    ti <- paste0("ARC-EVI Agreement for Pixels with R > ", corr.value)
  # plot filename  
    fn <- paste0(out.path, "VIspatialcorrelation_",corr.value,"_.pdf")

  ggmap(outmap) + geom_point(aes(x=Longitude, y=Latitude, color = Agreement), data = arc.evi) + scale_color_gradient(low = "red", high = "green") + facet_wrap(~ Window, nrow = 1) + labs(title = ti)    
  ggsave(file = fn, scale = 1.5)
