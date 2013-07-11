ILO July 2013 Report 
===

http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/.est_prcp/X/38.725/VALUE/Y/7.851/VALUE/X/removeGRID/Y/removeGRID/T/1.0/monthlyAverage/SOURCES/.USGS/.LandDAAC/.MODIS/.version_005/.EAF/.EVI/X/38.525/0.05/38.925/GRID/Y/7.651/0.05/8.051/GRID/T/1.0/monthlyAverage/T/-1/1/3/shiftdatashort/[T]/correlate/[X/Y]average/T_lag+exch+table-+text+text+-table++.csv



Tasks 
=== 
<<<<<<< HEAD
* Automate procedure of reviewing across correlation values - consolidate results into a single dataframe that includes correlation value as a column  - need to make graphic of this
* Create a matrix and visual output that would show the XY average of each site's correlation 
* Would be great to incorporate farmers' worst years from the sitetable, but it appears that the latest 'worst' year is 2002 - this cannot be usable from an EVI perspective.  
=======
* Speak with Dan about how multiple years with same precip value were handled in rankings - e.g. MayTeklia with early window tie of 3 
>>>>>>> b25e53a17f70995378837642622062fe1f369e80
* Can create a way to ensure that number of plots on a single row is not excessive, instead of standard nrow = 1 argument
* Identify why spatial correlation late window is not working.. 
* Determine if reason for discrepancy against earlier tables is because of imported contract parameters 
* Can save the downloaded EVI files to an output folder and create a flag to avoid lengthy downloads 
* Better integration between evi.regrid and evi.corr.regrid is needed - further work needed to operationalize evi.regrid
* Create an "agree" function which is common to several sections of code - can leverage on the existence of the vi.compare function 
* Plots for scaled VI box - agreement with worst years 
* Review the code and determine what can be turned into functions 
* Work on stage 2 of grid scaling: first I work on getting the total grid to function properly
* Identify which windows were used for calculating early/late ARC2 in Phase One report 
* Presumable scoping problems in grabbing df in the `arc.vi.vis' function in "vi_functions.R" - have searched the S/O forums but not found an answer that works
* Consider working with stamen maps, code example below 


Completed Tasks
=== 
* Created a new .csv that contains all worst years to be called by spatial correlation section 
* Send brief write-up of Chris/Marshall's work to Chris for review - 7/2 
* Modify ggmap output to not deal with baselayer - use same syntax as for woreda/site maps 



Report Text
===
* Currently maintained in Google docs 


Code Examples
=== 
Stamen Maps: map <- get_stamenmap(bbox = c(left = lon.r[1], bottom = (lat.r[1]+3), top = (lat.r[2]-1), right = lon.r[2]), zoom = 12, maptype = "toner")
ggmap(map) + geom_point(aes(x= Longitude, y = Latitude, size = RVar), data = vi.mat, colour = 'red') 
