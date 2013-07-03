ILO
===

ILO_Materials


Tasks 
=== 
* Determine if reason for discrepancy against earlier tables is because of imported contract parameters 
* Can save the downloaded EVI files to an output folder and create a flag to avoid lengthy downloads 
* Save ggplot figures as external PDFs - takes too long to load as internal plots 
* Better integration between evi.regrid and evi.corr.regrid is needed - further work needed to operationalize evi.regrid
* Create an "agree" function which is common to several sections of code - can leverage on the existence of the vi.compare function 
* Plots for scaled VI box - agreement with worst years 
* Review the code and determine what can be turned into functions 
* Work on stage 2 of grid scaling: first I work on getting the total grid to function properly
* Identify which windows were used for calculating early/late ARC2 in Phase One report 
* Testing MODIS data under new IRI-DL servers
* Presumable scoping problems in grabbing df in the `arc.vi.vis' function in "vi_functions.R" - have searched the S/O forums but not found an answer that works
* Consider working with stamen maps, code example below 


Completed Tasks
=== 
* Send brief write-up of Chris/Marshall's work to Chris for review - 7/2 - awaiting response from Marshall  



Report Text
===
* Currently maintained in Google docs 


Code Examples
=== 
Stamen Maps: map <- get_stamenmap(bbox = c(left = lon.r[1], bottom = (lat.r[1]+3), top = (lat.r[2]-1), right = lon.r[2]), zoom = 12, maptype = "toner")
ggmap(map) + geom_point(aes(x= Longitude, y = Latitude, size = RVar), data = vi.mat, colour = 'red') 
