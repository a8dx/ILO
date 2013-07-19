ILO July 2013 Report 
===



Tasks 
=== 
* Now squelcing resizing (i.e. Box average not occurring 
* Agreement with worst years according to Ethiopian years in site table, adding 7 years to each  
* Run across multiple correlation values then automate procedure of reviewing across correlation values - consolidate results into a single dataframe that includes correlation value as a column  - need to make graphic of this 
* Create an "agree" function which is common to several sections of code - can leverage on the existence of the vi.compare function 
* Presumable scoping problems in grabbing df in the `arc.vi.vis' function in "vi_functions.R" - have searched the S/O forums but not found an answer that works
* Consider working with stamen maps, code example below 


Completed Tasks
=== 
* Provide NA value if have NA for worst years?  Now defaults to 0 which we don't want.  
* Create a matrix and visual output that would show the XY average of each site's correlation 
* Should include b.s in the file naming convention of the spatial correlation subsection.  
* Need to fix the `grid' command in Ingrid code to boxAverage 
* Created a new .csv that contains all worst years to be called by spatial correlation section 
* Send brief write-up of Chris/Marshall's work to Chris for review - 7/2 
* Modify ggmap output to not deal with baselayer - use same syntax as for woreda/site maps 



Report Text
===
* Currently maintained on Github 


Code Examples
=== 
Stamen Maps: map <- get_stamenmap(bbox = c(left = lon.r[1], bottom = (lat.r[1]+3), top = (lat.r[2]-1), right = lon.r[2]), zoom = 12, maptype = "toner")
ggmap(map) + geom_point(aes(x= Longitude, y = Latitude, size = RVar), data = vi.mat, colour = 'red') 
