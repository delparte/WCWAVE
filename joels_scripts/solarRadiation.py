#Import necessary modules
import arcpy
from arcpy.sa import *
import numpy
import os

#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
station_locations = arcpy.GetParameterAsText(1)
data_table = arcpy.GetParameterAsText(2)
date_time = arcpy.GetParameter(3)
out_raster = arcpy.GetParameterAsText(4)

#Set up workspace
#output cell size should be the same as elevation raster cell size
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
scratchGDB = arcpy.env.scratchGDB
arcpy.env.overwriteOutput = True

#RUN SIMULATED CLEAR-SKY CALCULATIONS
#make a copy of the station location feature class in the scratchspace
arcpy.CopyFeatures_management(station_locations, scratchGDB + "/station_locations_sr")

#set up area solar radiation tool parameters and run the tool
arcpy.AddMessage("Running area solar radiation tool")
global_radiation_raster = scratchGDB + "/global_radiation_raster"
skySize = 200
arcpy.gp.AreaSolarRadiation_sa(elevation_raster, global_radiation_raster, "",skySize, date_time,"","1","NOINTERVAL","1","FROM_DEM","32","8","8","UNIFORM_SKY","0.3","0.5","#","#","#")

#CORRECT SIMULATED VALUES TO OBSERVED DATA
arcpy.AddMessage("Correcting simulated radiation values")

#Extract simulated global radiation values to station location feature class
arcpy.MakeFeatureLayer_management(station_locations, "station_locations_temp")
arcpy.gp.ExtractValuesToPoints_sa("station_locations_temp", global_radiation_raster, scratchGDB + "/simPoints", "NONE", "VALUE_ONLY")

#Join the data table to "simPoints"
arcpy.JoinField_management(scratchGDB + "/simPoints", "Site_Key", data_table, "site_key", "in_solar_radiation")

#select non-null rows
arcpy.TableSelect_analysis(scratchGDB + "/simPoints","in_memory/selectedPoints","(in_solar_radiation IS NOT NULL and in_solar_radiation > -500) and RASTERVALU > 0")

#add "ratio" field to station locations
arcpy.AddField_management("in_memory/selectedPoints","ratio","FLOAT","#","#","#","#","NULLABLE","NON_REQUIRED","#")

#calculate "ratio" field (observed radiation / simulated radiation)
arcpy.CalculateField_management("in_memory/selectedPoints","ratio","!in_solar_radiation!/ !RASTERVALU!","PYTHON_9.3","#")

#convert "ratio" field to numpy array
na = arcpy.da.TableToNumPyArray("in_memory/selectedPoints", "ratio")

#calculate average ratio
dMeanRatio = numpy.mean(na["ratio"])
dMeanRatio2 = numpy.asscalar(dMeanRatio)

#multiply simulated raster by average ratio
output_raster = Raster(global_radiation_raster) * dMeanRatio2
output_raster.save(out_raster)

