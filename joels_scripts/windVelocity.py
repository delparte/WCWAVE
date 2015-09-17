#Import necessary modules
import arcpy
import datetime
import os
import subprocess
from arcpy.sa import *

#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
date_time = arcpy.GetParameterAsText(1)
station_file = arcpy.GetParameterAsText(2)
outRaster = arcpy.GetParameterAsText(3)

ninjaPath = "C:/WindNinja/WindNinja-2.5.2/bin/WindNinja_cli.exe"
if not os.path.exists(ninjaPath):
    ninjaPath = "C:/WindNinja/WindNinja-2.4.0/bin/WindNinja_cli.exe"

#Setup workspace
#output cell size should be the same as elevation raster cell size
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
scratchGDB = arcpy.env.scratchGDB
arcpy.env.overwriteOutput = True
sr = arcpy.Describe(elevation_raster).spatialReference

arcpy.AddMessage("Setting up WindNinja parameters")
#parse "date_time" string to fill in values in args below
dateParts = date_time.split("/")
strMonth = dateParts[0]
strDay = dateParts[1]
dateParts2 = dateParts[2].split(" ")
strYear = dateParts2[0]
timeParts = dateParts2[1].split(":")
strHour = timeParts[0]
strMinute = timeParts[1]
if dateParts2[2] == "PM" and strHour != "12":
    intHour = int(strHour) + 12
    strHour = str(intHour)
elif dateParts2[2] == "AM" and strHour == "12":
    strHour = str(00)

#"args" lists the WindNinja parameters that are required for the program to run
args = [ninjaPath,
"--initialization_method", "pointInitialization",
"--num_threads", "4", #number of threads to use in model run
"--elevation_file", elevation_raster, #elevation raster (cannot contain any "no-data" values)
"--match_points", "false", #match simulations to points (simulation fails if set to true)
"--year", strYear,
"--month", strMonth,
"--day", strDay,
"--hour", strHour,
"--minute", strMinute,
"--mesh_resolution", output_cell_size, #Resolution of model calculations
"--vegetation", "brush", #Vegetation type (can be 'grass', 'brush', or 'trees')
"--time_zone", "America/Boise", #time zone of target simulation
"--diurnal_winds", "true", #consider diurnal cycles in calculations
"--write_goog_output", "false", #write kml output (boolean: true/false)
"--write_shapefile_output", "false", #write shapefile output (boolean: true/false)
"--write_farsite_atm", "false", #write fire behavior file (boolean: true/false)
"--write_ascii_output", "true", #write ascii file output (this should always be set to true)
"--units_mesh_resolution", "m", #units of resolution of model calculations (should be "m" for meters)
"--units_output_wind_height", "m", #units of output wind height
"--output_speed_units", "mps",
"--output_wind_height", "3",
"--wx_station_filename", station_file] #weather station csv file used in point initialization method

arcpy.AddMessage("elevation file: " + elevation_raster)
arcpy.AddMessage("day: " + strDay)
arcpy.AddMessage("month: " + strMonth)
arcpy.AddMessage("year: " + strYear)
arcpy.AddMessage("hour: " + strHour)
arcpy.AddMessage("minute: " + strMinute)
arcpy.AddMessage("resolution: " + output_cell_size)

if len(strDay) < 2:
    strDay ="0" + strDay
if len(strMonth) < 2:
    strMonth ="0" + strMonth

path2ASCII = elevation_raster[:-4] + "_point_" + strMonth + "-" + strDay + "-" + strYear + "_" + strHour + strMinute + "_" + output_cell_size + "m" +"_non_neutral_stability_vel.asc"
arcpy.AddMessage("ASCII Path: " + path2ASCII)

#run the WindNinja_cli.exe (output is written to same location as elevation raster)
arcpy.AddMessage("Calling WindNinja command line interface")
runfile = subprocess.Popen(args, stdout = subprocess.PIPE, bufsize = -1)
runfile.wait()
output = runfile.stdout.read()
if output is None:
    arcpy.AddMessage("Results: None returned\n")
else:
    arcpy.AddMessage("Results:\n" + output)

#If the above execution was successful, there should be a "*_vel.asc" file in the same directory
#as the elevation raster. Get a handle on this file and create a new grid with the appropriate
#projection

#convert ascii file to new grid
arcpy.AddMessage("Converting WindNinja output to new raster")
arcpy.ASCIIToRaster_conversion(in_ascii_file=path2ASCII, out_raster=outRaster, data_type="FLOAT")

arcpy.DefineProjection_management(outRaster, sr)










