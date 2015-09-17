#Import necessary modules
import arcpy
from arcpy.sa import *
import numpy
import os, sys
import json
from scipy import stats



#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
ll_interp_values = json.loads(arcpy.GetParameter(1).JSON)
ul_interp_values = json.loads(arcpy.GetParameter(2).JSON)
density_interp_values = json.loads(arcpy.GetParameter(3).JSON)
outDensity = arcpy.GetParameterAsText(4)
outULTemp = arcpy.GetParameterAsText(5)
outSnowTemp = arcpy.GetParameterAsText(6)

#Set up workspace
scratchWS = arcpy.env.scratchWorkspace
scratchGDB = arcpy.env.scratchGDB
#output cell size and processing extent should be the same as elevation raster
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
arcpy.env.extent = elevation_raster
arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "75%"
arcpy.Delete_management("in_memory")

if len(density_interp_values['features']) <= 1:
    #Density Equation: y = -0.0395(elevation) + 405.26
    arcpy.AddMessage("Creating snow density raster from linear interpolation")
    snow_density_raster = -0.0395 * Raster(elevation_raster) + 405.26
    snow_density_raster.save(outDensity)
else:
    arcpy.AddMessage("Creating snow density raster from linear interpolation")
    lsElevation = []
    lsDensity = []
    for rec in density_interp_values['features']:
        lsElevation.append(rec['attributes']['Elevation'])
        lsDensity.append(rec['attributes']['Density'])

    lr_results = stats.linregress(lsElevation,lsDensity)
    slope1 = lr_results[0]
    intercept1 = lr_results[1]
    snow_density_raster = slope1 * Raster(elevation_raster) + intercept1
    snow_density_raster.save(outDensity)

if len(ul_interp_values['features']) <= 1:
    #Upper Layer Temperature Equation: y = -0.0008(elevation) + 0.1053
    arcpy.AddMessage("Creating upper snow layer temperature raster from linear interpolation")
    upper_layer_temperature = -0.0008 * Raster(elevation_raster) + 0.1053
    upper_layer_temperature.save(outULTemp)
else:
    arcpy.AddMessage("Creating upper snow layer temperature raster from linear interpolation")
    lsElevation = []
    lsTemperature = []
    for rec in ul_interp_values['features']:
        lsElevation.append(rec['attributes']['Elevation'])
        lsTemperature.append(rec['attributes']['Temperature'])

    lr_results = stats.linregress(lsElevation,lsTemperature)
    slope2 = lr_results[0]
    intercept2 = lr_results[1]
    upper_layer_temperature = slope2 * Raster(elevation_raster) + intercept2
    upper_layer_temperature.save(outULTemp)

if len(ll_interp_values['features']) <= 1:
    #lower layer temperature equation: y = -0.0008(elevation) + 1.3056
    arcpy.AddMessage("Creating lower snow layer temperature raster from linear interpolation")
    lower_layer_temperature = -0.0008 * Raster(elevation_raster) + 1.3056
else:
    arcpy.AddMessage("Creating lower snow layer temperature raster from linear interpolation")
    lsElevation = []
    lsTemperature = []
    for rec in ul_interp_values['features']:
        lsElevation.append(rec['attributes']['Elevation'])
        lsTemperature.append(rec['attributes']['Temperature'])

    lr_results = stats.linregress(lsElevation,lsTemperature)
    slope3 = lr_results[0]
    intercept3 = lr_results[1]
    lower_layer_temperature = slope2 * Raster(elevation_raster) + intercept2

#average snowcover temperature is the average of the upper and lower layer temperatures
arcpy.AddMessage("Creating average snowcover temperature raster")
average_snowcover_temperature = arcpy.sa.CellStatistics([upper_layer_temperature, lower_layer_temperature],"MEAN","NODATA")
average_snowcover_temperature.save(outSnowTemp)















