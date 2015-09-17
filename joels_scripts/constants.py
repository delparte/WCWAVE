#Import necessary modules
import arcpy
from arcpy.sa import *
import numpy

#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
conRL = arcpy.GetParameter(1)
conRL_ouRaster = arcpy.GetParameterAsText(2)
conH2OSat = arcpy.GetParameter(3)
conH2OSat_outRaster = arcpy.GetParameterAsText(4)

#Set up workspace
scratchWS = arcpy.env.scratchWorkspace
scratchGDB = arcpy.env.scratchGDB
#output cell size and processing extent should be the same as elevation raster
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
arcpy.env.extent = elevation_raster
extent = arcpy.env.extent
arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "75%"
arcpy.Delete_management("in_memory")

#Get coordinate system information
desc = arcpy.Describe(elevation_raster)
coordSystem = desc.spatialReference


arcpy.AddMessage("Creating constant roughness length raster")
rlConstant = CreateConstantRaster(conRL, "FLOAT", output_cell_size, extent)
arcpy.DefineProjection_management(rlConstant, coordSystem)
rlConstant.save(conRL_ouRaster)

arcpy.AddMessage("Creating constant liquid water saturation raster")
waterConstant = CreateConstantRaster(conH2OSat, "FLOAT", output_cell_size, extent)
arcpy.DefineProjection_management(waterConstant, coordSystem)
waterConstant.save(conH2OSat_outRaster)