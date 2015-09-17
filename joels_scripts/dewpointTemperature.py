#Import necessary modules
import arcpy
from arcpy.sa import *
import numpy

#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
station_locations = arcpy.GetParameterAsText(1)
data_table = arcpy.GetParameterAsText(2)
outDPRaster = arcpy.GetParameterAsText(3)
outPercent = arcpy.GetParameterAsText(4)
outDensity = arcpy.GetParameterAsText(5)

#Set up workspace
scratchWS = arcpy.env.scratchWorkspace
scratchGDB = arcpy.env.scratchGDB
#output cell size and processing extent should be the same as elevation raster
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
arcpy.env.extent = elevation_raster
arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "75%"
#arcpy.Delete_management("in_memory")

#extract elevations to stations
arcpy.AddMessage("Extracting elevations")
fcStations_wElevation = scratchGDB + "/stations_wElevation"
arcpy.gp.ExtractValuesToPoints_sa(station_locations, elevation_raster, fcStations_wElevation, "NONE", "VALUE_ONLY")

#Calculate vapor pressure averages over the 3 hour time period (ignoring "no-data" values: -999 etc.)
#and join to stations_wElevation
arcpy.AddMessage("Calculating average dewpoint temperatures")
arcpy.MakeTableView_management(data_table, "table1", "dewpoint_temperature > -500")
arcpy.Statistics_analysis("table1", "in_memory/tableAverage", "dewpoint_temperature MEAN", "site_key")
arcpy.JoinField_management(fcStations_wElevation, "Site_Key", "in_memory/tableAverage", "site_key", "MEAN_dewpoint_temperature")

#Delete rows from stations_wElevation that have negative or null elevations
arcpy.AddMessage("Removing \"no-data\" rows from station locations (eg. negative elevations)")
cursor = arcpy.UpdateCursor(fcStations_wElevation)
for row in cursor:
    if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None":
        cursor.deleteRow(row)
del cursor
del row

#Add unique ID field to stations_wElevation for use in OLS function
arcpy.AddField_management(fcStations_wElevation, "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
arcpy.CalculateField_management(fcStations_wElevation,"Unique_ID","!OBJECTID!","PYTHON_9.3","")

#Run ordinary least squares on stats1
arcpy.AddMessage("Running ordinary least squares")
arcpy.CreateTable_management(scratchGDB, "coef_table")
ols = arcpy.OrdinaryLeastSquares_stats(fcStations_wElevation,"Unique_ID","in_memory/fcResid","MEAN_dewpoint_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
#arcpy.TableToTable_conversion(scratchGDB + "/coef_table", "in_memory", "coef_table")
#arcpy.Delete_management(scratchGDB + "/coef_table")
intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]

#Estimate air temperature at all stations using elevation, and slope and intercept from above
arcpy.AddMessage("Estimating dewpoint temperature at all stations")
cursor = arcpy.UpdateCursor(fcStations_wElevation)
for row in cursor:
    if str(row.getValue("MEAN_dewpoint_temperature")) == "None":
        row.setValue("MEAN_dewpoint_temperature", (row.getValue("RASTERVALU") * slope) + intercept)
        cursor.updateRow(row)
del cursor
del row

#Run empirical bayesian kriging
arcpy.AddMessage("Running Empirical Bayesian Kriging")
arcpy.EmpiricalBayesianKriging_ga(in_features=fcStations_wElevation, z_field="MEAN_dewpoint_temperature", out_ga_layer="#", out_raster=scratchGDB + "/raster3", \
    cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", overlap_factor="1", number_semivariograms="100", \
    search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", output_type="PREDICTION", quantile_value="0.5", \
    threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

arcpy.CopyRaster_management(scratchGDB + "/raster3", outDPRaster)

#Create precipitation properties rasters (% precip as snow, density of snow portion of precip)
arcpy.AddMessage("Creating precipitation properties rasters")
#Percent snow conditional:
arcpy.AddMessage("  Creating raster for percentage of precipitation that was snow")
inRas = Raster(scratchGDB + "/raster3")
outPercentSnowRaster = arcpy.sa.Con(inRas < -5.0, 1.0, Con((inRas >= -5.0) & (inRas < -3.0), 1.0,
                                                           Con((inRas >= -3.0) & (inRas < -1.5), 1.0,
                                                               Con((inRas >= -1.5) & (inRas < -0.5), 1.0,
                                                                   Con((inRas >= -0.5) & (inRas < 0.0), 0.75,
                                                                       Con((inRas >= 0.0) & (inRas < 0.5), 0.25,
                                                                           Con(inRas >= 0.5,0.0)))))))
arcpy.CopyRaster_management(outPercentSnowRaster,outPercent,pixel_type="32_BIT_FLOAT")

#Snow density conditional:
arcpy.AddMessage("  Creating raster for density of snow portion of precipitation")
outSnowDensityRaster = arcpy.sa.Con(inRas < -5.0, 75.0, Con((inRas >= -5.0) & (inRas < -3.0), 100.0,
                                                           Con((inRas >= -3.0) & (inRas < -1.5), 150.0,
                                                               Con((inRas >= -1.5) & (inRas < -0.5), 175.0,
                                                                   Con((inRas >= -0.5) & (inRas < 0), 200.0,
                                                                       Con((inRas >= 0.0) & (inRas < 0.5), 250.0,
                                                                           Con(inRas >= 0.5,0.0)))))))
arcpy.CopyRaster_management(outSnowDensityRaster,outDensity,pixel_type="32_BIT_FLOAT")
