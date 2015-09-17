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
outRaster = arcpy.GetParameterAsText(3)

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
arcpy.AddMessage("Calculating average vapor pressures")
arcpy.MakeTableView_management(data_table, "table1", "vapor_pressure > -500")
arcpy.Statistics_analysis("table1", "in_memory/tableAverage", "vapor_pressure MEAN", "site_key")
arcpy.JoinField_management(fcStations_wElevation, "Site_Key", "in_memory/tableAverage", "site_key", "MEAN_vapor_pressure")

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
ols = arcpy.OrdinaryLeastSquares_stats(fcStations_wElevation,"Unique_ID","in_memory/fcResid","MEAN_vapor_pressure","RASTERVALU", scratchGDB + "/coef_table","","")
#arcpy.TableToTable_conversion(scratchGDB + "/coef_table", "in_memory", "coef_table")
#arcpy.Delete_management(scratchGDB + "/coef_table")
intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]

#Estimate air temperature at all stations using elevation, and slope and intercept from above
arcpy.AddMessage("Estimating vapor pressure at all stations")
cursor = arcpy.UpdateCursor(fcStations_wElevation)
for row in cursor:
    if str(row.getValue("MEAN_vapor_pressure")) == "None":
        row.setValue("MEAN_vapor_pressure", (row.getValue("RASTERVALU") * slope) + intercept)
        cursor.updateRow(row)
del cursor
del row

#Run empirical bayesian kriging
arcpy.AddMessage("Running Empirical Bayesian Kriging")
arcpy.EmpiricalBayesianKriging_ga(in_features=fcStations_wElevation, z_field="MEAN_vapor_pressure", out_ga_layer="#", \
    out_raster=outRaster, cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", overlap_factor="1", \
    number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
    output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")
