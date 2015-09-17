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

#Copy station_locations
# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "station_locations"
arcpy.FeatureClassToFeatureClass_conversion(in_features="station_locations", out_path=scratchGDB, out_name="temp_stations3", where_clause="", field_mapping="""Site_Key "Site_Key" true true false 254 Text 0 0 ,First,#,station_locations,Site_Key,-1,-1;Climate "Climate" true true false 254 Text 0 0 ,First,#,station_locations,Climate,-1,-1;Soil_Micro "Soil_Micro" true true false 254 Text 0 0 ,First,#,station_locations,Soil_Micro,-1,-1;Precipitat "Precipitat" true true false 254 Text 0 0 ,First,#,station_locations,Precipitat,-1,-1;Streamflow "Streamflow" true true false 254 Text 0 0 ,First,#,station_locations,Streamflow,-1,-1;Eddy_Covar "Eddy_Covar" true true false 254 Text 0 0 ,First,#,station_locations,Eddy_Covar,-1,-1;Snow_Cours "Snow_Cours" true true false 254 Text 0 0 ,First,#,station_locations,Snow_Cours,-1,-1;Site_Name "Site_Name" true true false 254 Text 0 0 ,First,#,station_locations,Site_Name,-1,-1;lat_wgs84 "lat_wgs84" true true false 8 Double 0 0 ,First,#,station_locations,lat_wgs84,-1,-1;lon_wgs84 "lon_wgs84" true true false 8 Double 0 0 ,First,#,station_locations,lon_wgs84,-1,-1""", config_keyword="")

#Calculate soil temperature averages over the 3 hour time period (ignoring "no-data" values: -999 etc.)
#and join to temp_stations
arcpy.AddMessage("Calculating average soil temperature")
arcpy.MakeTableView_management(data_table, "table1", "st005 > -500")
arcpy.Statistics_analysis("table1", "in_memory/tableAverage", [["st005", "MEAN"], ["elevation", "MEAN"]], "site_key")
arcpy.JoinField_management(scratchGDB + "/temp_stations3", "Site_Key", "in_memory/tableAverage", "site_key", ["MEAN_st005", "MEAN_elevation"])

#Delete rows from stations_wElevation that have negative or null elevations
arcpy.AddMessage("Removing \"no-data\" rows from station locations (eg. negative elevations)")
cursor = arcpy.UpdateCursor(scratchGDB + "/temp_stations3")
for row in cursor:
    if row.getValue("MEAN_elevation") < 0 or row.getValue("MEAN_elevation") == "None":
        cursor.deleteRow(row)
del cursor
del row

#Add unique ID field to stations_wElevation for use in OLS function
arcpy.AddField_management(scratchGDB + "/temp_stations3", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
arcpy.CalculateField_management(scratchGDB + "/temp_stations3","Unique_ID","!OBJECTID!","PYTHON_9.3","")

#Run ordinary least squares on stats1
arcpy.AddMessage("Running ordinary least squares")
arcpy.CreateTable_management(scratchGDB, "coef_table")
ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/temp_stations3","Unique_ID","in_memory/fcResid","MEAN_st005","MEAN_elevation", scratchGDB + "/coef_table","","")
intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]

#Create final raster
output_raster = (Raster(elevation_raster) * slope + intercept)
output_raster.save(outRaster)