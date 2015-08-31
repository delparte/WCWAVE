import arcpy
from arcpy.sa import *
import os
import sys
import csv
import traceback
import numpy
from scipy import stats
import mysql.connector
from mysql.connector import errorcode
import datetime
import json
import shutil
import subprocess

#Assign input parameters
sFromDate = arcpy.GetParameterAsText(0)
sToDate = arcpy.GetParameterAsText(1)
iTimeStep = int(arcpy.GetParameterAsText(2))
sKrigMethod = arcpy.GetParameterAsText(3)
boolAllTools = arcpy.GetParameter(4)
boolAirTemperature = arcpy.GetParameter(5)

boolConstants = arcpy.GetParameter(6)
dRLConstant = arcpy.GetParameter(7)
dH2OConstant = arcpy.GetParameter(8)

boolDewPoint = arcpy.GetParameter(9)
boolPrecipMass = arcpy.GetParameter(10)
boolSnowDepth = arcpy.GetParameter(11)

boolSnowProperties = arcpy.GetParameter(12)
ll_interp_values = json.loads(arcpy.GetParameter(13).JSON)
ul_interp_values = json.loads(arcpy.GetParameter(14).JSON)
density_interp_values = json.loads(arcpy.GetParameter(15).JSON)

boolSoilTemperature = arcpy.GetParameter(16)
boolSolarRadiation = arcpy.GetParameter(17)
boolThermalRadiation = arcpy.GetParameter(18)
boolVaporPressure = arcpy.GetParameter(19)
boolWindSpeed = arcpy.GetParameter(20)
#tempOutput = arcpy.GetParameterAsText(16)

#Specify workspace
scratchWS = arcpy.env.scratchFolder
arcpy.AddMessage("Scratch Workspace: " + scratchWS)
#scratchWS = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\scratch'
#Make output folder to zip
dateNow = datetime.datetime.now()
sNow = dateNow.strftime("%Y%d%b_%H%M")
os.makedirs(scratchWS + "/Output_" + sNow)
outFolder = scratchWS + "/Output_" + sNow
arcpy.AddMessage("Output Folder: " + outFolder)
scratchGDB = arcpy.env.scratchGDB
arcpy.env.overwriteOutput = True

#SQL variables
tmp3 = "ta"


#Define Functions
def roundTime(dt, roundTo=60):
    seconds = (dt - dt.min).seconds
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

def airTemperature():
    #Caclulate average air temperatures over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "airTemperature_Table", "air_temperature > -500")
    arcpy.Statistics_analysis("airTemperature_Table", "in_memory/airTemperature_Table2", "air_temperature MEAN", "site_key")
    #Copy air temperature values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/airTemperature_Table2", "site_key", "MEAN_air_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for air temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row

    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_air_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Calculate residuals and add them to tempStations
    arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        row.setValue("residual", row.getValue("MEAN_air_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
        cursor.updateRow(row)
    del cursor
    del row
    arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
        out_raster=scratchGDB + "/airTemperature_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
        overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
        output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
    lsScratchData_Imd.append(scratchGDB + "/airTemperature_residual")

    output_raster = arcpy.Raster(scratchGDB + "/airTemperature_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/air_temperature_" + sTimeStamp + ".tif")

    return outFolder + "/air_temperature_" + sTimeStamp + ".tif"

def constants():
    #Get coordinate system information
    desc = arcpy.Describe(rc_elevation)
    coordSystem = desc.spatialReference
    rlConstant = CreateConstantRaster(dRLConstant, "FLOAT", output_cell_size)
    arcpy.DefineProjection_management(rlConstant, coordSystem)
    rlConstant.save(outFolder + "/roughness_length_" + sTimeStamp + ".tif")
    waterConstant = CreateConstantRaster(dH2OConstant, "FLOAT", output_cell_size)
    arcpy.DefineProjection_management(waterConstant, coordSystem)
    waterConstant.save(outFolder + "/H2O_saturation_" + sTimeStamp + ".tif")

    return outFolder + "/roughness_length_" + sTimeStamp + ".tif", outFolder + "/H2O_saturation_" + sTimeStamp + ".tif"

def dewPoint():
    #Caclulate average dew point temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "dewPoint_Table", "dewpoint_temperature > -500")
    arcpy.Statistics_analysis("dewPoint_Table", "in_memory/dewPoint_Table2", "dewpoint_temperature MEAN", "site_key")
    #Copy dew point temperature values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/dewPoint_Table2", "site_key", "MEAN_dewpoint_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for dew point temperatures
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_dewpoint_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row

    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_dewpoint_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Calculate residuals and add them to tempStations
    arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        row.setValue("residual", row.getValue("MEAN_dewpoint_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
        cursor.updateRow(row)
    del cursor
    del row
    arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
        out_raster=scratchGDB + "/dewPoint_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
        overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
        output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
    lsScratchData_Imd.append(scratchGDB + "/dewPoint_residual")
    #Add back elevation trends and save final raster
    output_raster = arcpy.Raster(scratchGDB + "/dewPoint_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif")

    #Create precipitation properties rasters (% precip as snow, density of snow portion of precip)
    #Percent snow conditional:
    inRas = Raster(outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif")
    outPercentSnowRaster = arcpy.sa.Con(inRas < -5.0, 1.0, Con((inRas >= -5.0) & (inRas < -3.0), 1.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 1.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 1.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0.0), 0.75,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 0.25,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outPercentSnowRaster, outFolder + "/percent_snow_" + sTimeStamp + ".tif",pixel_type="32_BIT_FLOAT")

    #Snow density conditional:
    outSnowDensityRaster = arcpy.sa.Con(inRas < -5.0, 75.0, Con((inRas >= -5.0) & (inRas < -3.0), 100.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 150.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 175.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0), 200.0,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 250.0,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outSnowDensityRaster,outFolder + "/precipitation_snow_density_" + sTimeStamp + ".tif",pixel_type="32_BIT_FLOAT")


    return outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif", outFolder + "/percent_snow_" + sTimeStamp + ".tif", outFolder + "/precipitation_snow_density_" + sTimeStamp + ".tif"

def precipMass():
    #Caclulate average precipitation values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/precipitation_table1", "precipitation_Table", "ppts > -500")
    arcpy.Statistics_analysis("precipitation_Table", "in_memory/precipitation_Table2", "ppts MEAN", "site_key")
    #Copy precipitation values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/precipitation_Table2", "site_key", "MEAN_ppts")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for precipitation
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_ppts") is None:
            cursor.deleteRow(row)
    del cursor
    del row

    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_ppts","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Calculate residuals and add them to tempStations
    arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        row.setValue("residual", row.getValue("MEAN_ppts") - ((slope * row.getValue("RASTERVALU")) + intercept))
        cursor.updateRow(row)
    del cursor
    del row
    arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
        out_raster=scratchGDB + "/precipitation_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
        overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
        output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
    #Add back elevation trends and save final raster
    output_raster = arcpy.Raster(scratchGDB + "/precipitation_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/precipitation_mass_" + sTimeStamp + ".tif")

    return outFolder + "/precipitation_mass_" + sTimeStamp + ".tif"

def snowDepth():
    #Caclulate average snow depth values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/snowDepth_table1", "snowDepth_Table", "snow_depth > -500")
    arcpy.Statistics_analysis("snowDepth_Table", "in_memory/snowDepth_Table2", "snow_depth MEAN", "site_key")
    #Copy snow depth values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/snowDepth_Table2", "site_key", "MEAN_snow_depth")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for snow depth
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_snow_depth") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if sKrigMethod == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_snow_depth","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_snow_depth") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/snowDepth_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratchGDB + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratchGDB + "/snowDepth_residual")
        lsScratchData_Imd.append(scratchGDB + "/snowDepth_residual")
        #Add back elevation trends and save final raster
        imdRaster = arcpy.Raster(scratchGDB + "/snowDepth_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster = Con(imdRaster > 0, imdRaster, 0)
        output_raster.save(outFolder + "/snow_depth_" + sTimeStamp + ".tif")
    else:
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_snow_depth","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_snow_depth") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
            out_raster=scratchGDB + "/snowDepth_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
        #Add back elevation trends and save final raster
        imdRaster = arcpy.Raster(scratchGDB + "/snowDepth_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster = Con(imdRaster > 0, imdRaster, 0)
        output_raster.save(outFolder + "/snow_depth_" + sTimeStamp + ".tif")

    return outFolder + "/snow_depth_" + sTimeStamp + ".tif"

def snowProperties():
    if len(density_interp_values['features']) <= 1:
        #Density Equation: y = -0.0395(elevation) + 405.26
        snow_density_raster = -0.0395 * Raster(rc_elevation) + 405.26
        snow_density_raster.save(outFolder + "/snow_density_" + sTimeStamp + ".tif")
    else:
        lsElevation = []
        lsDensity = []
        for rec in density_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsDensity.append(rec['attributes']['Density'])

        lr_results = stats.linregress(lsElevation,lsDensity)
        slope1 = lr_results[0]
        intercept1 = lr_results[1]
        snow_density_raster = slope1 * Raster(rc_elevation) + intercept1
        snow_density_raster.save(outFolder + "/snow_density_" + sTimeStamp + ".tif")

    if len(ul_interp_values['features']) <= 1:
        #Upper Layer Temperature Equation: y = -0.0008(elevation) + 0.1053
        upper_layer_temperature = -0.0008 * Raster(rc_elevation) + 0.1053
        upper_layer_temperature.save(outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif")
    else:
        lsElevation = []
        lsTemperature = []
        for rec in ul_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsTemperature.append(rec['attributes']['Temperature'])

        lr_results = stats.linregress(lsElevation,lsTemperature)
        slope2 = lr_results[0]
        intercept2 = lr_results[1]
        upper_layer_temperature = slope2 * Raster(rc_elevation) + intercept2
        upper_layer_temperature.save(outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif")

    if len(ll_interp_values['features']) <= 1:
        #lower layer temperature equation: y = -0.0008(elevation) + 1.3056
        lower_layer_temperature = -0.0008 * Raster(rc_elevation) + 1.3056
    else:
        lsElevation = []
        lsTemperature = []
        for rec in ul_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsTemperature.append(rec['attributes']['Temperature'])

        lr_results = stats.linregress(lsElevation,lsTemperature)
        slope3 = lr_results[0]
        intercept3 = lr_results[1]
        lower_layer_temperature = slope2 * Raster(rc_elevation) + intercept2

    #average snowcover temperature is the average of the upper and lower layer temperatures
    average_snowcover_temperature = arcpy.sa.CellStatistics([upper_layer_temperature, lower_layer_temperature],"MEAN","NODATA")
    average_snowcover_temperature.save(outFolder + "/average_snow_cover_temperature_" + sTimeStamp + ".tif")

    return outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif", outFolder + "/average_snow_cover_temperature_" + sTimeStamp + ".tif", outFolder + "/snow_density_" + sTimeStamp + ".tif"

def soilTemperature():
    #Set extent to full feature class (all stations)
    arcpy.env.extent = extFullFeatures
    #Caclulate average soil temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/soilTemperature_table1", "soilTemperature_Table", "stm005 > -500")
    arcpy.Statistics_analysis("soilTemperature_Table", "in_memory/soilTemperature_Table2", "stm005 MEAN", "site_key")
    #Copy snow depth values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/soilTemperature_Table2", "site_key", "MEAN_stm005")
    #Delete rows from station feature class that have null values for soil temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_stm005") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_stm005","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Create final raster
    arcpy.env.extent = extElevation
    output_raster = (Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/soil_temperature_" + sTimeStamp + ".tif")



    return outFolder + "/soil_temperature_" + sTimeStamp + ".tif"


def solarRadiation():
    #Caclulate average solar radiation values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "solarRadiation_Table", "in_solar_radiation > -500")
    arcpy.Statistics_analysis("solarRadiation_Table", "in_memory/solarRadiation_Table2", "in_solar_radiation MEAN", "site_key")
    #Copy solar radiation values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/solarRadiation_Table2", "site_key", "MEAN_in_solar_radiation")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for solar radiation
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_in_solar_radiation") is None:
            cursor.deleteRow(row)
    del cursor
    del row


    #set up area solar radiation tool parameters and run the tool
    #Set up time parameters
    day_of_year = dateIncrement.timetuple().tm_yday
    iSRStart = int(dateIncrement.strftime("%H"))
    iSREnd = iSRStart + iTimeStep
    inTWD = TimeWithinDay(day_of_year,iSRStart,iSREnd)

    global_radiation_raster = scratchGDB + "/global_radiation_raster"
    skySize = 200
    outGlobalRadiation = AreaSolarRadiation(rc_elevation, "", skySize, inTWD)
    #outGlobalRadiation = AreaSolarRadiation(r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\RC_DEM_100m_South', "", skySize, inTWD)
    outGlobalRadiation = outGlobalRadiation / iTimeStep

    #Correct global radiation raster for cloud conditions
    #Extract simulated global radiation values to station location feature class
    arcpy.AlterField_management(scratchGDB + "/tempStations","RASTERVALU","Elevation")
    arcpy.gp.ExtractValuesToPoints_sa(scratchGDB + "/tempStations", outGlobalRadiation, scratchGDB + "/simPoints", "NONE", "VALUE_ONLY")
    lsScratchData_Imd.append(scratchGDB + "/simPoints")
    #add "ratio" field to station feature class
    arcpy.AddField_management(scratchGDB + "/simPoints","ratio","FLOAT","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    #calculate "ratio" field (observed radiation / simulated radiation)
    arcpy.CalculateField_management(scratchGDB + "/simPoints","ratio","!MEAN_in_solar_radiation!/ !RASTERVALU!","PYTHON_9.3","#")

    #convert "ratio" field to numpy array
    na = arcpy.da.TableToNumPyArray(scratchGDB + "/simPoints", "ratio")

    #calculate average ratio
    dMeanRatio = numpy.mean(na["ratio"])
    dMeanRatio2 = numpy.asscalar(dMeanRatio)

    #multiply simulated raster by average ratio
    outGlobalRadiation_Corrected = outGlobalRadiation * dMeanRatio2
    outGlobalRadiation_Corrected.save(outFolder + "/solar_radiation_" + sTimeStamp + ".tif")

    return outFolder + "/solar_radiation_" + sTimeStamp + ".tif"


def thermalRadiation(inAirTemperature, inVaporPressure, inSurfaceTemp):
    #Constants and re-defined variables (See Marks and Dozier (1979), pg. 160)
    z = rc_elevation
    vf = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\RC_ViewFactor_10m_South'
    T_a = inAirTemperature
    vp = inVaporPressure

    #Select one station and its data to use as reference for elevation, pressure, and temperature
    #Caclulate average air temperature and vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "thermalRadiation_Table", "air_temperature > -500 AND vapor_pressure > -500")
    arcpy.Statistics_analysis("thermalRadiation_Table", "in_memory/thermalRadiation_Table2", [["air_temperature", "MEAN"], ["vapor_pressure", "MEAN"]], "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/thermalRadiation_Table2", "site_key", ["MEAN_air_temperature", "MEAN_vapor_pressure"])
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure or air temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None or row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Use update cursor to extract air temperature, vapor pressure, and elevation values for the first station, then delete that station
    P_m = 0.0
    T_m = 0.0
    z_m = 0.0
    T_s = inSurfaceTemp
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        z_m = row.getValue("RASTERVALU")
        P_m = row.getValue("MEAN_vapor_pressure")
        T_m = row.getValue("MEAN_air_temperature")
        cursor.deleteRow(row)
        break

    del cursor
    del row

    arcpy.AddMessage("P_m: " + str(P_m))
    arcpy.AddMessage("T_m: " + str(T_m))
    arcpy.AddMessage("z_m: " + str(z_m))
    arcpy.AddMessage("T_s: " + str(T_s))

    g = 9.8
    m = 0.0289
    R = 8.3143
    sigma = 5.6697*10**-8
    epsilon_s = 0.95
    gamma = -0.006

    #convert temperature parameters to Kelvin
    T_m = T_m + 274.15
    T_s = T_s + 274.15
    T_a = arcpy.sa.Float(Raster(T_a) + 274.15)

    #convert vapor pressure parameters to mb
    P_m = P_m * 0.01
    vp = arcpy.sa.Float(Raster(vp) * 0.01)

    #Correct air temperature and vapor pressure rasters (Marks and Dozier (1979), pg. 164)
    T_prime = T_a + (0.0065 * Raster(rc_elevation)) #(4) corrected air temperature
    e_sa = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from original air temperature (T_a)
    e_sprime = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from corrected air temperature (T_prime)
    rh = arcpy.sa.Float(vp / e_sa) #(5) relative humidity
    e_prime = arcpy.sa.Float(rh * e_sprime) #(6) corrected vapor pressure

    #Pressure at a given elevation (Marks and Dozier (1979), pg. 168-169)
    term1 = ((-g*m)/(R*gamma))
    delta_z = Raster(z) - z_m
    term2 = ((T_m + gamma * delta_z)) / T_m
    lnTerm = arcpy.sa.Ln(term2)
    expTerm = arcpy.sa.Exp(term1 * lnTerm)
    P_a = P_m * expTerm #(10) air pressure

    #effective emissivity (Marks and Dozier (1979), pg. 164)
    epsilon_a = arcpy.sa.Float((1.24 * (e_prime / T_prime)**(1/7)) * (P_a / 1013.0)) #(7)

    #Incoming longwave radiation (Marks and Dozier (1979), pg. 164)
    term3 = arcpy.sa.Float((epsilon_a * sigma * (T_a ** 4)) * vf)
    term4 = arcpy.sa.Float(epsilon_s * sigma * (T_s ** 4))
    term5 = (1 - Raster(vf))
    output_thermal_radiation = arcpy.sa.Float(term3 + (term4 * term5)) #(9)
    output_thermal_radiation.save(outFolder + "/thermal_radiation_" + sTimeStamp + ".tif")

    return outFolder + "/thermal_radiation_" + sTimeStamp + ".tif"

def vaporPressure():
    #Caclulate average vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "vaporPressure_Table", "vapor_pressure > -500")
    arcpy.Statistics_analysis("vaporPressure_Table", "in_memory/vaporPressure_Table2", "vapor_pressure MEAN", "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/vaporPressure_Table2", "site_key", "MEAN_vapor_pressure")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row

    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_vapor_pressure","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Calculate residuals and add them to tempStations
    arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        row.setValue("residual", row.getValue("MEAN_vapor_pressure") - ((slope * row.getValue("RASTERVALU")) + intercept))
        cursor.updateRow(row)
    del cursor
    arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
        out_raster=scratchGDB + "/vaporPressure_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
        overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
        output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
    lsScratchData_Imd.append(scratchGDB + "/vaporPressure_residual")
    #Add back elevation trends and save final raster
    output_raster = arcpy.Raster(scratchGDB + "/vaporPressure_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/vapor_pressure_" + sTimeStamp + ".tif")

    return outFolder + "/vapor_pressure_" + sTimeStamp + ".tif"

def windSpeed(inDateTime):

    #Caclulate average parameter values (wind speed, wind direction, air temperature) over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "wind_Table", "air_temperature > -500 AND wind_speed_average > -500 AND wind_direction > -500")
    arcpy.Statistics_analysis("wind_Table", "in_memory/wind_Table2", [["air_temperature", "MEAN"], ["wind_speed_average", "MEAN"], ["wind_direction", "MEAN"]], "site_key")
    #Copy parameter values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/wind_Table2", "site_key", ["MEAN_air_temperature", "MEAN_wind_speed_average", "MEAN_wind_direction"])
    #Delete rows from station feature class that have null values for air temperature, wind speed or wind direction
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None or row.getValue("MEAN_wind_speed_average") is None or row.getValue("MEAN_wind_direction") is None or row.getValue("Site_Key") == 'rmsp3':
            cursor.deleteRow(row)
    del cursor
    del row



    #Set up time parameters
    #ninjaPath = "C:/WindNinja/WindNinja-2.5.1/bin/WindNinja_cli.exe"
    sWindDate = inDateTime.split(" ")[0]
    sWindTime = inDateTime.split(" ")[1]
    lsWindDate = sWindDate.split("-")
    lsWindTime = sWindTime.split(":")
    sWindYear = lsWindDate[0]
    sWindMonth = lsWindDate[1]
    sWindDay = lsWindDate[2]
    sWindHour = lsWindTime[0]
    sWindMinute = lsWindTime[1]

    #Build station csv file from SQL data
    #open(scratchWS + "/wnStations.csv", 'w').close()
    #Add coordinates to station feature class
    arcpy.AddGeometryAttributes_management(scratchGDB + "/tempStations", "POINT_X_Y_Z_M")
    #Loop through stations in station feature class and write parameter values to a csv file
    with open(scratchWS + "/wnStations.csv", 'wb') as csvFile:
        a = csv.writer(csvFile)
        #write header row
        a.writerow(['Station_Name', 'Coord_Sys(PROJCS,GEOGCS)', 'Datum(WGS84,NAD83,NAD27)', 'Lat/YCoord', 'Lon/XCoord', 'Height', 'Height_Units(meters,feet)', 'Speed', 'Speed_Units(mph,kph,mps)', 'Direction(degrees)', 'Temperature', 'Temperature_Units(F,C)', 'Cloud_Cover(%)', 'Radius_of_Influence', 'Radius_of_Influence_Units(miles,feet,meters,km)'])
        srchCursor = arcpy.SearchCursor(scratchGDB + "/tempStations")
        for row in srchCursor:
            a.writerow([row.getValue("Site_Key"), 'PROJCS', 'NAD83', row.getValue("Point_Y"), row.getValue("Point_X"), '3', 'meters', row.getValue("MEAN_wind_speed_average"), 'mps', row.getValue("MEAN_wind_direction"), row.getValue("MEAN_air_temperature"), 'C', '0', '-1', 'miles'])

    #List arguments for WindNinja CLI
    args = []
##    args = [ninjaPath,
##    "--initialization_method", "pointInitialization",
##    "--elevation_file", r'C:\AA_Thesis_Project\ZZ_MySQL_Work\rc_elevation_filled.tif', #elevation raster (cannot contain any "no-data" values)
##    "--match_points", "false", #match simulations to points (simulation fails if set to true)
##    "--year", sWindYear,
##    "--month", sWindMonth,
##    "--day", sWindDay,
##    "--hour", sWindHour,
##    "--minute", sWindMinute,
##    "--mesh_resolution", output_cell_size, #Resolution of model calculations
##    "--vegetation", "brush", #Vegetation type (can be 'grass', 'brush', or 'trees')
##    "--time_zone", "America/Boise", #time zone of target simulation
##    "--diurnal_winds", "true", #consider diurnal cycles in calculations
##    "--write_goog_output", "false", #write kml output (boolean: true/false)
##    "--write_shapefile_output", "false", #write shapefile output (boolean: true/false)
##    "--write_farsite_atm", "false", #write fire behavior file (boolean: true/false)
##    "--write_ascii_output", "true", #write ascii file output (this should always be set to true)
##    "--ascii_out_resolution", "-1", #resolution of output (-1 means same as mesh_resolution)
##    "--units_ascii_out_resolution", "m",
##    "--units_mesh_resolution", "m", #units of resolution of model calculations (should be "m" for meters)
##    "--units_output_wind_height", "m", #units of output wind height
##    "--output_speed_units", "mps",
##    "--output_wind_height", "3",
##    "--wx_station_filename", scratchWS + "/wnStations.csv", #weather station csv file used in point initialization method
##    "--output_path", scratchWS] #path to output

    #run the WindNinja_cli.exe (output is written to same location as elevation raster)
    arcpy.AddMessage("Calling WindNinja command line interface")
    runfile = subprocess.Popen(args, stdout = subprocess.PIPE, bufsize = -1)
    runfile.wait()
    output = runfile.stdout.read()
    if output is None:
        arcpy.AddMessage("Results: None returned\n")
    else:
        arcpy.AddMessage("Results:\n" + output)

    #convert ascii file to new grid
    for file in os.listdir(scratchWS):
        if file.endswith("_vel.asc"):
            path2ASCII = scratchWS + "/" + file
    arcpy.ASCIIToRaster_conversion(in_ascii_file=path2ASCII, out_raster=outFolder + "/wind_speed_" + sTimeStamp + ".tif", data_type="FLOAT")

    #Get coordinate system information
    desc = arcpy.Describe(rc_elevation)
    coordSystem = desc.spatialReference
    arcpy.DefineProjection_management(outFolder + "/wind_speed_" + sTimeStamp + ".tif", coordSystem)

    return outFolder + "/wind_speed_" + sTimeStamp + ".tif"

def deleteScratch(inList):
    for path in inList:
        arcpy.Delete_management(path)

#START MAIN SCRIPT
#Calculate time range and number of time steps
dateFrom = roundTime(datetime.datetime.strptime(sFromDate, "%Y-%m-%d %H:%M:%S"), 60*60)
dateTo = roundTime(datetime.datetime.strptime(sToDate, "%Y-%m-%d %H:%M:%S"))
timeDelta = dateTo - dateFrom
iDays = timeDelta.days
iSeconds = timeDelta.seconds + (iDays * 86400)
iMinutes = iSeconds / 60
iHours = iMinutes / 60
iNumSteps = iHours / iTimeStep

#connect to MySQL database
try:
  cnx = mysql.connector.connect(user='root', password='',
                                host='localhost',
                                database='jd_data',
                                buffered=True)
except mysql.connector.Error as err:

    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        arcpy.AddMessage("Something is wrong with your user name or password")
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        arcpy.AddMessage("Database does not exist")
    else:
        arcpy.AddMessage(err)
else:
    arcpy.AddMessage("Connection successful")

#Initiate list of scratch data to delete at end of script
lsScratchData = []

#Initiate list of output .tif files
lsOutput = []

#Create master station feature class that will be copied for each gridding function
fcStations_wElevation = scratchGDB + "/stations_wElevation"
lsScratchData.append(fcStations_wElevation)
##arcpy.WFSToFeatureClass_conversion(input_WFS_server="https://miles-gis-sandbox.northwestknowledge.net/arcgis/services/WCWAVE/rc_data_mapService/MapServer/WFSServer", \
##    WFS_feature_type="station_locations", out_path=scratchGDB, out_name="station_locations")
station_locations = scratchGDB + "/station_locations"
lsScratchData.append(station_locations)
arcpy.CopyFeatures_management(r'C:\Users\johawesl\Desktop\JD_Work\Required_Data.gdb\station_locations_JD2', station_locations)
extFullFeatures = arcpy.Describe(scratchGDB + "/station_locations").extent
#station_locations = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\station_locations'

rc_elevation = r'C:\Users\johawesl\Desktop\JD_Work\Required_Data.gdb\JD_Clipped_10m_2'
arcpy.env.cellSize = rc_elevation
output_cell_size = arcpy.env.cellSize
extElevation = arcpy.Describe(rc_elevation).extent
arcpy.env.extent = extElevation
arcpy.gp.ExtractValuesToPoints_sa(station_locations, rc_elevation, fcStations_wElevation, "NONE", "VALUE_ONLY")

#Run separate gridding functions
delta = datetime.timedelta(hours=iTimeStep)
dateIncrement = dateFrom
while dateIncrement < dateTo:
    #Check which tables in the SQL database will have to be queried
    if any([boolAllTools, boolAirTemperature, boolDewPoint, boolVaporPressure, boolWindSpeed, boolSolarRadiation, boolThermalRadiation]):
        lsScratchData_Imd = []

        #Initiate parameter lists
        lsSiteKey = []
        lsDatetime = []
        lsAirTemperature = []
        lsVaporPressure = []
        lsDewpoint = []
        lsSolarRadiation = []
        lsWindSpeed = []
        lsWindDirection = []

        #Query climate table
        sFrom = dateIncrement.strftime("%Y-%m-%d %H:%M:%S")
        sTimeStamp = dateIncrement.strftime("%Y%m%d_%H")
        dateToTemp = dateIncrement + delta
        sTo = dateToTemp.strftime("%Y-%m-%d %H:%M:%S")
        sQuery = "SELECT * FROM weather WHERE date_time >= '" + sFrom + "' AND date_time < '" + sTo + "'"
        cur = cnx.cursor()
        cur.execute(sQuery)
        iNumReturn = cur.rowcount
        arcpy.AddMessage("sQuery: " + sQuery)
        arcpy.AddMessage("iNumReturn: " + str(iNumReturn))

        #Append query results to parameter lists
        for i in range(0,iNumReturn):
            row = cur.fetchone()
            lsSiteKey.append(row[0])
            lsDatetime.append(row[1])
            lsAirTemperature.append(row[8])
            lsVaporPressure.append(row[10])
            lsDewpoint.append(row[11])
            lsSolarRadiation.append(row[12])
            lsWindSpeed.append(row[13])
            lsWindDirection.append(row[14])
        cur.close()

        #Build climate table to pass to gridding tools
        arcpy.CreateTable_management(scratchGDB,"climate_table1")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "date_time", "DATE")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "air_temperature", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "vapor_pressure", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "dewpoint_temperature", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "in_solar_radiation", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "wind_speed_average", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/climate_table1", "wind_direction", "FLOAT")
        lsScratchData_Imd.append(scratchGDB + "/climate_table1")

        #Add values from parameter lists to climate table
        inCursor = arcpy.InsertCursor(scratchGDB + "/climate_table1")
        for j in range(0, iNumReturn):
            row = inCursor.newRow()
            row.setValue("site_key", lsSiteKey[j])
            row.setValue("date_time", lsDatetime[j])
            row.setValue("air_temperature", lsAirTemperature[j])
            row.setValue("vapor_pressure", lsVaporPressure[j])
            row.setValue("dewpoint_temperature", lsDewpoint[j])
            row.setValue("in_solar_radiation", lsSolarRadiation[j])
            row.setValue("wind_speed_average", lsWindSpeed[j])
            row.setValue("wind_direction", lsWindDirection[j])
            inCursor.insertRow(row)
        del inCursor
        del row

        #Run gridding functions
        if boolAirTemperature:
            pathAirTemperature = airTemperature()
            lsOutput.append(pathAirTemperature)
        if boolDewPoint:
            pathDewPointTemperature, pathPercentSnow, pathPrecipSnowDensity = dewPoint()
            lsOutput.append(pathDewPointTemperature)
            lsOutput.append(pathPercentSnow)
            lsOutput.append(pathPrecipSnowDensity)
        if boolVaporPressure:
            pathVaporPressure = vaporPressure()
            lsOutput.append(pathVaporPressure)
        if boolWindSpeed:
            pathWindSpeed = windSpeed(sFrom)
            lsOutput.append(pathWindSpeed)
        if boolSolarRadiation:
            pathSolarRadiation = solarRadiation()
            lsOutput.append(pathSolarRadiation)
        if boolThermalRadiation:
            #Query database for average air temperature for current day
            sFromTR = dateIncrement.strftime("%Y-%m-%d")
            sQuery2 = "SELECT AVG(NULLIF(" + tmp3 + ", -999)) FROM weather WHERE date_time >= '" + sFromTR + " 00:00:00" + "' AND date_time <= '" + sFromTR + " 23:00:00'"
            cur2 = cnx.cursor()
            cur2.execute(sQuery2)
            dRefTemp = cur2.fetchone()[0]
            cur2.close()

            pathThermalRadiation = thermalRadiation(pathAirTemperature, pathVaporPressure, dRefTemp)

            lsOutput.append(pathThermalRadiation)


        #Delete intermediate scratch
        deleteScratch(lsScratchData_Imd)

    if any([boolAllTools, boolPrecipMass]):
        lsScratchData_Imd = []

        #Initiate parameter lists
        lsSiteKey = []
        lsPPTS = []
        lsPPTU = []
        lsPPTA = []

        #Query precipitation table
        sFrom = dateIncrement.strftime("%Y-%m-%d %H:%M:%S")
        sTimeStamp = dateIncrement.strftime("%Y%m%d_%H")
        dateToTemp = dateIncrement + delta
        sTo = dateToTemp.strftime("%Y-%m-%d %H:%M:%S")
        sQuery = "SELECT * FROM precipitation WHERE date_time >= '" + sFrom + "' AND date_time < '" + sTo + "'"
        cur = cnx.cursor()
        cur.execute(sQuery)
        iNumReturn = cur.rowcount

        #Append query results to parameter lists
        for i in range(0,iNumReturn):
            row = cur.fetchone()
            lsSiteKey.append(row[0])
            lsPPTS.append(row[2])
            lsPPTU.append(row[3])
            lsPPTA.append(row[4])
        cur.close()

        #Build precipitation table to pass to gridding tools
        arcpy.CreateTable_management(scratchGDB,"precipitation_table1")
        arcpy.AddField_management(scratchGDB + "/precipitation_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratchGDB + "/precipitation_table1", "ppts", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/precipitation_table1", "pptu", "FLOAT")
        arcpy.AddField_management(scratchGDB + "/precipitation_table1", "ppta", "FLOAT")
        lsScratchData_Imd.append(scratchGDB + "/precipitation_table1")

        #Add values from parameter lists to precipitation table
        inCursor = arcpy.InsertCursor(scratchGDB + "/precipitation_table1")
        for j in range(0, iNumReturn):
            row = inCursor.newRow()
            row.setValue("site_key", lsSiteKey[j])
            row.setValue("ppts", lsPPTS[j])
            row.setValue("pptu", lsPPTU[j])
            row.setValue("ppta", lsPPTA[j])
            inCursor.insertRow(row)
        del inCursor
        del row

        #Run gridding functions
        if boolPrecipMass:
            pathPrecipMass = precipMass()
            lsOutput.append(pathPrecipMass)

        #Delete intermediate scratch
        deleteScratch(lsScratchData_Imd)

    if any([boolAllTools, boolSoilTemperature]):
        lsScratchData_Imd = []

        #Initiate parameter lists
        lsSiteKey = []
        lsSoilTemp = []

        #Query soil temperature table
        sFrom = dateIncrement.strftime("%Y-%m-%d %H:%M:%S")
        sTimeStamp = dateIncrement.strftime("%Y%m%d_%H")
        dateToTemp = dateIncrement + delta
        sTo = dateToTemp.strftime("%Y-%m-%d %H:%M:%S")
        sQuery = "SELECT * FROM soil_temperature WHERE date_time >= '" + sFrom + "' AND date_time < '" + sTo + "'"
        cur = cnx.cursor()
        cur.execute(sQuery)
        iNumReturn = cur.rowcount

        #Append query results to parameter lists
        for i in range(0,iNumReturn):
            row = cur.fetchone()
            lsSiteKey.append(row[0])
            lsSoilTemp.append(row[3]) #Column 3 is temperature at 5 cm depth
        cur.close()

        #Build soil temperature table to pass to gridding tools
        arcpy.CreateTable_management(scratchGDB,"soilTemperature_table1")
        arcpy.AddField_management(scratchGDB + "/soilTemperature_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratchGDB + "/soilTemperature_table1", "stm005", "FLOAT")
        lsScratchData_Imd.append(scratchGDB + "/soilTemperature_table1")

        #Add values from parameter lists to soil temperature table
        inCursor = arcpy.InsertCursor(scratchGDB + "/soilTemperature_table1")
        for j in range(0, iNumReturn):
            row = inCursor.newRow()
            row.setValue("site_key", lsSiteKey[j])
            row.setValue("stm005", lsSoilTemp[j])
            inCursor.insertRow(row)
        del inCursor
        del row

        #Run gridding functions
        if boolSoilTemperature:
            pathSoilTemperature = soilTemperature()
            lsOutput.append(pathSoilTemperature)

        #Delete intermediate scratch
        deleteScratch(lsScratchData_Imd)

    dateIncrement += delta




#Run initial condition functions once
sFrom = dateFrom.strftime("%Y-%m-%d %H:%M:%S")
sTimeStamp = dateFrom.strftime("%Y%m%d_%H")
delta = datetime.timedelta(hours=1)
dateToTemp = dateFrom + delta
sTo = dateToTemp.strftime("%Y-%m-%d %H:%M:%S")

#Initial snow depth:
if any([boolAllTools, boolSnowDepth]):
    lsScratchData_Imd = []

    #Initiate parameter lists
    lsSiteKey = []
    lsSnowDepth = []

    #Query snow depth table
    sQuery = "SELECT * FROM snow_depth WHERE date_time >= '" + sFrom + "' AND date_time < '" + sTo + "'"
    cur = cnx.cursor()
    cur.execute(sQuery)
    iNumReturn = cur.rowcount

    #Append query results to parameter lists
    for i in range(0,iNumReturn):
        row = cur.fetchone()
        lsSiteKey.append(row[0])
        lsSnowDepth.append(row[8])
    cur.close()

    #Build snow depth table to pass to gridding tools
    arcpy.CreateTable_management(scratchGDB,"snowDepth_table1")
    arcpy.AddField_management(scratchGDB + "/snowDepth_table1", "site_key", "TEXT")
    arcpy.AddField_management(scratchGDB + "/snowDepth_table1", "snow_depth", "FLOAT")
    lsScratchData_Imd.append(scratchGDB + "/snowDepth_table1")

    #Add values from parameter lists to snow depth table
    inCursor = arcpy.InsertCursor(scratchGDB + "/snowDepth_table1")
    for j in range(0, iNumReturn):
        row = inCursor.newRow()
        row.setValue("site_key", lsSiteKey[j])
        row.setValue("snow_depth", lsSnowDepth[j])
        inCursor.insertRow(row)
    del inCursor
    del row

    #Run gridding functions
    if boolSnowDepth:
        pathSnowDepth = snowDepth()
        lsOutput.append(pathSnowDepth)

    #Delete intermediate scratch
    deleteScratch(lsScratchData_Imd)

#Snow properties:
if boolSnowProperties:
    pathULSnowTemperature, pathAverageSnowTemperature, pathSnowDensity = snowProperties()
    lsOutput.append(pathULSnowTemperature)
    lsOutput.append(pathAverageSnowTemperature)
    lsOutput.append(pathSnowDensity)

#Constants
if boolConstants:
    pathRLConstant, pathH2OConstant = constants()
    lsOutput.append(pathRLConstant)
    lsOutput.append(pathH2OConstant)


#Close connection to database
cnx.close()

#Delete scratch data
deleteScratch(lsScratchData)

#Zip output folder and return as output zip file
shutil.make_archive(outFolder,'zip',outFolder)




#Set output parameter as list of output
arcpy.SetParameterAsText(21, outFolder + ".zip")


