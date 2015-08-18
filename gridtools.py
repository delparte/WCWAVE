'''
Copyright (c) 2015 Wesley Joel Johansen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

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
s_from_date = arcpy.GetParameterAsText(0)
s_to_date = arcpy.GetParameterAsText(1)
i_time_step = int(arcpy.GetParameterAsText(2))
s_krig_method = arcpy.GetParameterAsText(3)
bool_all_tools = arcpy.GetParameter(4)
bool_air_temperature = arcpy.GetParameter(5)

bool_constants = arcpy.GetParameter(6)
d_rl_constant = arcpy.GetParameter(7)
d_h20_constant = arcpy.GetParameter(8)

bool_dew_point = arcpy.GetParameter(9)
bool_precip_mass = arcpy.GetParameter(10)
bool_snow_depth = arcpy.GetParameter(11)

bool_snow_properties = arcpy.GetParameter(12)
ll_interp_values = json.loads(arcpy.GetParameter(13).JSON)
ul_interp_values = json.loads(arcpy.GetParameter(14).JSON)
density_interp_values = json.loads(arcpy.GetParameter(15).JSON)

bool_soil_temperature = arcpy.GetParameter(16)
bool_solar_radiation = arcpy.GetParameter(17)
bool_thermal_radiation = arcpy.GetParameter(18)
bool_VaporPressure = arcpy.GetParameter(19)
bool_wind_speed = arcpy.GetParameter(20)
#tempOutput = arcpy.GetParameterAsText(16)

#Specify workspace
scratch_ws = arcpy.env.scratchFolder
arcpy.AddMessage("Scratch Workspace: " + scratch_ws)
#scratch_ws = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\scratch'
#Make output folder to zip
date_now = datetime.datetime.now()
s_now = date_now.strftime("%Y%d%b_%H%M")
os.makedirs(scratch_ws + "/Output_" + s_now)
out_folder = scratch_ws + "/Output_" + s_now
arcpy.AddMessage("Output Folder: " + out_folder)
scratch_gdb = arcpy.env.scratchGDB
arcpy.env.overwriteOutput = True

#Define Functions
def roundTime(dt, roundTo=60):
    seconds = (dt - dt.min).seconds
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

def airTemperature():
    #Caclulate average air temperatures over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "airTemperature_Table", "air_temperature > -500")
    arcpy.Statistics_analysis("airTemperature_Table", "in_memory/airTemperature_Table2", "air_temperature MEAN", "site_key")
    #Copy air temperature values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/airTemperature_Table2", "site_key", "MEAN_air_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for air temperature
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if s_krig_method == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratch_gdb +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratch_gdb +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratch_gdb, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_air_temperature","RASTERVALU", scratch_gdb + "/coef_table","","")
        ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratch_gdb + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_air_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        arcpy.gp.Kriging_sa(scratch_gdb + "/tempStations", "residual", scratch_gdb + "/airTemperature_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        ls_scratch_data_imd.append(scratch_gdb + "/airTemperature_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratch_gdb + "/airTemperature_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(out_folder + "/air_temperature_" + s_time_stamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratch_gdb + "/tempStations", z_field="MEAN_air_temperature", out_ga_layer="#", \
            out_raster=out_folder + "/air_temperature_" + s_time_stamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="50", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

    return out_folder + "/air_temperature_" + s_time_stamp + ".tif"

def constants():
    #Get coordinate system information
    desc = arcpy.Describe(rc_elevation)
    coord_system = desc.spatialReference
    rl_constant = CreateConstantRaster(d_rl_constant, "FLOAT", output_cell_size)
    arcpy.DefineProjection_management(rlConstant, coord_system)
    rl_constant.save(out_folder + "/roughness_length_" + s_time_stamp + ".tif")
    water_constant = CreateConstantRaster(d_h20_constant, "FLOAT", output_cell_size)
    arcpy.DefineProjection_management(water_constant, coord_system)
    water_constant.save(out_folder + "/H2O_saturation_" + s_time_stamp + ".tif")

    return out_folder + "/roughness_length_" + s_time_stamp + ".tif", out_folder + "/H2O_saturation_" + s_time_stamp + ".tif"

def dewPoint():
    #Caclulate average dew point temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "dewPoint_Table", "dewpoint_temperature > -500")
    arcpy.Statistics_analysis("dewPoint_Table", "in_memory/dewPoint_Table2", "dewpoint_temperature MEAN", "site_key")
    #Copy dew point temperature values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/dewPoint_Table2", "site_key", "MEAN_dewpoint_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for dew point temperatures
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_dewpoint_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if s_krig_method == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratch_gdb +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratch_gdb +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratch_gdb, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_dewpoint_temperature","RASTERVALU", scratch_gdb + "/coef_table","","")
        ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratch_gdb + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_dewpoint_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratch_gdb + "/tempStations", "residual", scratch_gdb + "/dewPoint_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratch_gdb + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratch_gdb + "/dewPoint_residual")
        ls_scratch_data_imd.append(scratch_gdb + "/dewPoint_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratch_gdb + "/dewPoint_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(out_folder + "/dew_point_temperature_" + s_time_stamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratch_gdb + "/tempStations", z_field="MEAN_dewpoint_temperature", out_ga_layer="#", \
            out_raster=out_folder + "/dew_point_temperature_" + s_time_stamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

    #Create precipitation properties rasters (% precip as snow, density of snow portion of precip)
    #Percent snow conditional:
    inRas = Raster(out_folder + "/dew_point_temperature_" + s_time_stamp + ".tif")
    outPercentSnowRaster = arcpy.sa.Con(inRas < -5.0, 1.0, Con((inRas >= -5.0) & (inRas < -3.0), 1.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 1.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 1.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0.0), 0.75,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 0.25,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outPercentSnowRaster, out_folder + "/percent_snow_" + s_time_stamp + ".tif",pixel_type="32_BIT_FLOAT")

    #Snow density conditional:
    outSnowDensityRaster = arcpy.sa.Con(inRas < -5.0, 75.0, Con((inRas >= -5.0) & (inRas < -3.0), 100.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 150.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 175.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0), 200.0,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 250.0,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outSnowDensityRaster,out_folder + "/precipitation_snow_density_" + s_time_stamp + ".tif",pixel_type="32_BIT_FLOAT")


    return out_folder + "/dew_point_temperature_" + s_time_stamp + ".tif", out_folder + "/percent_snow_" + s_time_stamp + ".tif", out_folder + "/precipitation_snow_density_" + s_time_stamp + ".tif"

def precipMass():
    #Caclulate average precipitation values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/precipitation_table1", "precipitation_Table", "ppts > -500")
    arcpy.Statistics_analysis("precipitation_Table", "in_memory/precipitation_Table2", "ppts MEAN", "site_key")
    #Copy precipitation values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/precipitation_Table2", "site_key", "MEAN_ppts")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for precipitation
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_ppts") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if s_krig_method == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratch_gdb +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratch_gdb +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratch_gdb, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_ppts","RASTERVALU", scratch_gdb + "/coef_table","","")
        ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratch_gdb + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_ppts") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratch_gdb + "/tempStations", "residual", scratch_gdb + "/precipitation_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratch_gdb + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratch_gdb + "/precipitation_residual")
        ls_scratch_data_imd.append(scratch_gdb + "/precipitation_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratch_gdb + "/precipitation_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(out_folder + "/precipitation_mass_" + s_time_stamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratch_gdb + "/tempStations", z_field="MEAN_ppts", out_ga_layer="#", \
            out_raster=out_folder + "/precipitation_mass_" + s_time_stamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

    return out_folder + "/precipitation_mass_" + s_time_stamp + ".tif"

def snowDepth():
    #Caclulate average snow depth values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/snowDepth_table1", "snowDepth_Table", "snow_depth > -500")
    arcpy.Statistics_analysis("snowDepth_Table", "in_memory/snowDepth_Table2", "snow_depth MEAN", "site_key")
    #Copy snow depth values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/snowDepth_Table2", "site_key", "MEAN_snow_depth")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for snow depth
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_snow_depth") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if s_krig_method == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratch_gdb +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratch_gdb +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratch_gdb, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_snow_depth","RASTERVALU", scratch_gdb + "/coef_table","","")
        ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratch_gdb + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_snow_depth") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratch_gdb + "/tempStations", "residual", scratch_gdb + "/snowDepth_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratch_gdb + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratch_gdb + "/snowDepth_residual")
        ls_scratch_data_imd.append(scratch_gdb + "/snowDepth_residual")
        #Add back elevation trends and save final raster
        imdRaster = arcpy.Raster(scratch_gdb + "/snowDepth_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster = Con(imdRaster > 0, imdRaster, 0)
        output_raster.save(out_folder + "/snow_depth_" + s_time_stamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratch_gdb + "/tempStations", z_field="MEAN_snow_depth", out_ga_layer="#", \
            out_raster=scratch_gdb + "/imdRaster", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")
        output_raster = Con(Raster(scratch_gdb + "/imdRaster") > 0, Raster(scratch_gdb + "/imdRaster"), 0)
        output_raster.save(out_folder + "/snow_depth_" + s_time_stamp + ".tif")

    return out_folder + "/snow_depth_" + s_time_stamp + ".tif"

def snowProperties():
    if len(density_interp_values['features']) <= 1:
        #Density Equation: y = -0.0395(elevation) + 405.26
        snow_density_raster = -0.0395 * Raster(rc_elevation) + 405.26
        snow_density_raster.save(out_folder + "/snow_density_" + s_time_stamp + ".tif")
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
        snow_density_raster.save(out_folder + "/snow_density_" + s_time_stamp + ".tif")

    if len(ul_interp_values['features']) <= 1:
        #Upper Layer Temperature Equation: y = -0.0008(elevation) + 0.1053
        upper_layer_temperature = -0.0008 * Raster(rc_elevation) + 0.1053
        upper_layer_temperature.save(out_folder + "/active_snow_layer_temperature_" + s_time_stamp + ".tif")
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
        upper_layer_temperature.save(out_folder + "/active_snow_layer_temperature_" + s_time_stamp + ".tif")

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
    average_snowcover_temperature.save(out_folder + "/average_snow_cover_temperature_" + s_time_stamp + ".tif")

    return out_folder + "/active_snow_layer_temperature_" + s_time_stamp + ".tif", out_folder + "/average_snow_cover_temperature_" + s_time_stamp + ".tif", out_folder + "/snow_density_" + s_time_stamp + ".tif"

def soilTemperature():
    #Set extent to full feature class (all stations)
    arcpy.env.extent = ext_full_features
    #Caclulate average soil temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/soilTemperature_table1", "soilTemperature_Table", "stm010 > -500")
    arcpy.Statistics_analysis("soilTemperature_Table", "in_memory/soilTemperature_Table2", "stm010 MEAN", "site_key")
    #Copy soil temperature values to station feature class
    fcStations4Soil = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\station_locations_soil'
    arcpy.CopyFeatures_management(fcStations4Soil, scratch_gdb +"/tempStations4soil")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations4soil")
    arcpy.JoinField_management(scratch_gdb +"/tempStations4soil", "Site_Key", "in_memory/soilTemperature_Table2", "site_key", "MEAN_stm010")
    #Delete rows from station feature class that have null values for soil temperature
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations4soil")
    for row in cursor:
        if row.getValue("MEAN_stm010") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratch_gdb +"/tempStations4soil", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratch_gdb +"/tempStations4soil","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratch_gdb, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations4soil","Unique_ID","in_memory/fcResid","MEAN_stm010","Elevation", scratch_gdb + "/coef_table","","")
    ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
    #Create final raster
    arcpy.env.extent = ext_elevation
    output_raster = (Raster(rc_elevation) * slope + intercept)
    output_raster.save(out_folder + "/soil_temperature_" + s_time_stamp + ".tif")



    return out_folder + "/soil_temperature_" + s_time_stamp + ".tif"


def solarRadiation():
    #Caclulate average solar radiation values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "solarRadiation_Table", "in_solar_radiation > -500")
    arcpy.Statistics_analysis("solarRadiation_Table", "in_memory/solarRadiation_Table2", "in_solar_radiation MEAN", "site_key")
    #Copy solar radiation values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/solarRadiation_Table2", "site_key", "MEAN_in_solar_radiation")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for solar radiation
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_in_solar_radiation") is None:
            cursor.deleteRow(row)
    del cursor
    del row


    #set up area solar radiation tool parameters and run the tool
    #Set up time parameters
    day_of_year = date_increment.timetuple().tm_yday
    iSRStart = int(date_increment.strftime("%H"))
    iSREnd = iSRStart + i_time_step
    inTWD = TimeWithinDay(day_of_year,iSRStart,iSREnd)

    global_radiation_raster = scratch_gdb + "/global_radiation_raster"
    skySize = 200
    outGlobalRadiation = AreaSolarRadiation(rc_elevation, "", skySize, inTWD)
    #outGlobalRadiation = AreaSolarRadiation(r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\RC_DEM_100m_South', "", skySize, inTWD)
    outGlobalRadiation = outGlobalRadiation / i_time_step

    #Correct global radiation raster for cloud conditions
    #Extract simulated global radiation values to station location feature class
    arcpy.AlterField_management(scratch_gdb + "/tempStations","RASTERVALU","Elevation")
    arcpy.gp.ExtractValues_toPoints_sa(scratch_gdb + "/tempStations", outGlobalRadiation, scratch_gdb + "/simPoints", "NONE", "VALUE_ONLY")
    ls_scratch_data_imd.append(scratch_gdb + "/simPoints")
    #add "ratio" field to station feature class
    arcpy.AddField_management(scratch_gdb + "/simPoints","ratio","FLOAT","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    #calculate "ratio" field (observed radiation / simulated radiation)
    arcpy.CalculateField_management(scratch_gdb + "/simPoints","ratio","!MEAN_in_solar_radiation!/ !RASTERVALU!","PYTHON_9.3","#")

    #convert "ratio" field to numpy array
    na = arcpy.da.TableToNumPyArray(scratch_gdb + "/simPoints", "ratio")

    #calculate average ratio
    dMeanRatio = numpy.mean(na["ratio"])
    dMeanRatio2 = numpy.asscalar(dMeanRatio)

    #multiply simulated raster by average ratio
    outGlobalRadiation_Corrected = outGlobalRadiation * dMeanRatio2
    outGlobalRadiation_Corrected.save(out_folder + "/solar_radiation_" + s_time_stamp + ".tif")

    return out_folder + "/solar_radiation_" + s_time_stamp + ".tif"


def thermalRadiation(in_air_temperature, in_vapor_pressure, in_surface_temperature):
    #Constants and re-defined variables (See Marks and Dozier (1979), pg. 160)
    z = rc_elevation
    vf = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\RC_ViewFactor_10m_South'
    T_a = in_air_temperature
    vp = in_vapor_pressure

    #Select one station and its data to use as reference for elevation, pressure, and temperature
    #Caclulate average air temperature and vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "thermalRadiation_Table", "air_temperature > -500 AND vapor_pressure > -500")
    arcpy.Statistics_analysis("thermalRadiation_Table", "in_memory/thermalRadiation_Table2", [["air_temperature", "MEAN"], ["vapor_pressure", "MEAN"]], "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/thermalRadiation_Table2", "site_key", ["MEAN_air_temperature", "MEAN_vapor_pressure"])
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure or air temperature
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None or row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Use update cursor to extract air temperature, vapor pressure, and elevation values for the first station, then delete that station
    P_m = 0.0
    T_m = 0.0
    z_m = 0.0
    T_s = in_surface_temperature
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
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
    output_thermal_radiation.save(out_folder + "/thermal_radiation_" + s_time_stamp + ".tif")

    return out_folder + "/thermal_radiation_" + s_time_stamp + ".tif"

def vaporPressure():
    #Caclulate average vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "vaporPressure_Table", "vapor_pressure > -500")
    arcpy.Statistics_analysis("vaporPressure_Table", "in_memory/vaporPressure_Table2", "vapor_pressure MEAN", "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/vaporPressure_Table2", "site_key", "MEAN_vapor_pressure")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Check if interpolation method is Empirical Bayesian or Detrended
    if s_krig_method == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratch_gdb +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratch_gdb +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratch_gdb, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratch_gdb + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_vapor_pressure","RASTERVALU", scratch_gdb + "/coef_table","","")
        ls_scratch_data_imd.append(scratch_gdb + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratch_gdb + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratch_gdb + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_vapor_pressure") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratch_gdb + "/tempStations", "residual", scratch_gdb + "/vaporPressure_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratch_gdb + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratch_gdb + "/vaporPressure_residual")
        ls_scratch_data_imd.append(scratch_gdb + "/vaporPressure_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratch_gdb + "/vaporPressure_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(out_folder + "/vapor_pressure_" + s_time_stamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratch_gdb + "/tempStations", z_field="MEAN_vapor_pressure", out_ga_layer="#", \
            out_raster=out_folder + "/vapor_pressure_" + s_time_stamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

    return out_folder + "/vapor_pressure_" + s_time_stamp + ".tif"

def windSpeed(inDateTime):

    #Caclulate average parameter values (wind speed, wind direction, air temperature) over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratch_gdb + "/climate_table1", "wind_Table", "air_temperature > -500 AND wind_speed_average > -500 AND wind_direction > -500")
    arcpy.Statistics_analysis("wind_Table", "in_memory/wind_Table2", [["air_temperature", "MEAN"], ["wind_speed_average", "MEAN"], ["wind_direction", "MEAN"]], "site_key")
    #Copy parameter values to station feature class
    arcpy.CopyFeatures_management(fc_stations_welevation, scratch_gdb +"/tempStations")
    ls_scratch_data_imd.append(scratch_gdb +"/tempStations")
    arcpy.JoinField_management(scratch_gdb +"/tempStations", "Site_Key", "in_memory/wind_Table2", "site_key", ["MEAN_air_temperature", "MEAN_wind_speed_average", "MEAN_wind_direction"])
    #Delete rows from station feature class that have null values for air temperature, wind speed or wind direction
    cursor = arcpy.UpdateCursor(scratch_gdb + "/tempStations")
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
    #open(scratch_ws + "/wnStations.csv", 'w').close()
    #Add coordinates to station feature class
    arcpy.AddGeometryAttributes_management(scratch_gdb + "/tempStations", "POINT_X_Y_Z_M")
    #Loop through stations in station feature class and write parameter values to a csv file
    with open(scratch_ws + "/wnStations.csv", 'wb') as csvFile:
        a = csv.writer(csvFile)
        #write header row
        a.writerow(['Station_Name', 'Coord_Sys(PROJCS,GEOGCS)', 'Datum(WGS84,NAD83,NAD27)', 'Lat/YCoord', 'Lon/XCoord', 'Height', 'Height_Units(meters,feet)', 'Speed', 'Speed_Units(mph,kph,mps)', 'Direction(degrees)', 'Temperature', 'Temperature_Units(F,C)', 'Cloud_Cover(%)', 'Radius_of_Influence', 'Radius_of_Influence_Units(miles,feet,meters,km)'])
        srchCursor = arcpy.SearchCursor(scratch_gdb + "/tempStations")
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
##    "--wx_station_filename", scratch_ws + "/wnStations.csv", #weather station csv file used in point initialization method
##    "--output_path", scratch_ws] #path to output

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
    for file in os.listdir(scratch_ws):
        if file.endswith("_vel.asc"):
            path2ASCII = scratch_ws + "/" + file
    arcpy.ASCIIToRaster_conversion(in_ascii_file=path2ASCII, out_raster=out_folder + "/wind_speed_" + s_time_stamp + ".tif", data_type="FLOAT")

    #Get coordinate system information
    desc = arcpy.Describe(rc_elevation)
    coord_system = desc.spatialReference
    arcpy.DefineProjection_management(out_folder + "/wind_speed_" + s_time_stamp + ".tif", coord_system)

    return out_folder + "/wind_speed_" + s_time_stamp + ".tif"

def deleteScratch(in_list):
    for path in in_list:
        arcpy.Delete_management(path)

#START MAIN SCRIPT
#Calculate time range and number of time steps
date_from = roundTime(datetime.datetime.strptime(s_from_date, "%Y-%m-%d %H:%M:%S"), 60*60)
date_to = roundTime(datetime.datetime.strptime(s_to_date, "%Y-%m-%d %H:%M:%S"))
time_delta = date_to - date_from
i_days = time_delta.days
i_seconds = time_delta.seconds + (i_days * 86400)
i_minutes = i_seconds / 60
i_hours = i_minutes / 60
i_num_steps = i_hours / i_time_step

#connect to MySQL database
try:
  cnx = mysql.connector.connect(user='root', password='',
                                host='localhost',
                                database='rc_data',
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
ls_scratch_data = []

#Initiate list of output .tif files
ls_output = []

#Create master station feature class that will be copied for each gridding function
fc_stations_welevation = scratch_gdb + "/stations_wElevation"
ls_scratch_data.append(fc_stations_welevation)
arcpy.WFSToFeatureClass_conversion(input_WFS_server="https://miles-gis-sandbox.northwestknowledge.net/arcgis/services/WCWAVE/rc_data_mapService/MapServer/WFSServer", \
    WFS_feature_type="station_locations", out_path=scratch_gdb, out_name="station_locations")
station_locations = scratch_gdb + "/station_locations"
ls_scratch_data.append(station_locations)
ext_full_features = arcpy.Describe(scratch_gdb + "/station_locations").extent
#station_locations = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\station_locations'

rc_elevation = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\Required_Data.gdb\RC_DEM_10m_South'
arcpy.env.cellSize = rc_elevation
output_cell_size = arcpy.env.cellSize
ext_elevation = arcpy.Describe(rc_elevation).extent
arcpy.env.extent = ext_elevation
arcpy.gp.ExtractValues_toPoints_sa(station_locations, rc_elevation, fc_stations_welevation, "NONE", "VALUE_ONLY")

#Run separate gridding functions
delta = datetime.timedelta(hours=i_time_step)
date_increment = date_from
while date_increment < date_to:
    #Check which tables in the SQL database will have to be queried
    if any([bool_all_tools, bool_air_temperature, bool_dew_point, bool_VaporPressure, bool_wind_speed, bool_solar_radiation, bool_thermal_radiation]):
        ls_scratch_data_imd = []

        #Initiate parameter lists
        ls_site_key = []
        ls_date_time = []
        ls_air_temperature = []
        ls_vapor_pressure = []
        ls_dew_point = []
        ls_solar_radiation = []
        ls_wind_speed = []
        ls_wind_direction = []

        #Query climate table
        s_from = date_increment.strftime("%Y-%m-%d %H:%M:%S")
        s_time_stamp = date_increment.strftime("%Y%m%d_%H")
        date_to_temp = date_increment + delta
        s_to = date_to_temp.strftime("%Y-%m-%d %H:%M:%S")
        s_query = "SELECT * FROM climate WHERE date_time >= '" + s_from + "' AND date_time < '" + s_to + "'"
        cur = cnx.cursor()
        cur.execute(s_query)
        i_num_return = cur.rowcount
        arcpy.AddMessage("s_query: " + s_query)
        arcpy.AddMessage("i_num_return: " + str(i_num_return))

        #Append query results to parameter lists
        for i in range(0,i_num_return):
            row = cur.fetchone()
            ls_site_key.append(row[0])
            ls_date_time.append(row[1])
            ls_air_temperature.append(row[9])
            ls_vapor_pressure.append(row[11])
            ls_dew_point.append(row[12])
            ls_solar_radiation.append(row[13])
            ls_wind_speed.append(row[14])
            ls_wind_direction.append(row[15])
        cur.close()

        #Build climate table to pass to gridding tools
        arcpy.CreateTable_management(scratch_gdb,"climate_table1")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "date_time", "DATE")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "air_temperature", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "vapor_pressure", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "dewpoint_temperature", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "in_solar_radiation", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "wind_speed_average", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/climate_table1", "wind_direction", "FLOAT")
        ls_scratch_data_imd.append(scratch_gdb + "/climate_table1")

        #Add values from parameter lists to climate table
        in_cursor = arcpy.InsertCursor(scratch_gdb + "/climate_table1")
        for j in range(0, i_num_return):
            row = in_cursor.newRow()
            row.setValue("site_key", ls_site_key[j])
            row.setValue("date_time", ls_date_time[j])
            row.setValue("air_temperature", ls_air_temperature[j])
            row.setValue("vapor_pressure", ls_vapor_pressure[j])
            row.setValue("dewpoint_temperature", ls_dew_point[j])
            row.setValue("in_solar_radiation", ls_solar_radiation[j])
            row.setValue("wind_speed_average", ls_wind_speed[j])
            row.setValue("wind_direction", ls_wind_direction[j])
            in_cursor.insertRow(row)
        del in_cursor
        del row

        #Run gridding functions
        if bool_air_temperature:
            pathAirTemperature = airTemperature()
            ls_output.append(pathAirTemperature)
        if bool_dew_point:
            pathDewPointTemperature, pathPercentSnow, pathPrecipSnowDensity = dewPoint()
            ls_output.append(pathDewPointTemperature)
            ls_output.append(pathPercentSnow)
            ls_output.append(pathPrecipSnowDensity)
        if bool_VaporPressure:
            pathVaporPressure = vaporPressure()
            ls_output.append(pathVaporPressure)
        if bool_wind_speed:
            pathWindSpeed = windSpeed(s_from)
            ls_output.append(pathWindSpeed)
        if bool_solar_radiation:
            pathSolarRadiation = solarRadiation()
            ls_output.append(pathSolarRadiation)
        if bool_thermal_radiation:
            #Query database for average air temperature for current day
            s_fromTR = date_increment.strftime("%Y-%m-%d")
            s_query2 = "SELECT AVG(NULLIF(tmp3, -999)) FROM climate WHERE date_time >= '" + s_fromTR + " 00:00:00" + "' AND date_time <= '" + s_fromTR + " 23:00:00'"
            cur2 = cnx.cursor()
            cur2.execute(s_query2)
            dRefTemp = cur2.fetchone()[0]
            cur2.close()

            pathThermalRadiation = thermalRadiation(pathAirTemperature, pathVaporPressure, dRefTemp)

            ls_output.append(pathThermalRadiation)


        #Delete intermediate scratch
        deleteScratch(ls_scratch_data_imd)

    if any([bool_all_tools, bool_precip_mass]):
        ls_scratch_data_imd = []

        #Initiate parameter lists
        ls_site_key = []
        lsPPTS = []
        lsPPTU = []
        lsPPTA = []

        #Query precipitation table
        s_from = date_increment.strftime("%Y-%m-%d %H:%M:%S")
        s_time_stamp = date_increment.strftime("%Y%m%d_%H")
        date_to_temp = date_increment + delta
        s_to = date_to_temp.strftime("%Y-%m-%d %H:%M:%S")
        s_query = "SELECT * FROM precipitation WHERE date_time >= '" + s_from + "' AND date_time < '" + s_to + "'"
        cur = cnx.cursor()
        cur.execute(s_query)
        i_num_return = cur.rowcount

        #Append query results to parameter lists
        for i in range(0,i_num_return):
            row = cur.fetchone()
            ls_site_key.append(row[0])
            lsPPTS.append(row[2])
            lsPPTU.append(row[3])
            lsPPTA.append(row[4])
        cur.close()

        #Build precipitation table to pass to gridding tools
        arcpy.CreateTable_management(scratch_gdb,"precipitation_table1")
        arcpy.AddField_management(scratch_gdb + "/precipitation_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratch_gdb + "/precipitation_table1", "ppts", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/precipitation_table1", "pptu", "FLOAT")
        arcpy.AddField_management(scratch_gdb + "/precipitation_table1", "ppta", "FLOAT")
        ls_scratch_data_imd.append(scratch_gdb + "/precipitation_table1")

        #Add values from parameter lists to precipitation table
        in_cursor = arcpy.InsertCursor(scratch_gdb + "/precipitation_table1")
        for j in range(0, i_num_return):
            row = in_cursor.newRow()
            row.setValue("site_key", ls_site_key[j])
            row.setValue("ppts", lsPPTS[j])
            row.setValue("pptu", lsPPTU[j])
            row.setValue("ppta", lsPPTA[j])
            in_cursor.insertRow(row)
        del in_cursor
        del row

        #Run gridding functions
        if bool_precip_mass:
            pathPrecipMass = precipMass()
            ls_output.append(pathPrecipMass)

        #Delete intermediate scratch
        deleteScratch(ls_scratch_data_imd)

    if any([bool_all_tools, bool_soil_temperature]):
        ls_scratch_data_imd = []

        #Initiate parameter lists
        ls_site_key = []
        lsSoilTemp = []

        #Query soil temperature table
        s_from = date_increment.strftime("%Y-%m-%d %H:%M:%S")
        s_time_stamp = date_increment.strftime("%Y%m%d_%H")
        date_to_temp = date_increment + delta
        s_to = date_to_temp.strftime("%Y-%m-%d %H:%M:%S")
        s_query = "SELECT * FROM soil_temperature WHERE date_time >= '" + s_from + "' AND date_time < '" + s_to + "'"
        cur = cnx.cursor()
        cur.execute(s_query)
        i_num_return = cur.rowcount

        #Append query results to parameter lists
        for i in range(0,i_num_return):
            row = cur.fetchone()
            ls_site_key.append(row[0])
            lsSoilTemp.append(row[4]) #Column 4 is temperature at 10 cm depth
        cur.close()

        #Build soil temperature table to pass to gridding tools
        arcpy.CreateTable_management(scratch_gdb,"soilTemperature_table1")
        arcpy.AddField_management(scratch_gdb + "/soilTemperature_table1", "site_key", "TEXT")
        arcpy.AddField_management(scratch_gdb + "/soilTemperature_table1", "stm010", "FLOAT")
        ls_scratch_data_imd.append(scratch_gdb + "/soilTemperature_table1")

        #Add values from parameter lists to soil temperature table
        in_cursor = arcpy.InsertCursor(scratch_gdb + "/soilTemperature_table1")
        for j in range(0, i_num_return):
            row = in_cursor.newRow()
            row.setValue("site_key", ls_site_key[j])
            row.setValue("stm010", lsSoilTemp[j])
            in_cursor.insertRow(row)
        del in_cursor
        del row

        #Run gridding functions
        if bool_soil_temperature:
            pathSoilTemperature = soilTemperature()
            ls_output.append(pathSoilTemperature)

        #Delete intermediate scratch
        deleteScratch(ls_scratch_data_imd)

    date_increment += delta




#Run initial condition functions once
s_from = date_from.strftime("%Y-%m-%d %H:%M:%S")
s_time_stamp = date_from.strftime("%Y%m%d_%H")
delta = datetime.timedelta(hours=1)
date_to_temp = date_from + delta
s_to = date_to_temp.strftime("%Y-%m-%d %H:%M:%S")

#Initial snow depth:
if any([bool_all_tools, bool_snow_depth]):
    ls_scratch_data_imd = []

    #Initiate parameter lists
    ls_site_key = []
    ls_snow_depth = []

    #Query snow depth table
    s_query = "SELECT * FROM snow_depth WHERE date_time >= '" + s_from + "' AND date_time < '" + s_to + "'"
    cur = cnx.cursor()
    cur.execute(s_query)
    i_num_return = cur.rowcount

    #Append query results to parameter lists
    for i in range(0,i_num_return):
        row = cur.fetchone()
        ls_site_key.append(row[0])
        ls_snow_depth.append(row[9])
    cur.close()

    #Build snow depth table to pass to gridding tools
    arcpy.CreateTable_management(scratch_gdb,"snowDepth_table1")
    arcpy.AddField_management(scratch_gdb + "/snowDepth_table1", "site_key", "TEXT")
    arcpy.AddField_management(scratch_gdb + "/snowDepth_table1", "snow_depth", "FLOAT")
    ls_scratch_data_imd.append(scratch_gdb + "/snowDepth_table1")

    #Add values from parameter lists to snow depth table
    in_cursor = arcpy.InsertCursor(scratch_gdb + "/snowDepth_table1")
    for j in range(0, i_num_return):
        row = in_cursor.newRow()
        row.setValue("site_key", ls_site_key[j])
        row.setValue("snow_depth", ls_snow_depth[j])
        in_cursor.insertRow(row)
    del in_cursor
    del row

    #Run gridding functions
    if bool_snow_depth:
        path_snow_depth = snowDepth()
        ls_output.append(path_snow_depth)

    #Delete intermediate scratch
    deleteScratch(ls_scratch_data_imd)

#Snow properties:
if bool_snow_properties:
    pathULSnowTemperature, pathAverageSnowTemperature, pathSnowDensity = snowProperties()
    ls_output.append(pathULSnowTemperature)
    ls_output.append(pathAverageSnowTemperature)
    ls_output.append(pathSnowDensity)

#Constants
if bool_constants:
    pathRLConstant, pathH2OConstant = constants()
    ls_output.append(pathRLConstant)
    ls_output.append(pathH2OConstant)


#Close connection to database
cnx.close()

#Delete scratch data
deleteScratch(ls_scratch_data)

#Zip output folder and return as output zip file
shutil.make_archive(out_folder,'zip',out_folder)




#Set output parameter as list of output
arcpy.SetParameterAsText(21, out_folder + ".zip")


