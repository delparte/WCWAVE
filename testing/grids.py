import arcpy
from arcpy.sa import *
import os
import sys
import csv
import traceback
import numpy
#from scipy import stats
import mysql.connector
from mysql.connector import errorcode
import datetime
import json
import shutil
import subprocess

data = {
    'watershed' : arcpy.GetParameterAsText(0),
    'from_date' : arcpy.GetParameterAsText(1),
    'to_date' : arcpy.GetParameterAsText(2),
    'time_step' : arcpy.GetParameterAsText(3),
    'kriging_Method' : arcpy.GetParameterAsText(4),
    'bool_all_tools' : arcpy.GetParameter(5),
    'bool_air_temp' : arcpy.GetParameter(6),
    'bool_constants' : arcpy.GetParameter(7),
    'rl_constant' : arcpy.GetParameter(8),
    'h20_constant' : arcpy.GetParameter(9),
    'bool_dew_point' : arcpy.GetParameter(10),
    'bool_precip_mass' : arcpy.GetParameter(11),
    'bool_snow_depth' : arcpy.GetParameter(12),
    'bool_snow_properties' : arcpy.GetParameter(13),
    'll_interp_values' : json.loads(arcpy.GetParameter(14).JSON),
    'ul_interp_values' : json.loads(arcpy.GetParameter(15).JSON),
    'bool_soil_temperature' : arcpy.GetParameter(17),
    'bool_solar_radiation' : arcpy.GetParameter(18),
    'bool_thermal_radiation' : arcpy.GetParameter(19),
    'bool_vapor_pressure' : arcpy.GetParameter(20),
    'bool_wind_speed' : arcpy.GetParameter(21)
}

#Set up workspaces
scratchWS = arcpy.env.scratchFolder
arcpy.env.workspace = scratchWS
arcpy.env.scratchWorkspace = scratchWS
arcpy.AddMessage("Scratch Workspace: " + scratchWS)
scratchGDB = arcpy.env.scratchGDB
arcpy.env.overwriteOutput = True

#Add to data dict
data['scratch_ws'] = scratchWS
data['scratch_gdb'] = scratchGDB

date_now = datetime.datetime.now()
s_now = date_now.strftime('%Y%d%b_%H%M%S')
os.makedirs(scratchWS + '/Output_' + s_now)
outFolder = scratchWS + '/Output_' + s_now
arcpy.AddMessage('Output Folder: ' + outFolder)

#Define Functions
def roundTime(dt, roundTo=60):
    seconds = (dt - dt.min).seconds
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)


