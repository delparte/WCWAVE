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

## data = {'watershed' : arcpy.GetParameterAsText(0),
##     'from_date' : arcpy.GetParameterAsText(1),
##     'to_date' : arcpy.GetParameterAsText(2),
##     'time_step' : int(arcpy.GetParameterAsText(3)),
##     'kriging_Method' : arcpy.GetParameterAsText(4),
##     'bool_all_tools' : arcpy.GetParameter(5),
##     'bool_air_temperature' : arcpy.GetParameter(6),
##     'bool_constants' : arcpy.GetParameter(7),
##     'rl_constant' : arcpy.GetParameter(8),
##     'h20_constant' : arcpy.GetParameter(9),
##     'bool_dew_point' : arcpy.GetParameter(10),
##     'bool_precip_mass' : arcpy.GetParameter(11),
##     'bool_snow_depth' : arcpy.GetParameter(12),
##     'bool_snow_properties' : arcpy.GetParameter(13),
##     'll_interp_values' : json.loads(arcpy.GetParameter(14).JSON),
##     'ul_interp_values' : json.loads(arcpy.GetParameter(15).JSON),
##     'bool_soil_temperature' : arcpy.GetParameter(17),
##     'bool_solar_radiation' : arcpy.GetParameter(18),
##     'bool_thermal_radiation' : arcpy.GetParameter(19),
##     'bool_vapor_pressure' : arcpy.GetParameter(20),
##     'bool_wind_speed' : arcpy.GetParameter(21)
## }

data = {'ll_interp_values': {u'fieldAliases': {u'Elevation': u'Elevation', u'Temperature': u'Temperature', u'OBJECTID': u'OBJECTID'}, u'fields': [{u'alias': u'OBJECTID', u'type': u'esriFieldTypeOID', u'name': u'OBJECTID'}, {u'alias': u'Elevation', u'type': u'esriFieldTypeSingle', u'name': u'Elevation'}, {u'alias': u'Temperature', u'type': u'esriFieldTypeSingle', u'name': u'Temperature'}], u'displayFieldName': u'', u'features': []},
    'bool_air_temperature': True, 
    'bool_vapor_pressure': False, 
    'to_date': u'2008-01-01 13:00:00', 
    'time_step': 1, 
    'bool_soil_temperature': False, 
    'rl_constant': 0.005, 
    'from_date': u'2008-01-01 12:00:00', 
    'bool_solar_radiation': False, 
    'bool_all_tools': False, 
    'h20_constant': 0.2, 
    'db': 'jd_data', 
    'bool_dew_point': False, 
    'bool_precip_mass': False, 
    'bool_wind_speed': False, 
    'kriging_Method': u'Empirical Bayesian', 
    'bool_thermal_radiation': False, 
    'bool_constants': False, 
    'bool_snow_properties': False, 
    'watershed': u'Reynolds Creek', 
    'bool_snow_depth': False, 
    'ul_interp_values': {u'fieldAliases': {u'Elevation': u'Elevation', u'Temperature': u'Temperature', u'OBJECTID': u'OBJECTID'}, u'fields': [{u'alias': u'OBJECTID', u'type': u'esriFieldTypeOID', u'name': u'OBJECTID'}, {u'alias': u'Elevation', u'type': u'esriFieldTypeSingle', u'name': u'Elevation'}, {u'alias': u'Temperature', u'type': u'esriFieldTypeSingle', u'name': u'Temperature'}], u'displayFieldName': u'', u'features': []}
    }

#Set up workspaces
scratchWS = arcpy.env.scratchFolder
arcpy.env.workspace = scratchWS
arcpy.env.scratchWorkspace = scratchWS
arcpy.AddMessage('Scratch Workspace: ' + scratchWS)
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

def selectWatershed(watershed):
    ''' Initialize all Relevant Data from the Geodatabase based on chosen watershed '''
    stations = '' # Feature class of station meta data/locations
    elev_tiff = '' # Needed for wind speed
    dem = '' # Needed for almost all
    view_factor = '' # Needed for Thermal radiation
    db = ''
    
    if watershed == 'Johnston Draw':
        arcpy.AddMessage('Johnston Draw Watershed')
        stations = r'C:\ReynoldsCreek\Relevant_Data.gdb\station_locations_JD'
        elev_tiff = r'C:\ReynoldsCreek\jd_elevation_filled.tif'
        dem = r'C:\ReynoldsCreek\Relevant_Data.gdb\JD_DEM_10m'
        view_factor = r'C:\ReynoldsCreek\Relevant_Data.gdb\JD_ViewFactor_10m' 
        db = 'jd_data'
    
    elif watershed == 'Reynolds Creek':
        arcpy.AddMessage('Reynolds Creek Watershed')
        stations = r'C:\ReynoldsCreek\Relevant_Data.gdb\station_locations'
        elev_tiff = r'C:\ReynoldsCreek\rc_elevation_filled.tif'
        dem = r'C:\ReynoldsCreek\Relevant_Data.gdb\RC_DEM_10m_South'
        view_factor = r'C:\ReynoldsCreek\Relevant_Data.gdb\RC_ViewFactor_10M_South'
        db = 'rc_data'
    return stations, elev_tiff, dem, view_factor, db

def ConnectDB(db):
    '''connect to MySQL database'''
    try:
      cnx = mysql.connector.connect(user='root', password='',
                                    host='localhost',
                                    database=db,
                                    buffered=True)
    except mysql.connector.Error as err:

        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            arcpy.AddMessage('Something is wrong with your user name or password')
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            arcpy.AddMessage('Database does not exist')
        else:
            arcpy.AddMessage(err)
    else:
        arcpy.AddMessage('Connection successful')
    return cnx

def ParameterList(param_dict, rows):
    '''Append all data to the end of the parameter list'''
    count = rows.rowcount
    for i in range(0, count):
        row = rows.fetchone()
        if data['watershed'] == 'Johnston Draw':
            param_dict['site_key'].append(row[0])
            param_dict['date_time'].append(row[1])
            param_dict['air_temperature'].append(row[8])
            param_dict['vapor_pressure'].append(row[10])
            param_dict['dew_point'].append(row[11])
            param_dict['solar_radiation'].append(row[12])
            param_dict['wind_speed'].append(row[13])
            param_dict['wind_direction'].append(row[14])
        elif data['watershed'] == 'Reynolds Creek':
            param_dict['site_key'].append(row[0])
            param_dict['date_time'].append(row[1])
            param_dict['air_temperature'].append(row[9])
            param_dict['vapor_pressure'].append(row[11])
            param_dict['dew_point'].append(row[12])
            param_dict['solar_radiation'].append(row[13])
            param_dict['wind_speed'].append(row[14])
            param_dict['wind_direction'].append(row[15])
    return param_dict

def BuildClimateTable(params, num):
    arcpy.management.CreateTable(data['scratch_gdb'], 'climate_table')
    table = data['scratch_gdb'] + '/climate_table'
    keys = [] # Holds data types collected (wind speed, air temperature, etc) to add to table
    for key in params:
        if key == 'site_key':
            ftype = 'TEXT'
        elif key == 'date_time':
            ftype = 'DATE'
        else:
            ftype = 'FLOAT'
        arcpy.management.AddField(in_table = table,
            field_name = key,
            field_type = ftype)
        keys.append(key)
    in_cursor = arcpy.InsertCursor(table)
    # Sorry about this thing.  I meant it to be able to scale according to the
    # type of data colelcted (ie. to make it simple to add different data to the
    # list.  
    for j in range(0, num):
        row = in_cursor.newRow()
        for k in range(0, len(keys)):
            # ---Explained---
            # keys[x] = site_key, air_temperature, etc.
            # params[keys[k][j] = value (ie -2.5)
            row.setValue(keys[k], params[keys[k]][j])
        in_cursor.insertRow(row)
    del in_cursor
    del row
    return table

def DeleteScratchData(in_list):
    for path in in_list:
        arcpy.management.Delete(path)

# Main Function --- Figure out a way to be run as script or as tool
#======================================================================
if __name__ == '__main__':
    #Checkout needed Extentions
    arcpy.CheckOutExtension('Spatial')
    
    data['from_date'] = roundTime(datetime.datetime.strptime(data['from_date'], '%Y-%m-%d %H:%M:%S'), 60*60)
    
    data['to_date'] = roundTime(datetime.datetime.strptime(data['to_date'], '%Y-%m-%d %H:%M:%S'))
    
    return_ws = selectWatershed(data['watershed'])
    data.update({'stations' : return_ws[0], 
        'elev_tiff' : return_ws[1], 
        'dem' : return_ws[2], 
        'view_factor' : return_ws[3], 
        'db' : return_ws[4] 
        })

    # Connect to database
    db_cnx = ConnectDB(data['db'])
    
    # Scratch and output lists
    ls_scratch_data = []
    ls_output = []
    
    #Master stations feature class to be copied for each gridding function
    data.update({'fc_stations_elev': data['scratch_gdb'] + '/stations_wElev'})
    ls_scratch_data.append(data['fc_stations_elev'])
    data.update({'station_locations' : data['scratch_gdb'] + '/station_locations'})
    arcpy.management.CopyFeatures(data['stations'], data['station_locations'])
    ls_scratch_data.append(data['station_locations'])
    data['ext_features'] = arcpy.Describe(data['station_locations']).extent
    
    arcpy.env.cellSize = data['dem']
    arcpy.AddMessage(arcpy.Describe(data['dem']).extent)
    data.update({'output_cell_size' : arcpy.env.cellSize,
        'ext_elev' : arcpy.Describe(data['dem']).extent
        })
    
    arcpy.sa.ExtractValuesToPoints(in_point_features = data['station_locations'],
            in_raster = data['dem'],
            out_point_features = data['fc_stations_elev'], 
            interpolate_values = 'NONE', 
            add_attributes = 'VALUE_ONLY')

    delta = datetime.timedelta(hours=data['time_step'])
    date_increment = data['from_date']
    while date_increment < data['to_date']:
        if any([data['bool_all_tools'], data['bool_air_temperature'],
            data['bool_dew_point'], data['bool_vapor_pressure'],
            data['bool_wind_speed'], data['bool_solar_radiation'],
            data['bool_thermal_radiation']]):
            # Run climate data 
            ls_scratch_data_imd = []

            # Paramter lists
            parameters = {'site_key' : [],
                    'date_time' : [],
                    'air_temperature' : [],
                    'vapor_pressure' : [],
                    'dew_point' : [],
                    'solar_radiation' : [],
                    'wind_speed' : [],
                    'wind_direction' : []
                    }
            
            # Query climage (weather) table
            from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
            time_stamp = date_increment.strftime('%Y%m%d_%H')
            to_date_temp = date_increment + delta
            to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
            query = ("SELECT * FROM weather WHERE "\
                "date_time >= '" + from_date + "' "\
                "AND date_time < '" + to_date + "'")
            cur = db_cnx.cursor()
            cur.execute(query)
            i_num_return = cur.rowcount
            print('Query: ' + query)
            print('Row Count: ' + str(i_num_return))
            #Build parameter lists into dictionary
            parameters = ParameterList(parameters, cur)
            cur.close()
            
            # Build Climate table
            climate_table = BuildClimateTable(parameters, i_num_return)
            ls_scratch_data_imd.append(climate_table)
            # Run interpolation tools
            if data['bool_air_temperature']:
                print('Air Temperature')
            
            DeleteScratchData(ls_scratch_data_imd)
        
        
        date_increment += delta
    DeleteScratchData(ls_scratch_data)
