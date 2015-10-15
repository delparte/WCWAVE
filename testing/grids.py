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

#Checkout needed Extentions
arcpy.CheckOutExtension('Spatial')
arcpy.CheckOutExtension('GeoStats')

## data = {'watershed' : arcpy.GetParameterAsText(0),
##     'from_date' : arcpy.GetParameterAsText(1),
##     'to_date' : arcpy.GetParameterAsText(2),
##     'time_step' : int(arcpy.GetParameterAsText(3)),
##     'kriging_method' : arcpy.GetParameterAsText(4),
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
    'to_date': u'2014-01-01 13:00:00', 
    'time_step': 1, 
    'bool_soil_temperature': False, 
    'rl_constant': 0.005, 
    'from_date': u'2014-01-01 12:00:00', 
    'bool_solar_radiation': False, 
    'bool_all_tools': False, 
    'h20_constant': 0.2, 
    'db': 'jd_data', 
    'bool_dew_point': False, 
    'bool_precip_mass': False, 
    'bool_wind_speed': False, 
    'kriging_method': u'Empirical Bayesian Kriging', 
    'bool_thermal_radiation': False, 
    'bool_constants': False, 
    'bool_snow_properties': False, 
    'watershed': u'Johnston Draw', 
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
data['out_folder'] = outFolder
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
        return cnx
    except mysql.connector.Error as err:

        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            arcpy.AddMessage('Something is wrong with your user name or password')
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            arcpy.AddMessage('Database does not exist')
        else:
            arcpy.AddMessage(err)
    else:
        arcpy.AddMessage('Connection successful')

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

def DataTable(parameter, data_table):
    ''' Create paramater scratch table to be used for interpolation '''
    scratch_data = []
    temp_table1 = parameter + '_table'
    out = arcpy.management.MakeTableView(in_table = data_table,
            out_view = temp_table1,
            where_clause = parameter + ' > -500')
    temp_table2 = 'in_memory/' + parameter + '_table2'
    scratch_data.append(temp_table1)
    out_mem = arcpy.analysis.Statistics(in_table = temp_table1,
            out_table = temp_table2,
            statistics_fields = parameter + ' MEAN',
            case_field = 'site_key')
    # Copy stats to tempStations feature class
    temp_stations = arcpy.management.CopyFeatures(in_features = data['fc_stations_elev'],
            out_feature_class = data['scratch_gdb'] + '/tempStations')
    # Join stats to temp stations feature class
    arcpy.management.JoinField(in_data = temp_stations, 
            in_field = 'Site_Key',
            join_table = temp_table2,
            join_field = 'site_key',
            fields = 'MEAN_' + parameter)
    
    # Delete rows from feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(temp_stations)
    for row in cursor:
        if (row.getValue('RASTERVALU') < 0 or
                row.getValue('RASTERVALU') == 'None' or
                row.getValue('RASTERVALU') is None ):
            cursor.deleteRow(row)
        else:
            row.setValue('RASTERVALU', round(row.getValue('RASTERVALU'), 2))
            cursor.updateRow(row)

    del cursor
    del row
    # Delete rows from feature class that have null values for paramter
    cursor = arcpy.UpdateCursor(temp_stations)
    for row in cursor:
        if row.isNull('MEAN_' + parameter):
            cursor.deleteRow(row)
        else:
            row.setValue('MEAN_' + parameter, round(row.getValue('MEAN_' + parameter), 2))
            cursor.updateRow(row)
    del cursor
    del row
    DeleteScratchData(scratch_data)
    return temp_stations

def DetrendedMethod(parameter, data_table, date_stamp):
    print('Detrended Kriging')
    # Add unique ID field to temporary data table for use in OLS function
    arcpy.management.AddField(in_table = data_table,
            field_name = 'Unique_ID',
            field_type = 'SHORT',
            field_is_nullable = 'NULLABLE',
            field_is_required = 'NON_REQUIRED')
    arcpy.management.CalculateField(in_table = data_table,
            field = 'Unique_ID',
            expression = '!OBJECTID!',
            expression_type = 'PYTHON_9.3')
    #Run ordinary least squares of temporary data table
    coef_table = arcpy.management.CreateTable(data['scratch_gdb'], 'coef_table_' + parameter)
    ols = arcpy.stats.OrdinaryLeastSquares(Input_Feature_Class = data_table,
            Unique_ID_Field = 'Unique_ID',
            Output_Feature_Class = 'in_memory/fcResid',
            Dependent_Variable = 'MEAN_' + parameter,
            Explanatory_Variables = 'RASTERVALU',
            Coefficient_Output_Table = coef_table)
    intercept = list((row.getValue('Coef') for row in arcpy.SearchCursor(coef_table, fields='Coef')))[0]
    slope = list((row.getValue('Coef') for row in arcpy.SearchCursor(coef_table, fields='Coef')))[1]
    #Calculate residuals and add them to temporary data table
    arcpy.management.AddField(in_table = data_table,
            field_name = 'residual', 
            field_type = 'DOUBLE',
            field_is_nullable = 'NULLABLE',
            field_is_required = 'NON_REQUIRED')
    cursor = arcpy.UpdateCursor(data_table)
    for row in cursor:
        row_math = row.getValue('MEAN_' + parameter) - ((slope * row.getValue('RASTERVALU')) + intercept)
        row.setValue('residual', row_math)
        cursor.updateRow(row)
    del cursor
    del row
    #Run ordinary kriging on residuals
    k_model = KrigingModelOrdinary('SPHERICAL', 460, 3686, .1214, .2192)
    radius = RadiusFixed(10000, 1)
    outKrig = Kriging(in_point_features = data_table, 
            z_field = 'residual',
            kriging_model = k_model,
            cell_size = data['output_cell_size'],
            search_radius = radius)
    resid_raster = data['scratch_gdb'] + '/' + parameter
    outKrig.save(resid_raster)
    return_raster = arcpy.Raster(resid_raster) + (arcpy.Raster(data['dem']) * slope + intercept)
    return_raster.save(data['out_folder'] + '/' + parameter + '_' + str(date_stamp) + '.tif')
    #Delete scratch/residual data.
    del outKrig
    del k_model
    del radius
    arcpy.management.Delete(resid_raster)
    return return_raster

def EBKMethod(parameter, data_table, date_stamp):
    print('Empirical Bayesian Kriging')
    arcpy.ga.EmpiricalBayesianKriging(in_features = data_table,
            z_field = 'MEAN_' + parameter,
            out_raster = data['scratch_gdb'] + '/' + parameter,
            cell_size = data['output_cell_size'],
            transformation_type = 'EMPIRICAL',
            max_local_points = '100',
            overlap_factor = '1',
            number_semivariograms = '100',
            search_neighborhood = 'NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2',
            output_type = 'PREDICTION',
            quantile_value = '0.5',
            threshold_type = 'EXCEED',
            semivariogram_model_type='WHITTLE_DETRENDED')
    outExtract = ExtractByMask(data['scratch_gdb'] + '/' + parameter, data['dem'])
    outExtract.save(data['out_folder'] + '/' + parameter + '_' + str(date_stamp) + '.tif')
    arcpy.management.Delete(data['scratch_gdb'] + '/' + parameter)
    return data['out_folder'] + '/' + parameter + '_' + str(date_stamp) + '.tif'

def AirTemperature(clim_tab, date_stamp):
    print('Air Temperature')
    param = 'air_temperature'
    scratch_table = DataTable(param, clim_tab)
    arcpy.management.CopyRows(scratch_table, data['scratch_gdb'] + '/temp_ta')
    
    #Interpolate using the chosen method
    if data['kriging_method'] == 'Detrended':
        raster = DetrendedMethod(param, scratch_table, date_stamp)
        #raster.save(data['out_folder'] + '/' + param + '.tif')
    elif data['kriging_method'] == 'Combined':
        raster = CombinedMethod(param, scratch_table, date_stamp)
    else:
        raster = EBKMethod(param, scratch_table, date_stamp)
    
    #Delete tempStations when done.
    arcpy.management.Delete(scratch_table)

def DewPoint(clim_tab, date_stamp):
    print('Dewpoint Temperature')
    param = 'dew_point'
    scratch_table = DataTable(param, clim_tab)
    arcpy.management.CopyRows(scratch_table, data['scratch_gdb'] + '/temp_dp')
    
    #Interpolate using the chosen method
    if data['kriging_method'] == 'Detrended':
        raster = DetrendedMethod(param, scratch_table, date_stamp)
        raster.save(data['out_folder'] + '/' + param + '.tif')
    elif data['kriging_method'] == 'Combined':
        raster = CombinedMethod(param, scratch_table, date_stamp)
    else:
        raster = EBKMethod(param, scratch_table, date_stamp)

    #Delete tempStations when done
    arcpy.management.Delete(scratch_table)

def DeleteScratchData(in_list):
    for path in in_list:
        arcpy.management.Delete(path)
    arcpy.management.Delete('in_memory')

# Main Function --- Figure out a way to be run as script or as tool
#======================================================================
def main():
    from_date_round = datetime.datetime.strptime(data['from_date'], '%Y-%m-%d %H:%M:%S')
    to_date_round = datetime.datetime.strptime(data['to_date'], '%Y-%m-%d %H:%M:%S')
    data['from_date'] = roundTime(from_date_round, 60*60)
    
    data['to_date'] = roundTime(to_date_round)
    
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
    arcpy.env.extent = data['ext_elev']
    
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
                path_air_temp = AirTemperature(climate_table, time_stamp)
                path_dew_point = DewPoint(climate_table, time_stamp)
            DeleteScratchData(ls_scratch_data_imd)
        
        
        date_increment += delta
    DeleteScratchData(ls_scratch_data)

if __name__ == '__main__':
    main() 
