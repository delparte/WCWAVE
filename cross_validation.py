import arcpy
import datetime
import time
import math
import numpy as np
import os

import matplotlib.pyplot as plt

import gridtools as grids

arcpy.env.parallelProcessingFactor = '100%'

def LeaveOneOutValue(raster, station_locations, site_toget):
    arcpy.AddMessage('Extracting modeled value for {0}'.format(site_toget))
    with arcpy.da.SearchCursor(in_table = station_locations, 
            field_names = ['Site_Key','SHAPE@XY']) as cursor:
        for row in cursor:
            if row[0] == site_toget:
                point = '{0} {1}'.format(row[1][0], row[1][1])
                result = arcpy.management.GetCellValue(raster, point)
                value = float(result.getOutput(0))
    return value

def ObservedValue(sites_data, parameter, site_toget):
    arcpy.AddMessage('Get observed value from table for {0}'.format(site_toget))
    index = 0
    for name in sites_data['site_key']:
        if name == str(site_toget):
            obs_value = sites_data[parameter][index]
        index += 1
    return obs_value

def GraphRegression(time_stamp='test', param_type = 'test', 
        label=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i','j','k'],
        x = [3.07, 3.49, 6.77, 3.87, 4.90, 4.64, 3.94, 4.93, 3.99, 4.46, 3.36], 
        y = [2.3, 3.1, 4.3, 5.3, 3.9, 5.2, 4.1, 5.0, 3.9, 5.3, 3.9] ):
    arcpy.AddMessage('Making {0} Graph'.format(param_type))
    abs_error = []
    sqr_error = []
    z = []
    for i in range(len(x)):
        z.append(y[i] - x[i]) # Observed minus modeled
        abs_error.append(abs(y[i]-x[i]))
        sqr_error.append((y[i]-x[i])**2)
        z[i] = abs(z[i]) * 100
##     print z
    mae = sum(abs_error)/len(abs_error)
    mae = round(mae, 2)
    rmse = math.sqrt(sum(sqr_error)/len(sqr_error))
    rmse = round(rmse, 2)
    
    #Select which color the points will be based on parameter type (air temp, dew point, etc.)
    if param_type == 'test' or param_type == 'T_a':
        point_color = 'gray'
    elif param_type == 'T_pp':
        point_color = 'green'
    elif param_type == 'e_a':
        point_color = 'red'
    elif param_type == 'm_pp':
        point_color = 'blue'
    elif param_type == 'z_s':
        point_color = 'yellow'
    #
    x_y_line = [min(x), max(x)]
    y_values = x_y_line
    #trendline = np.polyfit(x, y, 1)
    #p = np.poly1d(trendline)
    plt.plot(x_y_line, y_values, 'k-')
    
    plt.title('{0} for {1}; RMSE={2}; MAE={3}'.format(param_type, time_stamp,rmse, mae))
    plt.scatter(x,y,s=z, c=point_color)
    plt.xlabel('Modeled')
    plt.ylabel('Observed')
    i=0
    for point in label:
        plt.annotate(point, xy = (x[i],y[i]), xytext = (3,10),
                textcoords = 'offset points', 
                bbox=dict(boxstyle = 'round,pad=0.2', fc='yellow', alpha = 0.5))
        i+=1
##     plt.show()
    plt.savefig('{0}\{1}_scatter_{2}.png'.format(grids.outFolder,param_type, time_stamp))
    plt.clf()

def PrintDataToCSV(parameter, date_value, sites_names=['test'], modeled_values=[1], observed_values=[2], time_values=[3]):
    out_file_name = '{0}\Output_{1}.txt'.format(grids.outFolder, parameter)

    # If file does not exists write headers.
    if os.path.isfile(out_file_name):
        filename = open(out_file_name, 'a')
    else: 
        filename = open(out_file_name, 'a')
        filename.write('site_key,date,modeled_value,observed_value,time_to_create (s)\n')
    
    # Write data to file
    for i in range(len(sites_names)):
        string = '{0},{1},{2},{3},{4},\n'.format(sites_names[i],date_value, modeled_values[i],observed_values[i], time_values[i])
        filename.write(string)
    

def main():
    # Initialize datetime, watershed, and database variables
    # ===============================================================
    from_date_round = datetime.datetime.strptime(grids.data['from_date'], '%Y-%m-%d %H:%M:%S')
    to_date_round = datetime.datetime.strptime(grids.data['to_date'], '%Y-%m-%d %H:%M:%S')
    grids.data['from_date'] = grids.roundTime(from_date_round, 60*60)
    grids.data['to_date'] = grids.roundTime(to_date_round)
    
    return_ws = grids.selectWatershed(grids.data['watershed']) # Necessary rasters and vector data 
    grids.data.update({'stations' : return_ws[0],
        'elev_tiff' : return_ws[1],
        'dem' : return_ws[2],
        'view_factor' : return_ws[3],
        'db' : return_ws[4]
        })
    
    db_cnx = grids.ConnectDB(grids.data['db']) # Database connection
    ls_scratch_data = []
    ls_output = []
    
    # Master stations feature class to be copied for each gridding function
    #    Stations with elevation and station locations
    # =====================================================================
    grids.data.update({'fc_stations_elev': grids.data['scratch_gdb'] + '/stations_wElev'})
    ls_scratch_data.append(grids.data['fc_stations_elev'])
    grids.data.update({'station_locations' : grids.data['scratch_gdb'] + '/station_locations'})
    arcpy.management.CopyFeatures(grids.data['stations'], grids.data['station_locations'])
    ls_scratch_data.append(grids.data['station_locations'])
    grids.data['ext_features'] = arcpy.Describe(grids.data['station_locations']).extent
    
    arcpy.env.cellSize = grids.data['dem']
    grids.data.update({'output_cell_size' : arcpy.env.cellSize,
        'ext_elev' : arcpy.Describe(grids.data['dem']).extent
        })
    arcpy.env.extent = grids.data['ext_elev']

    # Extract elevations from dem and add to station locations 
    arcpy.sa.ExtractValuesToPoints(in_point_features = grids.data['station_locations'],
            in_raster = grids.data['dem'],
            out_point_features = grids.data['fc_stations_elev'],
            interpolate_values = 'NONE',
            add_attributes = 'VALUE_ONLY')

    # Setup for while loop increments. Create delta variable.  
    #    Create station lists to iterate through
    # ===============================================================
    delta = datetime.timedelta(hours=grids.data['time_step'])
    date_increment = grids.data['from_date']
    print grids.data['fc_stations_elev']
    scursor = arcpy.SearchCursor(grids.data['fc_stations_elev'])
    station_welevation = []
    for row in scursor: # Station list for iteration
        station = row.getValue('Site_Key')
        station_welevation.append(station)

    # While loop run for every time step between from data and to date.
    # ==================================================================
    while date_increment < grids.data['to_date']:
        arcpy.AddMessage('Current timestep: {0}\n'.format(date_increment)) # Print date of current time step
        
        ls_scratch_data_imd = []
        # Air Temp, dew point, or vapor pressure.
        # =======================================
        if any([grids.data['bool_air_temperature'], 
            grids.data['bool_dew_point'], 
            grids.data['bool_vapor_pressure']]):
            
            parameters = {'site_key' : [],
                    'date_time' : [],
                    'air_temperature' : [],
                    'vapor_pressure' : [],
                    'dew_point' : [],
                    'solar_radiation' : [],
                    'wind_speed' : [],
                    'wind_direction' : []
                    }
            
            # query climate (weather) table to get all sites with elevation and data
            # =====================================================================
            from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
            time_stamp = date_increment.strftime('%Y%m%d_%H')
            to_date_temp = date_increment + delta
            to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
            query = ('SELECT * FROM weather WHERE '\
                    'date_time >= ' + grids.data['sql_ph'] + ' '\
                    'AND date_time < ' + grids.data['sql_ph'] + '')
            cur = db_cnx.cursor()
            cur.execute(query, (from_date, to_date))
            rows = cur.fetchall()
            i_num_return = len(rows)
            parameters = grids.ParameterList(parameters, rows, table_type = 'climate')
            obs_parameters = parameters
            cur.close()
            
            # Populate a list of sites with necessary data.
            sites_list = []
            for st in parameters['site_key']:
                if st in station_welevation and st not in sites_list:
                    sites_list.append(st)
            ## arcpy.AddMessage('sites_lists : {0}'.format(sites_list))
            observed = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []}
            modeled = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []} 
            create_time = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []}
            
            # Iterate through sites_list leaving one site out per iteration.
            # ===============================================================
            for site in sites_list:
                arcpy.AddMessage('Leave out {0}'.format(site))
                parameters = {'site_key' : [],
                    'date_time' : [],
                    'air_temperature' : [],
                    'vapor_pressure' : [],
                    'dew_point' : [],
                    'solar_radiation' : [],
                    'wind_speed' : [],
                    'wind_direction' : []
                    }
                query = ('SELECT * FROM weather WHERE '\
                        'date_time >= ' + grids.data['sql_ph'] + ' '\
                        'AND date_time < ' + grids.data['sql_ph'] + ' '\
                        'AND Site_Key <> ' + grids.data['sql_ph'] + '')
                cur = db_cnx.cursor()
                cur.execute(query, (from_date, to_date, site))
                rows = cur.fetchall()
                i_num_return = len(rows)
                ##arcpy.AddMessage('Number of rows: '.format(i_num_return))
                parameters = grids.ParameterList(parameters,rows,table_type = 'climate')
                cur.close()
                
                # Build Climate table
                climate_table = grids.BuildClimateTable(parameters, i_num_return)
                ls_scratch_data_imd.append(climate_table)
                
                # Run interpolation tools
                #    Interpolates data, Retrieves observed data, 
                #    retrieves modeled data, and creation time.
                # ======================================================
                if grids.data['bool_air_temperature']:
                    start = time.time()
                    path_air_temp = grids.AirTemperature(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_air_temp, 
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'air_temperature', site)
                    modeled['air_temperature'].append(point_error)
                    observed['air_temperature'].append(obs_point)
                    end_air = time.time()
                    create_time['air_temperature'].append(end_air - start)
                    ls_scratch_data_imd.append(path_air_temp)
                if grids.data['bool_dew_point']:
                    start = time.time()
                    path_dew_point = grids.DewPoint(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_dew_point,
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'dew_point', site)
                    modeled['dew_point'].append(point_error)
                    observed['dew_point'].append(obs_point)
                    end_dew = time.time()
                    create_time['dew_point'].append(end_dew - start)
                    ls_scratch_data_imd.append(path_dew_point)
                if grids.data['bool_vapor_pressure']:
                    start = time.time()
                    path_vapor_pressure = grids.VaporPressure(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_vapor_pressure,
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'vapor_pressure', site)
                    modeled['vapor_pressure'].append(point_error)
                    observed['vapor_pressure'].append(obs_point)
                    end_dew = time.time()
                    create_time['vapor_pressure'].append(end_dew - start)
                    ls_scratch_data_imd.append(path_vapor_pressure)
                arcpy.AddMessage(' ') # Newline space for clarity
                # END INTERPOLATION
                # ========================================================

            # Graph results and print to Output_[PARAM].txt
            # =====================================================
            if grids.data['bool_air_temperature']:
                GraphRegression(time_stamp = time_stamp,
                        label = sites_list,
                        param_type = 'T_a',
                        x = modeled['air_temperature'],
                        y = observed['air_temperature'])
                PrintDataToCSV('air_temp',
                        date_increment,
                        sites_list,
                        modeled['air_temperature'],
                        observed['air_temperature'],
                        create_time['air_temperature'])
            if grids.data['bool_dew_point']:
                GraphRegression(time_stamp = time_stamp,
                        label = sites_list,
                        param_type = 'T_pp',
                        x = modeled['dew_point'],
                        y = observed['dew_point'])
                PrintDataToCSV('dew_point', 
                        date_increment, 
                        sites_list, 
                        modeled['dew_point'], 
                        observed['dew_point'], 
                        create_time['dew_point'])
            if grids.data['bool_vapor_pressure']:
                GraphRegression(time_stamp = time_stamp,
                        label = sites_list,
                        param_type = 'e_a',
                        x = modeled['vapor_pressure'],
                        y = observed['vapor_pressure'])
                PrintDataToCSV('vapor_pressure',
                        date_increment, 
                        sites_list, 
                        modeled['vapor_pressure'],
                        observed['vapor_pressure'],
                        create_time['vapor_pressure'])
            # END GRAPHING
            # ======================================================
        # Precipitation Mass
        # ========================================
        if any([grids.data['bool_precip_mass']]):
            parameters = {'site_key': [],
                'ppts' : [],
                'pptu' : [],
                'ppta' : []}

            # query precipitation table to get all sites with elevation and data
            # =====================================================================
            from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
            time_stamp = date_increment.strftime('%Y%m%d_%H')
            to_date_temp = date_increment + delta
            to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
            query = ('SELECT * FROM precipitation WHERE '\
                    'date_time >= ' + grids.data['sql_ph'] + ' '\
                    'AND date_time < ' + grids.data['sql_ph'] + '')
            cur = db_cnx.cursor()
            cur.execute(query, (from_date, to_date))
            rows = cur.fetchall()
            i_num_return = len(rows)
            parameters = grids.ParameterList(parameters, rows, table_type = 'precip')
            obs_parameters = parameters
            cur.close()

            # Populate a list of sites with necessary data.
            sites_list = []
            for st in parameters['site_key']:
                if st in station_welevation and st not in sites_list:
                    sites_list.append(st)
            print('Sites_list len: {0}'.format(len(sites_list)))
            observed = {'precip' : []}
            modeled = {'precip' : []}
            create_time = {'precip' : []}
            
            # Iterate through site_list leaving one site out per iteration.
            # ==============================================================
            for site in sites_list:
                arcpy.AddMessage('Leave out {0}'.format(site))
                parameters = {'site_key': [],
                    'ppts' : [],
                    'pptu' : [],
                    'ppta' : []}

                # query precipitation table to get all sites with elevation and data
                # =====================================================================
                from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
                time_stamp = date_increment.strftime('%Y%m%d_%H')
                to_date_temp = date_increment + delta
                to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
                query = ('SELECT * FROM precipitation WHERE '\
                        'date_time >= ' + grids.data['sql_ph'] + ' '\
                        'AND date_time < ' + grids.data['sql_ph'] + ' '\
                        'AND Site_Key <> ' + grids.data['sql_ph'] + '')
                cur = db_cnx.cursor()
                cur.execute(query, (from_date, to_date, site))
                rows = cur.fetchall()
                i_num_return = len(rows)
                parameters = grids.ParameterList(parameters, rows, table_type = 'precip')
                cur.close()
                 
                # Build precipitation table
                climate_table = grids.BuildClimateTable(parameters, i_num_return)
                ls_scratch_data_imd.append(climate_table)
                
                if grids.data['bool_precip_mass']:
                    start = time.time()
                    path_precip_mass = grids.PrecipitationMass(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_precip_mass, 
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'ppts', site)
                    modeled['precip'].append(point_error)
                    observed['precip'].append(obs_point)
                    end_precip = time.time()
                    ls_scratch_data_imd.append(path_precip_mass)
                    create_time['precip'].append(end_precip - start) 
                # END PRECIP INTERPOLATION
                # ===========================================================
                
            # Graph the results and print to Output_[PARAM].txt
            # =====================================================
            if grids.data['bool_precip_mass']:
                GraphRegression(time_stamp = time_stamp,
                        label = sites_list,
                        param_type = 'm_pp',
                        x = modeled['precip'],
                        y = observed['precip'])
                PrintDataToCSV('precip',
                        date_increment,
                        sites_list,
                        modeled['precip'],
                        observed['precip'],
                        create_time['precip'])
            # END Precip graphing and output
            # =======================================================
        if any([grids.data['bool_snow_depth']]):
            parameters = {'site_key' : [],
                    'zs': []}

            # query snow_depth table to get all sites with elevation and data
            # =====================================================================
            from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
            time_stamp = date_increment.strftime('%Y%m%d_%H')
            to_date_temp = date_increment + delta
            to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
            query = ('SELECT * FROM snow_depth WHERE '\
                    'date_time >= ' + grids.data['sql_ph'] + ' '\
                    'AND date_time < ' + grids.data['sql_ph'] + '')
            cur = db_cnx.cursor()
            cur.execute(query, (from_date, to_date))
            rows = cur.fetchall()
            i_num_return = len(rows)
            parameters = grids.ParameterList(parameters, rows, table_type = 'snow_depth')
            obs_parameters = parameters
            arcpy.AddMessage(obs_parameters)
            cur.close()
            
            # Populate a list of sites with necessary data.
            sites_list = []
            for st in parameters['site_key']:
                if st in station_welevation and st not in sites_list:
                    sites_list.append(st)
            print('Sites_list : {0}'.format(sites_list))
            observed = {'snow_depth' : []}
            modeled = {'snow_depth' : []}
            create_time = {'snow_depth' : []}

            for site in sites_list:
                arcpy.AddMessage('Leave out {0}'.format(site))
                parameters = {'site_key': [],
                    'zs' : []}

                # query precipitation table to get all sites with elevation and data
                # =====================================================================
                from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
                time_stamp = date_increment.strftime('%Y%m%d_%H')
                to_date_temp = date_increment + delta
                to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
                query = ('SELECT * FROM snow_depth WHERE '\
                        'date_time >= ' + grids.data['sql_ph'] + ' '\
                        'AND date_time < ' + grids.data['sql_ph'] + ' '\
                        'AND Site_Key <> ' + grids.data['sql_ph'] + '')
                cur = db_cnx.cursor()
                cur.execute(query, (from_date, to_date, site))
                rows = cur.fetchall()
                i_num_return = len(rows)
                parameters = grids.ParameterList(parameters, rows, table_type = 'snow_depth')
                cur.close()
                 
                # Build precipitation table
                climate_table = grids.BuildClimateTable(parameters, i_num_return)
                ls_scratch_data_imd.append(climate_table)
                
                if grids.data['bool_snow_depth']:
                    start = time.time()
                    path_snow_depth = grids.SnowDepth(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_snow_depth, 
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'zs', site)
                    modeled['snow_depth'].append(point_error)
                    observed['snow_depth'].append(obs_point)
                    end_snow = time.time()
                    ls_scratch_data_imd.append(path_snow_depth)
                    create_time['snow_depth'].append(end_snow - start) 
                # END PRECIP INTERPOLATION
                # ===========================================================
            # Graph the results and print to Output_[PARAM].txt
            # =====================================================
            if grids.data['bool_snow_depth']:
                GraphRegression(time_stamp = time_stamp,
                        label = sites_list,
                        param_type = 'z_s',
                        x = modeled['snow_depth'],
                        y = observed['snow_depth'])
                PrintDataToCSV('snow_depth',
                        date_increment,
                        sites_list,
                        modeled['snow_depth'],
                        observed['snow_depth'],
                        create_time['snow_depth'])
            # END Precip graphing and output
            # =======================================================

        date_increment += delta
        grids.DeleteScratchData(ls_scratch_data_imd)


if __name__ == '__main__':
##     # Initialize data variables for Reynolds Creek run.
##     # ========================================================
##     grids.data.update({'from_date' : u'2007-12-03 23:00:00',
##         'to_date' : u'2007-12-04 00:00:00',
##         'bool_air_temperature' : False,
##         'bool_dew_point': False,
##         'bool_vapor_pressure': False,
##         'bool_precip_mass': False,
##         'bool_snow_depth' : True,
##         'watershed' : 'Reynolds Creek',
##         'time_step' : 1,
##         'kriging_method' : 'Empirical Bayesian'
##         })
##
##     # Initialize data variables for Johnston Draw run.
##     # =========================================================
##     grids.data.update({'from_date' : u'2014-01-01 12:00:00',
##         'to_date' : u'2014-01-01 14:00:00',
##         'bool_air_temperature' : False,
##         'bool_dew_point': False,
##         'bool_vapor_pressure': False,
##         'bool_precip_mass': False,
##         'bool_snow_depth' : True,
##         'watershed' : 'Johnston Draw',
##         'time_step' : 1,
##         'kriging_method' : 'Empirical Bayesian'
##         })
##
##     # Initialize data variables for TESTING run.
##     # =========================================================
##     grids.data.update({'from_date' : u'2014-01-01 12:00:00',
##         'to_date' : u'2014-01-01 13:00:00',
##         'bool_air_temperature' : False,
##         'bool_dew_point': False,
##         'bool_vapor_pressure': False,
##         'bool_precip_mass': False,
##         'bool_snow_depth' : True,
##         'watershed' : 'TESTING',
##         'time_step' : 1,
##         'kriging_method' : 'Empirical Bayesian'
##         })
##
    # Initialize data variables for ArcMap run.
    # Comment and uncomment above to run as script
    # =========================================================
    grids.data.update({'watershed' : arcpy.GetParameterAsText(0),
        'from_date' : arcpy.GetParameterAsText(1),
        'to_date' : arcpy.GetParameterAsText(2),
        'time_step' : int(arcpy.GetParameterAsText(3)),
        'kriging_method' : arcpy.GetParameterAsText(4),
        'bool_air_temperature' : arcpy.GetParameter(5),
        'bool_dew_point' : arcpy.GetParameter(6),
        'bool_vapor_pressure': arcpy.GetParameter(7),
        'bool_snow_depth': arcpy.GetParameter(8),
        'bool_precip_mass': arcpy.GetParameter(9)
        })
    main()
    arcpy.AddMessage('\nAll Data output to {0}\n'.format(grids.data['out_folder']))
