import arcpy
import datetime
import time
import math
import numpy as np

import matplotlib.pyplot as plt

import sys
sys.path.append('..')
import gridtools as grids

arcpy.env.parallelProcessingFactor = '100%'

def LeaveOneOutValue(raster, station_locations, site_toget):
    arcpy.AddMessage('Extracting leave one out value for {0}'.format(site_toget))
    with arcpy.da.SearchCursor(in_table = station_locations, 
            field_names = ['Site_Key','SHAPE@XY']) as cursor:
        for row in cursor:
            if row[0] == site_toget:
                point = '{0} {1}'.format(row[1][0], row[1][1])
                result = arcpy.management.GetCellValue(raster, point)
                value = float(result.getOutput(0))
    return value

def ObservedValue(sites_data, parameter, site_toget):
    arcpy.AddMessage('Get observed value from table')
    index = 0
    for name in sites_data['site_key']:
        if name == site_toget:
            obs_value = sites_data[parameter][index]
        index += 1
    return obs_value

def GraphRegression(time_stamp='test', param_type = 'test', 
        label=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i','j','k'],
        x = [3.07, 3.49, 6.77, 3.87, 4.90, 4.64, 3.94, 4.93, 3.99, 4.46, 3.36], 
        y = [2.3, 3.1, 4.3, 5.3, 3.9, 5.2, 4.1, 5.0, 3.9, 5.3, 3.9] ):
    arcpy.AddMessage('Making {0} Graph'.format(param_type))
    arcpy.AddMessage('X: {0}'.format(x))
    arcpy.AddMessage('Y: {0}'.format(y))
    abs_error = []
    sqr_error = []
    z = []
    for i in range(len(x)):
        z.append(y[i] - x[i]) # Observed minus modeled
        abs_error.append(abs(y[i]-x[i]))
        sqr_error.append((y[i]-x[i])**2)
        z[i] = abs(z[i]) * 100
    print z
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
    elif param_type == 'ppts':
        point_color = 'blue'
    elif param_type == 'zs':
        point_color = 'yellow'

    trendline = np.polyfit(x, y, 1)
    p = np.poly1d(trendline)
    plt.plot(x, p(x), 'k-')
    
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

def PrintDataToCSV(parameter, sites_names=['test'], modeled_values=[1], observed_values=[2], time_values=[3]):
    filename = open('{0}\Output_{1}.txt'.format(grids.outFolder, parameter), 'a')
    filename.write('site_key,modeled_value,observed_value,time_to_create (s)\n')
    for i in range(len(sites_names)):
        string = '{0},{1},{2},{3},\n'.format(sites_names[i],modeled_values[i],observed_values[i], time_values[i])
        filename.write(string)
    

def main():
    from_date_round = datetime.datetime.strptime(grids.data['from_date'], '%Y-%m-%d %H:%M:%S')
    to_date_round = datetime.datetime.strptime(grids.data['to_date'], '%Y-%m-%d %H:%M:%S')
    grids.data['from_date'] = grids.roundTime(from_date_round, 60*60)
    grids.data['to_date'] = grids.roundTime(to_date_round)
    
    return_ws = grids.selectWatershed(grids.data['watershed'])
    grids.data.update({'stations' : return_ws[0],
        'elev_tiff' : return_ws[1],
        'dem' : return_ws[2],
        'view_factor' : return_ws[3],
        'db' : return_ws[4]
        })
    
    db_cnx = grids.ConnectDB(grids.data['db'])
    # Scratch and output lists
    ls_scratch_data = []
    ls_output = []
    
    #Master stations feature class to be copied for each gridding function
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
    
    arcpy.sa.ExtractValuesToPoints(in_point_features = grids.data['station_locations'],
            in_raster = grids.data['dem'],
            out_point_features = grids.data['fc_stations_elev'],
            interpolate_values = 'NONE',
            add_attributes = 'VALUE_ONLY')
    delta = datetime.timedelta(hours=grids.data['time_step'])
    date_increment = grids.data['from_date']
    print grids.data['fc_stations_elev']
    scursor = arcpy.SearchCursor(grids.data['fc_stations_elev'])
    station_welevation = []
    for row in scursor:
        station = row.getValue('Site_Key')
        station_welevation.append(station)
    while date_increment < grids.data['to_date']:
        arcpy.AddMessage(date_increment)
        if any([grids.data['bool_air_temperature']]):
            ls_scratch_data_imd = []
            
            parameters = {'site_key' : [],
                    'date_time' : [],
                    'air_temperature' : [],
                    'vapor_pressure' : [],
                    'dew_point' : [],
                    'solar_radiation' : [],
                    'wind_speed' : [],
                    'wind_direction' : []
                    }
            
            # query climate (weather) table
            from_date = date_increment.strftime('%Y-%m-%d %H:%M:%S')
            time_stamp = date_increment.strftime('%Y%m%d_%H')
            to_date_temp = date_increment + delta
            to_date = to_date_temp.strftime('%Y-%m-%d %H:%M:%S')
            query = ('SELECT * FROM weather WHERE '\
                    'date_time >= "' + from_date + '" '\
                    'AND date_time < "' + to_date + '"')
            cur = db_cnx.cursor()
            cur.execute(query)
            i_num_return = cur.rowcount
            parameters = grids.ParameterList(parameters, cur, table_type = 'climate')
            obs_parameters = parameters
            print obs_parameters
            cur.close()
            
            # get sites then reset parameters dictionary
            sites_list = []
            for st in parameters['site_key']:
                if st in station_welevation:
                    sites_list.append(st)
            print 'sites_lists len: {0}'.format(len(sites_list))
            observed = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []}
            modeled = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []} 
            create_time = {'air_temperature': [], 'dew_point': [], 'vapor_pressure': []}
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
                        'date_time >= "' + from_date + '" '\
                        'AND date_time < "' + to_date + '" '\
                        'AND Site_Key <> "' + site + '"')
                cur = db_cnx.cursor()
                cur.execute(query)
                i_num_return = cur.rowcount
                parameters = grids.ParameterList(parameters,cur,table_type = 'climate')
                cur.close()
                
                # Build Climate table
                climate_table = grids.BuildClimateTable(parameters, i_num_return)
                ls_scratch_data_imd.append(climate_table)
                # Run interpolation tools
                start = time.time()
                if grids.data['bool_air_temperature']:
                    path_air_temp = grids.AirTemperature(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_air_temp, 
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'air_temperature', site)
                    modeled['air_temperature'].append(point_error)
                    observed['air_temperature'].append(obs_point)
                end_air = time.time()
                create_time['air_temperature'].append(end_air - start)
                start = time.time()
                if grids.data['bool_dew_point']:
                    path_dew_point = grids.DewPoint(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_dew_point,
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'dew_point', site)
                    modeled['dew_point'].append(point_error)
                    observed['dew_point'].append(obs_point)
                end_dew = time.time()
                create_time['dew_point'].append(end_dew - start)
                if grids.data['bool_vapor_pressure']:
                    path_vapor_pressure = grids.DewPoint(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_vapor_pressure,
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'vapor_pressure', site)
                    modeled['vapor_pressure'].append(point_error)
                    observed['vapor_pressure'].append(obs_point)
                end_dew = time.time()
                create_time['vapor_pressure'].append(end_dew - start)

            GraphRegression(time_stamp = time_stamp,
                    label = sites_list,
                    param_type = 'T_a',
                    x = modeled['air_temperature'],
                    y = observed['air_temperature'])
##             GraphRegression(time_stamp = time_stamp,
##                     label = sites_list,
##                     param_type = 'T_pp',
##                     x = modeled['dew_point'],
##                     y = observed['dew_point'])
##             GraphRegression(time_stamp = time_stamp,
##                     label = sites_list,
##                     param_type = 'e_a',
##                     x = modeled['vapor_pressure'],
##                     y = observed['vapor_pressure'])

##         GraphRegression()
            PrintDataToCSV('air_temp', sites_list, modeled['air_temperature'], observed['air_temperature'], create_time['air_temperature'])
            
        date_increment += delta


if __name__ == '__main__':
##     grids.data.update({'from_date' : u'2007-12-01 00:00:00',
##         'to_date' : u'2007-12-03 00:00:00',
##         'bool_air_temperature' : True,
##         'bool_dew_point': False,
##         'bool_vapor_pressure': False,
##         'watershed' : 'Reynolds Creek',
##         'time_step' : 1,
##         'kriging_method' : 'Empirical Bayesian'
##         })
    grids.data.update({'from_date' : u'2014-01-01 00:00:00',
        'to_date' : u'2014-01-01 01:00:00',
        'bool_air_temperature' : True,
        'bool_dew_point': False,
        'bool_vapor_pressure': False,
        'watershed' : 'Johnston Draw',
        'time_step' : 1,
        'kriging_method' : 'Empirical Bayesian'
        })
    main()
