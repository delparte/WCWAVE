import grids
import arcpy
import datetime
import time

import matplotlib.pyplot as plt

def LeaveOneOutValue(raster, station_locations, site_toget):
    print ('Extracting leave one out value for {0}'.format(site_toget))
    with arcpy.da.SearchCursor(in_table = station_locations, 
            field_names = ['Site_Key','SHAPE@XY']) as cursor:
        for row in cursor:
            if row[0] == site_toget:
                point = '{0} {1}'.format(row[1][0], row[1][1])
                result = arcpy.management.GetCellValue(raster, point)
                value = float(result.getOutput(0))
    return value

def ObservedValue(sites_data, parameter, site_toget):
    print('Get observed value from table')
    index = 0
    for name in sites_data['site_key']:
        if name == site_toget:
            obs_value = sites_data[parameter][index]
        index += 1
    return obs_value

def GraphRegression(time_stamp, label=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i','j','k'],
        x = [3.07, 3.49, 6.77, 3.87, 4.90, 4.64, 3.94, 4.93, 3.99, 4.46, 3.36], 
        y = [2.3, 3.1, 4.3, 5.3, 3.9, 5.2, 4.1, 5.0, 3.9, 5.3, 3.9] ):
    print('Making Graph')
    print('X: {0}'.format(x))
    print('Y: {0}'.format(y))
    abs_error = []
    sqr_error = []
    z = []
    for i in range(len(x)):
        z.append(y[i] - x[i])
        abs_error.append(abs(y[i]-x[i]))
        sqr_error.append(abs(y[i]-x[i])**2)
        z[i] = (z[i]**2) * 150 + 1
    mae = sum(abs_error)/len(abs_error)
    mae = round(mae, 2)
    rmse = sum(sqr_error)/len(sqr_error)
    rmse = round(rmse, 2)
    plt.title('Data for {0}; RMSE={1}; MAE={2}'.format(time_stamp,rmse, mae))
    plt.scatter(x,y,s=z)
    plt.xlabel('Modeled')
    plt.ylabel('Observed')
    i=0
    for point in label:
        plt.annotate(point, xy = (x[i],y[i]), xytext = (3,10),
                textcoords = 'offset points', 
                bbox=dict(boxstyle = 'round,pad=0.2', fc='yellow', alpha = 0.5))
        i+=1
##     plt.show()
    plt.savefig('{0}\scatter_{1}.png'.format(grids.outFolder,time_stamp))
    plt.clf()

def PrintDataToCSV(sites_names=['test'], modeled_values=[1], observed_values=[2], time_values=[3]):
    filename = open(grids.outFolder + '\Output.txt', 'a')
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
            observed = []
            modeled = []
            create_time = []
            for site in sites_list:
                print('Leave out {0}'.format(site))
                start = time.time()
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
                if grids.data['bool_air_temperature']:
                    path_air_temp = grids.AirTemperature(climate_table, time_stamp)
                    point_error = LeaveOneOutValue(raster = path_air_temp, 
                            site_toget = site,
                            station_locations = grids.data['station_locations'])
                    obs_point = ObservedValue(obs_parameters, 'air_temperature', site)
                    modeled.append(point_error)
                    observed.append(obs_point)
                end_air = time.time()
                create_time.append(end_air - start)
##             print create_time
            GraphRegression(time_stamp = time_stamp,
                    label = sites_list,
                    x = modeled,
                    y = observed)
##         GraphRegression()
            PrintDataToCSV(sites_list, modeled, observed, create_time)
        date_increment += delta


if __name__ == '__main__':
    grids.data.update({'from_date' : u'2007-12-01 00:00:00',
        'to_date' : u'2007-12-02 00:00:00',
        'bool_air_temperature' : True,
        'watershed' : 'Reynolds Creek',
        'time_step' : 1,
        'kriging_method' : 'Empirical Bayesian'
        })
    main()
