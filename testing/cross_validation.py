import grids
import arcpy
import datetime

import matplotlib.pyplot as plt

def LeaveOneOutValue(raster, station_locations, site_toget):
    print ('Extracting leave one out value for {0}'.format(site_toget))
    with arcpy.da.SearchCursor(in_table = station_locations, 
            field_names = ['Site_Key','SHAPE@XY']) as cursor:
        for row in cursor:
            if row[0] == site_toget:
                point = '{0} {1}'.format(row[1][0], row[1][1])
                value = arcpy.management.GetCellValue(raster, point)
    return value

def ObservedValue(sites_data, parameter, site_toget):
    print('Get observed value from table')
    index = 0
    for name in sites_data['site_key']:
        if name == site_toget:
            obs_value = sites_data[parameter][index]
        index += 1
    return obs_value

def GraphRegression(x, y):
    print('Making Graph')
    print 'X: {0}'.format(x)
    print 'Y: {0}'.format(y)
    z = y-x
    print 'Z: {0}'.format(z)
    z = z**2
    plt.scatter(x,y,s=z, alpha=0.5)
    

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
    arcpy.AddMessage(arcpy.Describe(grids.data['dem']).extent)
    grids.data.update({'output_cell_size' : arcpy.env.cellSize,
        'ext_elev' : arcpy.Describe(grids.data['dem']).extent
        })
    arcpy.env.extent = grids.data['ext_elev']
    
    arcpy.sa.ExtractValuesToPoints(in_point_features = grids.data['station_locations'],
            in_raster = grids.data['dem'],
            out_point_features = grids.data['fc_stations_elev'],
            interpolate_values = 'NONE',
            add_attributes = 'VALUE_ONLY')
    print grids.data
    delta = datetime.timedelta(hours=grids.data['time_step'])
    date_increment = grids.data['from_date']
    
    while date_increment < grids.data['to_date']:
        if any([grids.data['bool_air_temperature']]):
            print date_increment
            
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
            cur.close()
            
            # get sites then reset parameters dictionary
            sites_list = parameters['site_key']
            observed = []
            modeled = []
            
            for site in sites_list:
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
            GraphRegression(modeled, observed)
        date_increment += delta
    print grids.data


if __name__ == '__main__':
    grids.data.update({'from_date' : u'2014-01-01 12:00:00',
        'to_date' : u'2014-01-01 13:00:00',
        'bool_air_temperature' : True,
        'watershed' : 'Johnston Draw',
        'time_step' : 1,
        'kriging_method' : 'Empirical Bayesian'
        })
    main()
