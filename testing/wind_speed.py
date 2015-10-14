def windSpeed(inDateTime, elevation):
    arcpy.AddMessage("Calculating Wind Speed") 
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

    #Elevation tiff
    elev_file = elevation
    #Set up time parameters
    ninjaPath = "C:/WindNinja/WindNinja-2.5.1/bin/WindNinja_cli.exe"
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
    args = [ninjaPath,
    "--initialization_method", "pointInitialization",
    "--elevation_file", elev_file, #elevation raster (cannot contain any "no-data" values)
    "--match_points", "false", #match simulations to points (simulation fails if set to true)
    "--year", sWindYear,
    "--month", sWindMonth,
    "--day", sWindDay,
    "--hour", sWindHour,
    "--minute", sWindMinute,
    "--mesh_resolution", output_cell_size, #Resolution of model calculations
    "--vegetation", "brush", #Vegetation type (can be 'grass', 'brush', or 'trees')
    "--time_zone", "America/Boise", #time zone of target simulation
    "--diurnal_winds", "true", #consider diurnal cycles in calculations
    "--write_goog_output", "false", #write kml output (boolean: true/false)
    "--write_shapefile_output", "false", #write shapefile output (boolean: true/false)
    "--write_farsite_atm", "false", #write fire behavior file (boolean: true/false)
    "--write_ascii_output", "true", #write ascii file output (this should always be set to true)
    "--ascii_out_resolution", "-1", #resolution of output (-1 means same as mesh_resolution)
    "--units_ascii_out_resolution", "m",
    "--units_mesh_resolution", "m", #units of resolution of model calculations (should be "m" for meters)
    "--units_output_wind_height", "m", #units of output wind height
    "--output_speed_units", "mps",
    "--output_wind_height", "3",
    "--wx_station_filename", scratchWS + "/wnStations.csv", #weather station csv file used in point initialization method
    "--output_path", scratchWS] #path to output

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
            lsScratchData_Imd.append(path2ASCII)
        elif ( file.endswith("_vel.prj") or file.endswith('_ang.asc') or
               file.endswith('_ang.prj') or file.endswith('cld.asc') or 
               file.endswith('_cld.prj') ):
            lsScratchData_Imd.append(scratchWS + '/' + file)
    arcpy.ASCIIToRaster_conversion(in_ascii_file=path2ASCII, out_raster=outFolder + "/wind_speed_" + sTimeStamp + ".tif", data_type="FLOAT")

    #Get coordinate system information
    desc = arcpy.Describe(rc_elevation)
    coordSystem = desc.spatialReference
    arcpy.DefineProjection_management(outFolder + "/wind_speed_" + sTimeStamp + ".tif", coordSystem)

    return outFolder + "/wind_speed_" + sTimeStamp + ".tif"


