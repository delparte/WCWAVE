def solarRadiation():
    arcpy.AddMessage("Calculating Solar Radiation") 
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
