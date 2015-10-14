def soilTemperature():
    arcpy.AddMessage("Calculating Soil Temperature") 
    #Set extent to full feature class (all stations)
    arcpy.env.extent = extFullFeatures
    #Caclulate average soil temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/soilTemperature_table1", "soilTemperature_Table", "stm005 > -500")
    arcpy.Statistics_analysis("soilTemperature_Table", "in_memory/soilTemperature_Table2", "stm005 MEAN", "site_key")
    #Copy snow depth values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/soilTemperature_Table2", "site_key", "MEAN_stm005")
    #Delete rows from station feature class that have null values for soil temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_stm005") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
    #Run ordinary least squares on tempStations
    arcpy.CreateTable_management(scratchGDB, "coef_table")
    ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_stm005","RASTERVALU", scratchGDB + "/coef_table","","")
    lsScratchData_Imd.append(scratchGDB + "/coef_table")
    intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
    slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
    #Create final raster
    arcpy.env.extent = extElevation
    output_raster = (Raster(rc_elevation) * slope + intercept)
    output_raster.save(outFolder + "/soil_temperature_" + sTimeStamp + ".tif")



    return outFolder + "/soil_temperature_" + sTimeStamp + ".tif"


