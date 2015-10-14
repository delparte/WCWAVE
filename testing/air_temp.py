def airTemperature():
    """ Kriging """
    arcpy.AddMessage("Calculating Air Temperature")
    #Caclulate average air temperatures over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "airTemperature_Table", "air_temperature > -500")
    arcpy.Statistics_analysis("airTemperature_Table", "in_memory/airTemperature_Table2", "air_temperature MEAN", "site_key")
    #Copy air temperature values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/airTemperature_Table2", "site_key", "MEAN_air_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for air temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    if sKrigMethod == "Combined":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_air_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_air_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
            out_raster=scratchGDB + "/airTemperature_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
        lsScratchData_Imd.append(scratchGDB + "/airTemperature_residual")
        
        outExtractByMask = ExtractByMask(scratchGDB + '/airTemperature_residual', rc_elevation)
        outExtractByMask.save(scratchGDB + '/airTemperature_scratch')
        lsScratchData_Imd.append(scratchGDB + '/airTemperature_scratch')

        output_raster = arcpy.Raster(scratchGDB + "/airTemperature_scratch") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/air_temperature_" + sTimeStamp + ".tif")
    #Check if interpolation method is Empirical Bayesian or Detrended
    elif sKrigMethod == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_air_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_air_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/airTemperature_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        lsScratchData_Imd.append(scratchGDB + "/airTemperature_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratchGDB + "/airTemperature_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/air_temperature_" + sTimeStamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="MEAN_air_temperature", out_ga_layer="#", \
            out_raster=outFolder + "/air_temperature_scratch_" + sTimeStamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="50", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")
        outExtractByMask = ExtractByMask(outFolder + '/air_temperature_scratch_' + sTimeStamp + '.tif', rc_elevation)
        outExtractByMask.save(outFolder + '/air_temperature_' + sTimeStamp + '.tif')
        lsScratchData_Imd.append(outFolder + "/air_temperature_scratch_" + sTimeStamp + ".tif")
    return outFolder + "/air_temperature_" + sTimeStamp + ".tif"
