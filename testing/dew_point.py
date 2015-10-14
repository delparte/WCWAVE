def dewPoint():
    """ Kriging """
    arcpy.AddMessage("Calculating Dew Point") 
    #Caclulate average dew point temperature values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "dewPoint_Table", "dewpoint_temperature > -500")
    arcpy.Statistics_analysis("dewPoint_Table", "in_memory/dewPoint_Table2", "dewpoint_temperature MEAN", "site_key")
    #Copy dew point temperature values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/dewPoint_Table2", "site_key", "MEAN_dewpoint_temperature")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for dew point temperatures
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_dewpoint_temperature") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    if sKrigMethod == "Combined":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_dewpoint_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_dewpoint_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
            out_raster=scratchGDB + "/dewPoint_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
        lsScratchData_Imd.append(scratchGDB + "/dewPoint_residual")
        
        outExtractByMask = ExtractByMask(scratchGDB + '/dewPoint_residual', rc_elevation)
        outExtractByMask.save(scratchGDB + '/dewPoint_scratch')
        lsScratchData_Imd.append(scratchGDB + '/dewPoint_scratch')
        
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratchGDB + "/dewPoint_scratch") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif")
    #Check if interpolation method is Empirical Bayesian or Detrended
    elif sKrigMethod == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_dewpoint_temperature","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_dewpoint_temperature") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/dewPoint_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratchGDB + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratchGDB + "/dewPoint_residual")
        lsScratchData_Imd.append(scratchGDB + "/dewPoint_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratchGDB + "/dewPoint_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="MEAN_dewpoint_temperature", out_ga_layer="#", \
            out_raster=outFolder + "/dew_point_temperature_scratch_" + sTimeStamp + ".tif", cell_size=output_cell_size, transformation_type="EMPIRICAL", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")
        outExtractByMask = ExtractByMask(outFolder + '/dew_point_temperature_scratch_' + sTimeStamp + '.tif', rc_elevation)
        outExtractByMask.save(outFolder + '/dew_point_temperature_' + sTimeStamp + '.tif')
        lsScratchData_Imd.append(outFolder + "/dew_point_temperature_scratch_" + sTimeStamp + ".tif")
    
    #Create precipitation properties rasters (% precip as snow, density of snow portion of precip)
    #Percent snow conditional:
    inRas = Raster(outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif")
    outPercentSnowRaster = arcpy.sa.Con(inRas < -5.0, 1.0, Con((inRas >= -5.0) & (inRas < -3.0), 1.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 1.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 1.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0.0), 0.75,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 0.25,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outPercentSnowRaster, outFolder + "/percent_snow_" + sTimeStamp + ".tif",pixel_type="32_BIT_FLOAT")

    #Snow density conditional:
    outSnowDensityRaster = arcpy.sa.Con(inRas < -5.0, 75.0, Con((inRas >= -5.0) & (inRas < -3.0), 100.0,
                                                               Con((inRas >= -3.0) & (inRas < -1.5), 150.0,
                                                                   Con((inRas >= -1.5) & (inRas < -0.5), 175.0,
                                                                       Con((inRas >= -0.5) & (inRas < 0), 200.0,
                                                                           Con((inRas >= 0.0) & (inRas < 0.5), 250.0,
                                                                               Con(inRas >= 0.5,0.0)))))))
    arcpy.CopyRaster_management(outSnowDensityRaster,outFolder + "/precipitation_snow_density_" + sTimeStamp + ".tif",pixel_type="32_BIT_FLOAT")


    return outFolder + "/dew_point_temperature_" + sTimeStamp + ".tif", outFolder + "/percent_snow_" + sTimeStamp + ".tif", outFolder + "/precipitation_snow_density_" + sTimeStamp + ".tif"


