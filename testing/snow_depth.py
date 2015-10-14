def snowDepth():
    """ Kriging """
    arcpy.AddMessage("Calculating Snow Depth") 
    #Caclulate average snow depth values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/snowDepth_table1", "snowDepth_Table", "snow_depth > -500")
    arcpy.Statistics_analysis("snowDepth_Table", "in_memory/snowDepth_Table2", "snow_depth MEAN", "site_key")
    #Copy snow depth values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/snowDepth_Table2", "site_key", "MEAN_snow_depth")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for snow depth
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    values = list()
    for row in cursor:
        if row.getValue("MEAN_snow_depth") is None:
            cursor.deleteRow(row)
        else:
            values.append(row.getValue('MEAN_snow_depth'))
    del cursor
    del row
    average = numpy.mean(values) 
    count = int(arcpy.GetCount_management(scratchGDB + '/tempStations').getOutput(0))
    if count >= 10 and average > 0: 
        #Check if interpolation method is Empirical Bayesian or Detrended
        if sKrigMethod == "Detrended":
            #Add unique ID field to tempStations for use in OLS function
            arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
            #Run ordinary least squares on tempStations
            arcpy.CreateTable_management(scratchGDB, "coef_table")
            ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_snow_depth","RASTERVALU", scratchGDB + "/coef_table","","")
            lsScratchData_Imd.append(scratchGDB + "/coef_table")
            intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
            slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
            #Calculate residuals and add them to tempStations
            arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
            for row in cursor:
                row.setValue("residual", row.getValue("MEAN_snow_depth") - ((slope * row.getValue("RASTERVALU")) + intercept))
                cursor.updateRow(row)
            del cursor
            del row
            #Run ordinary kriging on residuals
            #arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/snowDepth_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
            outKrig = Kriging(scratchGDB + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
            outKrig.save(scratchGDB + "/snowDepth_residual")
            lsScratchData_Imd.append(scratchGDB + "/snowDepth_residual")
            #Add back elevation trends and save final raster
            imdRaster = arcpy.Raster(scratchGDB + "/snowDepth_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
            output_raster = Con(imdRaster > 0, imdRaster, 0)
            output_raster.save(outFolder + "/snow_depth_" + sTimeStamp + ".tif")
        else:
            #Add unique ID field to tempStations for use in OLS function
            arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
            #Run ordinary least squares on tempStations
            arcpy.CreateTable_management(scratchGDB, "coef_table")
            ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_snow_depth","RASTERVALU", scratchGDB + "/coef_table","","")
            lsScratchData_Imd.append(scratchGDB + "/coef_table")
            intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
            slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
            #Calculate residuals and add them to tempStations
            arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
            for row in cursor:
                row.setValue("residual", row.getValue("MEAN_snow_depth") - ((slope * row.getValue("RASTERVALU")) + intercept))
                cursor.updateRow(row)
            del cursor
            del row
            arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
                out_raster=scratchGDB + "/snowDepth_residual", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
                overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
                output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
            #Add back elevation trends and save final raster
            imdRaster = arcpy.Raster(scratchGDB + "/snowDepth_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
            output_raster = Con(imdRaster > 0, imdRaster, 0)
            output_raster.save(outFolder + "/snow_depth_" + sTimeStamp + ".tif")
    else: 
        if count < 10:
            arcpy.AddMessage('Not enough data for snow depth.  Try a different time step')
        if average == 0:
            arcpy.AddMessage('No snow on the ground. Try a different time step if needed.')

    return outFolder + "/snow_depth_" + sTimeStamp + ".tif"


