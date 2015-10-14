def vaporPressure():
    """ Kriging """
    arcpy.AddMessage("Calculating Vapor Pressure") 
    #Caclulate average vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "vaporPressure_Table", "vapor_pressure > -500")
    arcpy.Statistics_analysis("vaporPressure_Table", "in_memory/vaporPressure_Table2", "vapor_pressure MEAN", "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/vaporPressure_Table2", "site_key", "MEAN_vapor_pressure")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    
    if sKrigMethod == "Combined":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_vapor_pressure","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_vapor_pressure") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
            out_raster=scratchGDB + "/vaporPressure_scratch", cell_size=output_cell_size, transformation_type="NONE", max_local_points="100", \
            overlap_factor="1", number_semivariograms="100", search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
            output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="THIN_PLATE_SPLINE")
        lsScratchData_Imd.append(scratchGDB + "/vaporPressure_scratch")

        # Mask the output of EBK to size of input grid
        arcpy.AddMessage("Masking")
        outExtractByMask = ExtractByMask(scratchGDB + '/vaporPressure_scratch', rc_elevation)
        outExtractByMask.save(scratchGDB + '/vaporPressure_residual')

        lsScratchData_Imd.append(scratchGDB + "/vaporPressure_residual")
        
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratchGDB + "/vaporPressure_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/vapor_pressure_" + sTimeStamp + ".tif")
    elif sKrigMethod == "Detrended":
        #Add unique ID field to tempStations for use in OLS function
        arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
        #Run ordinary least squares on tempStations
        arcpy.CreateTable_management(scratchGDB, "coef_table")
        ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_vapor_pressure","RASTERVALU", scratchGDB + "/coef_table","","")
        lsScratchData_Imd.append(scratchGDB + "/coef_table")
        intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
        slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
        #Calculate residuals and add them to tempStations
        arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
        for row in cursor:
            row.setValue("residual", row.getValue("MEAN_vapor_pressure") - ((slope * row.getValue("RASTERVALU")) + intercept))
            cursor.updateRow(row)
        del cursor
        del row
        #Run ordinary kriging on residuals
        #arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/vaporPressure_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
        outKrig = Kriging(scratchGDB + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
        outKrig.save(scratchGDB + "/vaporPressure_residual")
        lsScratchData_Imd.append(scratchGDB + "/vaporPressure_residual")
        #Add back elevation trends and save final raster
        output_raster = arcpy.Raster(scratchGDB + "/vaporPressure_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
        output_raster.save(outFolder + "/vapor_pressure_" + sTimeStamp + ".tif")
    else:
        arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="MEAN_vapor_pressure",\
                out_ga_layer="#", out_raster=outFolder + "/vapor_pressure_scratch_" + sTimeStamp + ".tif", cell_size=output_cell_size,\
                transformation_type="EMPIRICAL", max_local_points="100", overlap_factor="1", number_semivariograms="100",\
                search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", output_type="PREDICTION",\
                quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", semivariogram_model_type="WHITTLE_DETRENDED")

        # Mask the output of EBK to size of input grid
        arcpy.AddMessage("Masking")
        outExtractByMask = ExtractByMask(outFolder + '/vapor_pressure_scratch_' + sTimeStamp + '.tif', rc_elevation)
        outExtractByMask.save(outFolder + '/vapor_pressure_' + sTimeStamp + '.tif')
        lsScratchData_Imd.append(outFolder + "/vapor_pressure_scratch_" + sTimeStamp + ".tif")
        
    return outFolder + "/vapor_pressure_" + sTimeStamp + ".tif"


