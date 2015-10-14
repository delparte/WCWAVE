def precipMass():
    """ Kriging """
    arcpy.AddMessage("Calculating Precipitation Mass") 
    #Caclulate average precipitation values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/precipitation_table1", "precipitation_Table", "ppts > -500")
    arcpy.Statistics_analysis("precipitation_Table", "in_memory/precipitation_Table2", "ppts MEAN", "site_key")
    #Copy precipitation values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/precipitation_Table2", "site_key", "MEAN_ppts")
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for precipitation
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_ppts") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Add unique ID field to tempStations for use in OLS function
    arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")

    if sWatershed == "Johnston Draw":
        cursor = arcpy.SearchCursor(scratchGDB + "/tempStations")
        x = []
        y = []
        for row in cursor:
            x.append(row.getValue('RASTERVALU'))
            y.append(row.getValue('MEAN_ppts'))
        del cursor
        del row
        A = numpy.vstack([x,numpy.ones(len(x))]).T
        slope, intercept = numpy.linalg.lstsq(A, y)[0]
        arcpy.AddMessage('Slope ' + str(slope) + ', Intercept ' + str(intercept))
        if slope != 0.0 and intercept != 0.0:
            #Create final raster
            arcpy.env.extent = extElevation
            output_raster = (Raster(rc_elevation) * slope + intercept)
            output_raster.save(outFolder + "/precipitation_mass_" + sTimeStamp + ".tif")
        else:
            return
    else:
        #Check if interpolation method is Empirical Bayesian or Detrended
        if sKrigMethod == "Detrended":
            #Add unique ID field to tempStations for use in OLS function
            arcpy.AddField_management(scratchGDB +"/tempStations", "Unique_ID", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            arcpy.CalculateField_management(scratchGDB +"/tempStations","Unique_ID","!OBJECTID!","PYTHON_9.3","")
            #Run ordinary least squares on tempStations
            arcpy.CreateTable_management(scratchGDB, "coef_table")
            ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_ppts","RASTERVALU", scratchGDB + "/coef_table","","")
            lsScratchData_Imd.append(scratchGDB + "/coef_table")
            intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
            slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]
            #Calculate residuals and add them to tempStations
            arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
            for row in cursor:
                row.setValue("residual", row.getValue("MEAN_ppts") - ((slope * row.getValue("RASTERVALU")) + intercept))
                cursor.updateRow(row)
            del cursor
            del row
            #Run ordinary kriging on residuals
            #arcpy.gp.Kriging_sa(scratchGDB + "/tempStations", "residual", scratchGDB + "/precipitation_residual", "Linear 37.061494", output_cell_size, "VARIABLE 100", "")
            outKrig = Kriging(scratchGDB + "/tempStations", "residual", KrigingModelOrdinary("SPHERICAL", 460, 3686, .1214, .2192), output_cell_size, RadiusFixed(10000, 1))
            outKrig.save(scratchGDB + "/precipitation_residual")
            lsScratchData_Imd.append(scratchGDB + "/precipitation_residual")
            #Add back elevation trends and save final raster
            output_raster = arcpy.Raster(scratchGDB + "/precipitation_residual") + (arcpy.Raster(rc_elevation) * slope + intercept)
            output_raster.save(outFolder + "/precipitation_mass_" + sTimeStamp + ".tif")
        elif sKrigMethod == "Combined":
            #Run ordinary least squares on tempStations
            arcpy.CreateTable_management(scratchGDB, "coef_table")
            #Check if there was precipitation - I would like to change this check to occur before trying to run OLS
            try:
                ols = arcpy.OrdinaryLeastSquares_stats(scratchGDB + "/tempStations","Unique_ID","in_memory/fcResid","MEAN_ppts","RASTERVALU", scratchGDB + "/coef_table","","")
            except arcpy.ExecuteError:
                msgs = arcpy.GetMessages(2)
                #arcpy.AddMessage(msgs)
                if 'Zero variance' in msgs:
                    arcpy.AddMessage("No precipitation")
                else: 
                    arcpy.AddMessage(msgs)
                return
            lsScratchData_Imd.append(scratchGDB + "/coef_table")
            intercept = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[0]
            slope = list((row.getValue("Coef") for row in arcpy.SearchCursor(scratchGDB + "/coef_table",fields="Coef")))[1]

            #Calculate residuals and add them to tempStations
            arcpy.AddField_management(scratchGDB + "/tempStations", "residual", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
            cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
            for row in cursor:
                row.setValue("residual", row.getValue("MEAN_ppts") - ((slope * row.getValue("RASTERVALU")) + intercept))
                cursor.updateRow(row)
            del cursor
            del row
            arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="residual", out_ga_layer="#", \
                out_raster=scratchGDB + "/precipitation_residual", cell_size=output_cell_size, transformation_type="NONE",\
                max_local_points="100", overlap_factor="1", number_semivariograms="100", \
                search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
                output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", \
                semivariogram_model_type="THIN_PLATE_SPLINE")

            outExtractByMask = ExtractByMask(scratchGDB + '/precipitation_residual', rc_elevation)
            outExtractByMask.save(scratchGDB + '/precipitation_scratch')
            lsScratchData_Imd.append(scratchGDB + '/precipitation_residual')
            lsScratchData_Imd.append(scratchGDB + '/precipitation_scratch')
            #Add back elevation trends and save final raster
            output_raster = arcpy.Raster(scratchGDB + "/precipitation_scratch") + (arcpy.Raster(rc_elevation) * slope + intercept)
            output_raster.save(outFolder + "/precipitation_mass_" + sTimeStamp + ".tif")
        else:
            arcpy.EmpiricalBayesianKriging_ga(in_features=scratchGDB + "/tempStations", z_field="MEAN_ppts", out_ga_layer="#", \
                out_raster=outFolder + "/precipitation_mass_scratch_" + sTimeStamp + ".tif", cell_size=output_cell_size, \
                transformation_type="EMPIRICAL", max_local_points="100", overlap_factor="1", number_semivariograms="100", \
                search_neighborhood="NBRTYPE=SmoothCircular RADIUS=10000.9518700025 SMOOTH_FACTOR=0.2", \
                output_type="PREDICTION", quantile_value="0.5", threshold_type="EXCEED", probability_threshold="", \
                semivariogram_model_type="WHITTLE_DETRENDED")

            # Mask the output of EBK to size of input grid
            arcpy.AddMessage("Masking")
            outExtractByMask = ExtractByMask(outFolder + '/precipitation_mass_scratch_' + sTimeStamp + '.tif', rc_elevation)
            outExtractByMask.save(outFolder + '/precipitation_mass_' + sTimeStamp + '.tif')
            lsScratchData_Imd.append(outFolder + "/precipitation_mass_scratch_" + sTimeStamp + ".tif")
    
    return outFolder + "/precipitation_mass_" + sTimeStamp + ".tif"
