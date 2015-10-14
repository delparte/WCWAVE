def snowProperties():
    arcpy.AddMessage("Calculating Snow Properties") 
    if len(density_interp_values['features']) <= 1:
        #Density Equation: y = -0.0395(elevation) + 405.26
        snow_density_raster = -0.0395 * Raster(rc_elevation) + 405.26
        snow_density_raster.save(outFolder + "/snow_density_" + sTimeStamp + ".tif")
    else:
        lsElevation = []
        lsDensity = []
        for rec in density_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsDensity.append(rec['attributes']['Density'])

        lr_results = stats.linregress(lsElevation,lsDensity)
        slope1 = lr_results[0]
        intercept1 = lr_results[1]
        snow_density_raster = slope1 * Raster(rc_elevation) + intercept1
        snow_density_raster.save(outFolder + "/snow_density_" + sTimeStamp + ".tif")

    if len(ul_interp_values['features']) <= 1:
        #Upper Layer Temperature Equation: y = -0.0008(elevation) + 0.1053
        upper_layer_temperature = -0.0008 * Raster(rc_elevation) + 0.1053
        upper_layer_temperature.save(outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif")
    else:
        lsElevation = []
        lsTemperature = []
        for rec in ul_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsTemperature.append(rec['attributes']['Temperature'])

        lr_results = stats.linregress(lsElevation,lsTemperature)
        slope2 = lr_results[0]
        intercept2 = lr_results[1]
        upper_layer_temperature = slope2 * Raster(rc_elevation) + intercept2
        upper_layer_temperature.save(outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif")

    if len(ll_interp_values['features']) <= 1:
        #lower layer temperature equation: y = -0.0008(elevation) + 1.3056
        lower_layer_temperature = -0.0008 * Raster(rc_elevation) + 1.3056
    else:
        lsElevation = []
        lsTemperature = []
        for rec in ul_interp_values['features']:
            lsElevation.append(rec['attributes']['Elevation'])
            lsTemperature.append(rec['attributes']['Temperature'])

        lr_results = stats.linregress(lsElevation,lsTemperature)
        slope3 = lr_results[0]
        intercept3 = lr_results[1]
        lower_layer_temperature = slope2 * Raster(rc_elevation) + intercept2

    #average snowcover temperature is the average of the upper and lower layer temperatures
    average_snowcover_temperature = arcpy.sa.CellStatistics([upper_layer_temperature, lower_layer_temperature],"MEAN","NODATA")
    average_snowcover_temperature.save(outFolder + "/average_snow_cover_temperature_" + sTimeStamp + ".tif")

    return outFolder + "/active_snow_layer_temperature_" + sTimeStamp + ".tif", outFolder + "/average_snow_cover_temperature_" + sTimeStamp + ".tif", outFolder + "/snow_density_" + sTimeStamp + ".tif"


