def thermalRadiation(inAirTemperature, inVaporPressure, inSurfaceTemp, v_factor):
    arcpy.AddMessage("Calculating Thermal Radiation") 
    #Constants and re-defined variables (See Marks and Dozier (1979), pg. 160)
    z = rc_elevation
    vf = v_factor
    T_a = inAirTemperature
    vp = inVaporPressure

    #Select one station and its data to use as reference for elevation, pressure, and temperature
    #Caclulate average air temperature and vapor pressure values over the n-hour time step (ignore "no-data" values: -999)
    arcpy.MakeTableView_management(scratchGDB + "/climate_table1", "thermalRadiation_Table", "air_temperature > -500 AND vapor_pressure > -500")
    arcpy.Statistics_analysis("thermalRadiation_Table", "in_memory/thermalRadiation_Table2", [["air_temperature", "MEAN"], ["vapor_pressure", "MEAN"]], "site_key")
    #Copy vapor pressure values to station feature class
    arcpy.CopyFeatures_management(fcStations_wElevation, scratchGDB +"/tempStations")
    lsScratchData_Imd.append(scratchGDB +"/tempStations")
    arcpy.JoinField_management(scratchGDB +"/tempStations", "Site_Key", "in_memory/thermalRadiation_Table2", "site_key", ["MEAN_air_temperature", "MEAN_vapor_pressure"])
    #Delete rows from station feature class that have negative or null elevations
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("RASTERVALU") < 0 or row.getValue("RASTERVALU") == "None" or row.getValue("RASTERVALU") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Delete rows from station feature class that have null values for vapor pressure or air temperature
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        if row.getValue("MEAN_air_temperature") is None or row.getValue("MEAN_vapor_pressure") is None:
            cursor.deleteRow(row)
    del cursor
    del row
    #Use update cursor to extract air temperature, vapor pressure, and elevation values for the first station, then delete that station
    P_m = 0.0
    T_m = 0.0
    z_m = 0.0
    T_s = inSurfaceTemp
    cursor = arcpy.UpdateCursor(scratchGDB + "/tempStations")
    for row in cursor:
        z_m = row.getValue("RASTERVALU")
        P_m = row.getValue("MEAN_vapor_pressure")
        T_m = row.getValue("MEAN_air_temperature")
        cursor.deleteRow(row)
        break

    del cursor
    del row

    arcpy.AddMessage("P_m: " + str(P_m))
    arcpy.AddMessage("T_m: " + str(T_m))
    arcpy.AddMessage("z_m: " + str(z_m))
    arcpy.AddMessage("T_s: " + str(T_s))

    g = 9.8
    m = 0.0289
    R = 8.3143
    sigma = 5.6697*10**-8
    epsilon_s = 0.95
    gamma = -0.006

    #convert temperature parameters to Kelvin
    T_m = T_m + 274.15
    T_s = T_s + 274.15
    T_a = arcpy.sa.Float(Raster(T_a) + 274.15)

    #convert vapor pressure parameters to mb
    P_m = P_m * 0.01
    vp = arcpy.sa.Float(Raster(vp) * 0.01)

    #Correct air temperature and vapor pressure rasters (Marks and Dozier (1979), pg. 164)
    T_prime = T_a + (0.0065 * Raster(rc_elevation)) #(4) corrected air temperature
    e_sa = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from original air temperature (T_a)
    e_sprime = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from corrected air temperature (T_prime)
    rh = arcpy.sa.Float(vp / e_sa) #(5) relative humidity
    e_prime = arcpy.sa.Float(rh * e_sprime) #(6) corrected vapor pressure

    #Pressure at a given elevation (Marks and Dozier (1979), pg. 168-169)
    term1 = ((-g*m)/(R*gamma))
    delta_z = Raster(z) - z_m
    term2 = ((T_m + gamma * delta_z)) / T_m
    lnTerm = arcpy.sa.Ln(term2)
    expTerm = arcpy.sa.Exp(term1 * lnTerm)
    P_a = P_m * expTerm #(10) air pressure

    #effective emissivity (Marks and Dozier (1979), pg. 164)
    epsilon_a = arcpy.sa.Float((1.24 * (e_prime / T_prime)**(1/7)) * (P_a / 1013.0)) #(7)

    #Incoming longwave radiation (Marks and Dozier (1979), pg. 164)
    term3 = arcpy.sa.Float((epsilon_a * sigma * (T_a ** 4)) * vf)
    term4 = arcpy.sa.Float(epsilon_s * sigma * (T_s ** 4))
    term5 = (1 - Raster(vf))
    output_thermal_radiation = arcpy.sa.Float(term3 + (term4 * term5)) #(9)
    output_thermal_radiation.save(outFolder + "/thermal_radiation_" + sTimeStamp + ".tif")

    return outFolder + "/thermal_radiation_" + sTimeStamp + ".tif"


