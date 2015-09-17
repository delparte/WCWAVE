#Import necessary modules
import arcpy
from arcpy.sa import *
import numpy

#Check-out necessary extensions
arcpy.CheckOutExtension('Spatial')

#Set input parameters
elevation_raster = arcpy.GetParameterAsText(0)
view_factor_raster = arcpy.GetParameterAsText(1)
air_temperature_raster = arcpy.GetParameterAsText(2)
vapor_pressure_raster = arcpy.GetParameterAsText(3)
reference_air_pressure = arcpy.GetParameter(4)
reference_air_temperature = arcpy.GetParameter(5)
reference_elevation = arcpy.GetParameter(6)
surface_air_temperature = arcpy.GetParameter(7)
outRaster = arcpy.GetParameterAsText(8)

#Set up workspace
scratchWS = arcpy.env.scratchWorkspace
scratchGDB = arcpy.env.scratchGDB
#output cell size and processing extent should be the same as elevation raster
arcpy.env.cellSize = elevation_raster
output_cell_size = arcpy.env.cellSize
arcpy.env.extent = elevation_raster
arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "75%"
arcpy.Delete_management("in_memory")

#Constants and re-defined variables (See Marks and Dozier (1979), pg. 160)
z = elevation_raster
vf = view_factor_raster
T_a = air_temperature_raster
vp = vapor_pressure_raster
P_m = reference_air_pressure
T_m = reference_air_temperature
z_m = reference_elevation
T_s = surface_air_temperature
g = 9.8
m = 0.0289
R = 8.3143
sigma = 5.6697*10**-8
epsilon_s = 0.95
gamma = -0.006

#convert temperature parameters to Kelvin
arcpy.AddMessage("Converting temperature parameters to Kelvin")
T_m = T_m + 274.15
T_s = T_s + 274.15
T_a = arcpy.sa.Float(Raster(T_a) + 274.15)

#Correct air temperature and vapor pressure rasters (Marks and Dozier (1979), pg. 164)
arcpy.AddMessage("Correcting air temperature and vapor pressure rasters")
T_prime = T_a + (0.0065 * Raster(elevation_raster)) #(4) corrected air temperature
e_sa = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from original air temperature (T_a)
e_sprime = arcpy.sa.Float(6.11 * 10**((7.5*arcpy.sa.Float(T_a))/(237.3 + arcpy.sa.Float(T_a)))) #saturated vapor pressure from corrected air temperature (T_prime)
rh = arcpy.sa.Float(vp / e_sa) #(5) relative humidity
e_prime = arcpy.sa.Float(rh * e_sprime) #(6) corrected vapor pressure

#Pressure at a given elevation (Marks and Dozier (1979), pg. 168-169)
arcpy.AddMessage("Calculating air pressures at given elevations")
term1 = ((-g*m)/(R*gamma))
delta_z = Raster(z) - z_m
term2 = ((T_m + gamma * delta_z)) / T_m
lnTerm = arcpy.sa.Ln(term2)
expTerm = arcpy.sa.Exp(term1 * lnTerm)
P_a = P_m * expTerm #(10) air pressure

#effective emissivity (Marks and Dozier (1979), pg. 164)
arcpy.AddMessage("Calculating effective emissivity")
epsilon_a = arcpy.sa.Float((1.24 * (e_prime / T_prime)**(1/7)) * (P_a / 1013.0)) #(7)

#Incoming longwave radiation (Marks and Dozier (1979), pg. 164)
arcpy.AddMessage("Calculating incoming longwave radiation")
term3 = arcpy.sa.Float((epsilon_a * sigma * (T_a ** 4)) * vf)
term4 = arcpy.sa.Float(epsilon_s * sigma * (T_s ** 4))
term5 = (1 - Raster(vf))
output_thermal_radiation = arcpy.sa.Float(term3 + (term4 * term5)) #(9)
arcpy.AddMessage("Creating final raster")
output_thermal_radiation.save(outRaster)
