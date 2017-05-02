# WCWAVE Climate Station Interpolation Toolkit (CSIT)

CSIT is a set of tools for creating spatially distributed climate data that can be used to force climate models. You can find more information about these tools at [the website][1]

## Installation


This repository includes two toolboxes that can be used in an ArcMap for Desktop installation (requires 10.3+). Just download the repository and extract it in a place that can be connected to ArcMap.  There are three tools included in the toolboxes. 

* Gridding tools Toolbox
  1. Climate Data Gridding Tools
  2. Reynolds Creek FTP to SQLite
* Cross Validation Toolkbox
  1. Cross Validation

The "Climate Data Gridding Tools" is the main tool used to create spatially distributed data. The "Reynolds Creek FTP to SQLite" tool can be used to convert [Reynolds Creek FTP Data][2] into the SQLite format required for the tool. "Cross Validation" is the tool used to run leave-one-out cross validation for multiple timesteps and can create simple xy graphs to compare the accuracy of the interpolation process.

The tool requires an elevation raster and a shapefile/featureclass of the station locations. To support Wind Speed the elevation raster needs to be in ASCII or geotiff format.  In some cases a soil station feature class with elevation in the attribute table. Optionally to support thermal radiation a View Factor Raster is required. View factor is calculated uses the techniques from  Dozier and Coutcalt (1979) equation 30 [Vf = cos^2(Horizon Angle)] (A horizon angle raster can be calculated with open source GIS tools. Such as GRASS GIS or Whitebox GIS).

### Server Installation

*Check Back Later*

### Description of Files

build_tables (dir) - old scripts for building the MySQL tables from ftp data
demo_data (dir) - Directory with sample data to test installations
joels_scripts (dir) - Origional scripts in separate files (used to create gridtools.py)
testing (dir) - directory used for testing scripts
utilites (dir) - ftp_2_sqlite.py script - used with ArcMap Toolbox
Cross_Validation.tbx - Cross validation toolbox for use in ArcMap
Gridding_Tools.tbx - Gridding Tools toolbox. Includes gridding tool and ftp2sqlite tool
cross_validation.py - script called from Cross Validation toolbox
gridtools.py - Gridding Tools script
gridtoolwrap.py - No longer supported. Used to create command line capable tool
runGriddingToolsJD.py - No longer supported. Gridding Tools on JD

## Naming convention

* initial:
   * z: elevation (m)
   * z_0: roughness_length (m)
   * z_s: snow_depth (m)
   * rho: snow_density (kg m^-3)
   * T_s_0: active_snow_layer_temperature (°C)
   * T_s: average_snow_covoer_temperature (°C)
   * h2o_sat: H2O_saturation (%)
* precip:
   * m_pp: precipitation_mass (kg m^-2 or mm)
   * %_snow: percent_snow (0 to 1.0)
   * rho_snow: precipitation_snow_density (kg m^-3)
   * T_pp: dew_point_temperature? (°C)
* input:
   * I_lw: thermal_radiation (W m^2)
   * T_a: air_temperature (°C)
   * e_a: vapor_pressure (Pa)
   * u: wind_speed (m s^-1)
   * T_g: soil_temperature (°C)
   * S_n: solar_radiation (W m^-2)
   
## References

Dozier, J., & Outcalt, S. I. (1979). An Approach toward Energy Balance Simulation over Rugged Terrain. Geographical Analysis, 11(1), 74. Retrieved from http://onlinelibrary.wiley.com/doi/10.1111/j.1538-4632.1979.tb00673.x/pdf

[1]: http://geoviz.geology.isu.edu/delparte_labs/VWCSIT/index.php
[2]: ftp://ftp.nwrc.ars.usda.gov/publicdatabase/reynolds-creek/
