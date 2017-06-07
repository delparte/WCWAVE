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

**Publish CSIT to geoprocessing server**

This guide assumes you have a working ArcGIS server instance and you have Publisher privileges (ArcGIS 10.4+ restricts publisher privileges to Administrator Accounts only.  Instructions to give publishers those privileges can be found [on the Esri resources website](http://server.arcgis.com/en/server/latest/administer/windows/change-geoprocessing-service-and-service-extension-publishing-privileges.htm)).

Geoprocessing services cannot have external executables in them (ie. WindNinja)  There is currently a workaround where these executables can be commented out and then uncommented after they are published.

Comment out `ninja_path` variable (~line 988) and the `args` array (Between comments ~line 1022-1046)

The service needs to be run successfully before it can be uploaded. Run the tool with minimal settings (Watershed: Reynolds Creek; From Date: 2008-01-01 12:00:00; To Date: 2008-01-01 13:00:00 Air Temperature: True).  After the tool has completed successfully open up the Results pane (ArcMap Menu: Geoprocessing > Results).  Under the "Current Session" dropdown right click on the "Climate Data Gridding Tools" and click "Share As > Geoprocessing Service"

Choose "Publish a service" (or "Overwrite an existing service").

[Screenshot]

Click the button to create a new connection, or select the connection you want to use.

_Create a new connection_ 

- Choose "Publish GIS services"
- Fill in the Server URL and Authentication information and click Finish
- On the Publish a Service Dialogue choose the new connection 

Change the Service Name and click next.

[Screenshot - obscure ip address]

Choose the folder to Publish the service to (Create a new one or use an existing one) and click "Continue".

_Service Editor Dialogue_

General - Check the information is correct on the general tab. 

[Screenshot - obscure ip]

Capabilities - Leave parameters as default.

Parameters - Change Message Level to "Info" (CSIT uses `arcpy.AddMessage()` to relay script data)

Pooling - Adjust Timeouts to your needs.

Processes - Leave as default.

Climate Data Gridding Tools - Adjust Descriptions as needed.

Item Description - Summary, Tags, and Descriptions are not added by default but are required.  Add whatever is not present.

Sharing - Change as need.

_Publishing the Service_

Double check everything using the "Analyze" Button.  Fix any errors. (A few warnings will show, Code 24046 warning can be ignored; 24032 is a warning that data needs to be uploaded and that can also be ignored.)

When all warnings and errors are resolved Click Publish. (If it fails to publish check the results window to see any errors or problems).


### Description of Files

- build_tables (dir) - old scripts for building the MySQL tables from ftp data
- demo_data (dir) - Directory with sample data to test installations
- joels_scripts (dir) - Origional scripts in separate files (used to create gridtools.py)
- testing (dir) - directory used for testing scripts
- utilites (dir) - ftp_2_sqlite.py script - used with Gridding Tools Toolbox
- Cross_Validation.tbx - Cross validation toolbox for use in ArcMap
- Gridding_Tools.tbx - Gridding Tools toolbox. Includes gridding tool and ftp2sqlite tool
- cross_validation.py - script called from Cross Validation toolbox
- gridtools.py - Gridding Tools script
- gridtoolwrap.py - No longer supported. Used to create command line capable tool
- runGriddingToolsJD.py - No longer supported. Gridding Tools on Johnston Draw (use gridtools.py now)

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
