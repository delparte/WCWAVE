# WCWAVE
Scripts for iSNOBAL Model

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
