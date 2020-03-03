# ATMS 597 Spring 2020 - Project 3

### Group C: Piyush Garg, Yang Lu, Sarah Szymborski

Create code using python `xarray` to organize and reduce climate data and plot with `cartopy`

#### Data source:
- Global Precipitaiton Climatology Project (GPCP) - Daily data, from https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-daily/access/
- NCEP Reanalysis, from https://journals.ametsoc.org/doi/pdf/10.1175/1520-0477(1996)077%3C0437%3ATNYRP%3E2.0.CO%3B2

#### Required inputs:
- Attributes, including wind speed, geopotential height, temperature, specific humidity, atmospheric column water vapor
- Levels, including surface, 250 hPa, 500 hPa, 850 hPa
- Global long term mean in JJA 1981-2010, extreams in Shanghai JJA 1981-2010, and global anomaly in JJA 1981-2010
- Plotting options

#### Outputs:
- Aggregated daily rainfall data from GPCP from 1981 to 2010
- From NCEP,
  - 250 hPa wind vectors and wind speed, 
  - 500 hPa wind vectors and geopotential height,
  - 850 hPa temperature, specific humidity, and winds,
  - skin temperature, and surface wind vectors (sig995 level), and
  - total atmospheric column water vapor.  

#### Caveats:
- The stream plots in Orthographic projection doesn't fit well

------
Last updated: 03/03/2020
