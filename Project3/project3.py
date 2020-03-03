# Install packages if missing
# !pip install netcdf4
# !pip install pydap
# !pip install wget
# !pip install metpy
# !pip install conda
# !conda install -c conda-forge metpy
# !apt-get -qq install libproj-dev proj-data proj-bin libgeos-dev
# !pip install Cython
# !pip install --upgrade --force-reinstall shapely --no-binary shapely
# !pip install cartopy

# Import libraries
import os
import xarray as xr
import wget
import glob
from bs4 import BeautifulSoup
import requests
import pandas as pd
import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import metpy.calc as mcalc
from metpy.units import units

%pylab inline 

import warnings
warnings.filterwarnings(action = "ignore", message = "^internal gelsd")

# Mount Drive
from google.colab import drive
drive.mount('/content/drive')

# Enter the directory in which you would be saving files and loading them
# We assume the aggregated data is stored in the datasets folder
PROJECT_ROOT_DIR = "./"
DATA_DIRECTORY = os.path.join(PROJECT_ROOT_DIR, "datasets/")
REFERENCE_PERIOD = "1997-2019"
REFERENCE_CITY = "Shanghai_JJA"
field_types = ["anomaly", "extreme", "ltm"]
IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "images/")
os.makedirs(IMAGES_PATH, exist_ok = True)

def load_data(variable_name, level, field_type):
    data_path = DATA_DIRECTORY + variable_name + '_' + level + "_" + field_type + '.nc'
    return xr.open_dataset(data_path)

def save_fig(fig_id, tight_layout=True, fig_extension="png", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)
    return None

# Part 1, aggregate daily rainfall data from GPCP
url = 'https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-daily/access/'
ext = 'nc'

def get_url_paths(url, ext='', params={}):
    response = requests.get(url, params=params)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]
    return parent

# Loop through all years and grab all of the datasets
years = np.arange(1996,2020)
datasets = []
for i in years:
    result = get_url_paths(url+'{i}/'.format(i=i),ext)
    print('working on {i} '.format(i=i))
    for j in range(len(result)):
        wget.download(result[j])
    files = glob.glob('gpcp*.nc')
    f = xr.open_mfdataset(files,concat_dim='time')
    var = xr.DataArray(f.precip.values,dims=['time','lat','lon'],
                    coords={'time':f.time.values,
                            'lat':f.latitude.values,
                            'lon':f.longitude.values})
    datasets.append(var)
    !rm gpcp*.nc
    #break
    
# Concatenate the datasets along time dimension
combined = xr.concat(datasets,dim='time')

# Convert to xarray dataset
combined_data = combined.to_dataset(name='precip')

# Convert to netCDF and save
combined_data.to_netcdf(DATA_DIRECTORY + '/GPCP_aggregate.nc', format='NETCDF4')

# Part 2, determine the 95% values of daily rainfall
# Create a function to find months JJA for Part 1/2
def is_jja(month):
    return (month >= 6) & (month <= 8)

# Open the combined dataset
combined_data = xr.open_dataset(DATA_DIRECTORY + 'GPCP_aggregate.nc')

# Slice data for June, July, and August only
jja_data = combined_data.sel(time=is_jja(combined_data['time.month']))

# Find data point: Shanghai, China lat, lon 31.2304° N, 121.4737° E
slat = 31.2304
slon = 121.4737

shanghai_jja = jja_data.sel(lon=slon, lat=slat, method="nearest")

# Find valid values (remove obvious incorrect values)
valid_ind = np.where((shanghai_jja.precip.values>=0.)&(shanghai_jja.precip.values<=1000.))

# Extract valid values
precip_shanghai = shanghai_jja.precip.values[valid_ind]

# Calculate 95 percentile
perc_95 = np.percentile(precip_shanghai,95)

# Plot CDF

# Plotting parameters
mpl.rcParams['xtick.major.size'] = 14
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 14
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 14
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 14
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.serif'] = ['Helvetica']
mpl.rc('xtick',labelsize=16)    # Formatting the x ticks
mpl.rc('ytick',labelsize=16)

# Plotting
counts, bin_edges = np.histogram (precip_shanghai, bins=50, density=True)
cdf = np.cumsum (counts)
plt.figure(figsize=(15,10))
plt.plot (bin_edges[1:], cdf/cdf[-1],'r',label='CDF')
plt.axvline(perc_95,c='b',ls='--',label='$95^{th}$ Percentile')
plt.xlabel('GPCP Daily Precipitation (mm)',fontsize=20)
plt.ylabel('Likelihood of Occurence',fontsize=20)
plt.legend(loc='upper left',fontsize=16)
fig_name = REFERENCE_CITY + "_" + REFERENCE_PERIOD + "_Daily_Precipitation(mm)_CDF&95_percentile"
save_fig(fig_name)
plt.show()

# Create Dataset for Shanghai significant points
shanghai_95th = shanghai_jja.where(shanghai_jja.precip>=perc_95,drop=True)

# Part 3, compute the global mean fields and JJA anomaly fields
# Season of JJA
JJA = [5, 6, 7]

# File strings for functions
filepath = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis'
fileperiod = '.1981-2010.ltm.nc'

# LONG TERM MEAN FUNCTION
def grab_LTM_data_sfc(variable):
    # Surface
    if variable == 'air':
        var = xr.open_dataset(filepath + '.derived/surface/air.sig995.mon' + str(fileperiod), engine='netcdf4').isel(nbnds=0)
    elif variable == 'uwnd':
        var = xr.open_dataset(filepath + '.derived/surface/uwnd.sig995.mon' + str(fileperiod), engine='netcdf4').isel(nbnds=0)
    elif variable == 'vwnd':
        var = xr.open_dataset(filepath + '.derived/surface/vwnd.sig995.mon' + str(fileperiod), engine='netcdf4').isel(nbnds=0)
    elif variable == 'pr_wtr':
        var = xr.open_dataset(filepath + '.derived/surface/pr_wtr.eatm.mon' + str(fileperiod), engine='netcdf4').isel(nbnds=0)
    else:
        print("Error. Please recheck argument inputs.")
    
    # Returns 2D Field
    return var

def grab_LTM_data(variable, lvl):
    var = xr.open_dataset(filepath + '.derived/pressure/' + str(variable) + '.mon'  + str(fileperiod), engine='netcdf4').isel(nbnds=0)
    x = var.sel(level=[lvl])
    
    # Returns 2D Field
    return x

# YEARLY FUNCTION
def grab_yearly_data_sfc(variable, year):
    # Surface
    if variable == 'air':
        var = xr.open_dataset(filepath + '.dailyavgs/surface/air.sig995.' + str(year) + '.nc', engine='netcdf4')
    elif variable == 'uwnd':
        var = xr.open_dataset(filepath + '.dailyavgs/surface/uwnd.sig995.' + str(year) + '.nc', engine='netcdf4')
    elif variable == 'vwnd':
        var = xr.open_dataset(filepath + '.dailyavgs/surface/vwnd.sig995.' + str(year) + '.nc', engine='netcdf4')
    elif variable == 'pr_wtr':
        var = xr.open_dataset(filepath + '.dailyavgs/surface/pr_wtr.eatm.' + str(year) + '.nc', engine='netcdf4')
    else:
        print("Error. Please recheck argument inputs.")
    
    # Returns 2D Field
    return var

def grab_yearly_data(variable, lvl, year):
    if variable == 'air':
        var = xr.open_dataset(filepath + '.dailyavgs/pressure/' + str(variable) + '.' + str(year) + '.nc', engine='netcdf4')
    elif variable == "shum":
        var = xr.open_dataset(filepath + '.dailyavgs/pressure/' + str(variable) + '.' + str(year) + '.nc', engine='netcdf4')
    else:
        var = xr.open_dataset(filepath + '.dailyavgs/pressure/' + str(variable) + '.' + str(year) + '.nc', engine='netcdf4')
    
    # Select level
    x = var.sel(level=[lvl])
    
    # Returns 2D Field
    return x

def calculate_means(variable, lvl, season):
    dataset = []
    dataset1 = []
    dataset2 = []
    t = pd.to_datetime(shanghai_95th.time.values)
    for i in range(len(shanghai_95th.time.values)):
        # Grab data for a specific year
        field = grab_yearly_data(variable, lvl, shanghai_95th.time.dt.year[i].values)

        # Select extreme precip day
        # Special cases: air and shum
        if variable == 'air':
            extreme = field.sel(time = field.time.dt.dayofyear == t[i].timetuple().tm_yday).squeeze()
            extreme['air'] = extreme['air'] - 273.15 # Change to Celsius from Kelvin
            extreme.assign_attrs({'units':'degC'})
        elif variable == 'shum':
            extreme = field.sel(time = field.time.dt.dayofyear == t[i].timetuple().tm_yday).squeeze()
            extreme['shum'] = extreme['shum']/1000 # Change to kg/kg from g/kg
            extreme.assign_attrs({'units':'kg/kg'})
        else:
            extreme = field.sel(time = field.time.dt.dayofyear == t[i].timetuple().tm_yday).squeeze()

        # Select seasonal long term mean
        ltm = grab_LTM_data(variable, lvl).isel(time=season).mean(dim='time').squeeze()

        # Calculate anomaly
        anomaly = extreme[variable] - ltm[variable]

        # Append
        dataset.append(anomaly)
        dataset1.append(extreme[variable])
        dataset2.append(ltm[variable])

    anomalies = xr.concat(dataset, dim='index').mean(dim="index")
    extremes = xr.concat(dataset1, dim='index').mean(dim="index")
    ltms = xr.concat(dataset2, dim='index').mean(dim="index")
    return anomalies, extremes, ltms

def calculate_means_sfc(variable, season):
    dataset = []
    dataset1 = []
    dataset2 = []
    t = pd.to_datetime(shanghai_95th.time.values)
    for i in range(len(shanghai_95th.time.values)):
        # Grab data for a specific year
        field = grab_yearly_data_sfc(variable, shanghai_95th.time.dt.year[i].values)

        # Select extreme precip day
        # Special case: air
        if variable == 'air':
            extreme = field.sel(time = field.time.dt.dayofyear == t[i].timetuple().tm_yday).squeeze()
            extreme['air'] = extreme['air'] - 273.15 # Change to Celsius from Kelvin
            extreme.assign_attrs({'units':'degC'})
        else:
            extreme = field.sel(time = field.time.dt.dayofyear == t[i].timetuple().tm_yday).squeeze()

        # Select seasonal long term mean
        ltm = grab_LTM_data_sfc(variable).isel(time=season).mean(dim='time').squeeze()

        # Calculate anomaly
        anomaly = extreme[variable] - ltm[variable]

        # Append
        dataset.append(anomaly)
        dataset1.append(extreme[variable])
        dataset2.append(ltm[variable])
          
    anomalies = xr.concat(dataset, dim='index').mean(dim="index")
    extremes = xr.concat(dataset1, dim='index').mean(dim="index")
    ltms = xr.concat(dataset2, dim='index').mean(dim="index")
    return anomalies, extremes, ltms

###### Upper levels
level = 850
variable = 'vwnd'
temp = calculate_means(variable, level, JJA)
print('finished calculating means')
which_mean = ['anomaly', 'extreme', 'ltm']
for a in range(len(temp)):
    combined_data = temp[a].to_dataset(name=variable)
    combined_data.to_netcdf(DATA_DIRECTORY + str(variable) + '_' + str(level) + '_' + str(which_mean[a]) + '.nc', format='NETCDF4')

###### Surface
variable = 'pr_wtr'
temp = calculate_means_sfc(variable, JJA)
print('finished calculating means')
which_mean = ['anomaly', 'extreme', 'ltm']
for a in range(len(temp)):
    combined_data = temp[a].to_dataset(name=variable)
    combined_data.to_netcdf(DATA_DIRECTORY + str(variable) + '_sfc_' + str(which_mean[a]) + '.nc', format='NETCDF4')

# Calculate wind speed
def calculate_wspd(level="sfc", field_type="ltm"):
    # For this, we need to calculate from the UWND and VWND components.
    uwnd = load_data("uwnd", level=level, field_type=field_type)
    vwnd = load_data("vwnd", level=level, field_type=field_type)
    lat, lon = uwnd["lat"].values, uwnd["lon"].values
    speed = mcalc.wind_speed(uwnd['uwnd'].values * units('m/s'),
                             vwnd['vwnd'].values * units('m/s')).m
    data = xr.DataArray(speed, coords=[lat, lon], dims=['lat', 'lon']).to_dataset(name = 'wspd')
    data.to_netcdf(DATA_DIRECTORY + "wspd" + '_' + str(level) + '_' + field_type + '.nc', 
                   format='NETCDF4')
    return None

levels = ["250", "850"]
for field_type in field_types:
    for level in levels:
        calculate_wspd(level=level, field_type=field_type)
        
# Part 4, plotting
def wind_vector_plot(uwnd=None, vwnd=None, level="sfc", field_type="ltm", 
                     projection="PlateCarree", plot_type="streamplot", 
                     figsize=(20, 20), plt_show=True):
    """
    Adopted from ATMS 597 M03 N02 jupyter notebok
    The crs will be a rotated
    pole CRS, meaning that the vectors will be unevenly spaced in
    regular PlateCarree space.
    """
    plt.figure(figsize=figsize)
    
    if not uwnd and not vwnd:
        uwnd = load_data("uwnd", level=level, field_type=field_type)
        vwnd = load_data("vwnd", level=level, field_type=field_type)
    
    x = uwnd["lon"].values - 180.
    y = uwnd["lat"].values
    u = uwnd["uwnd"].values
    v = vwnd["vwnd"].values
    
    if projection == "PlateCarree":
        ax = plt.axes(projection=ccrs.PlateCarree())
        transform = ccrs.PlateCarree()
        if plot_type == "streamplot":
            ax.streamplot(x, y, u, v, transform=transform)
        elif plot_type == "quiver":
            ax.quiver(x, y, u, v, transform=transform)
        elif plot_type == "barbs":
            ax.barbs(x, y, u, v, transform=transform)
        else:
            print("Unavailable plot type, please select from 'streamplot', 'quiver', or 'barbs'")
            return None
    elif projection == "Orthographic":
        crs = ccrs.RotatedPole(pole_longitude=slon, pole_latitude=slat)
        transform = ccrs.Orthographic(slon, slat)
        ax = plt.axes(projection=transform)
        lat = range(5, 35)
        lon = range(35, 65)
        if plot_type == "streamplot":
            ax.streamplot(x[lon], y[lat], u[lat, lon], v[lat, lon], transform=crs)
        elif plot_type == "quiver":
            ax.quiver(x[20:80], y[0:40], u[0:40,20:80], v[0:40,20:80], transform=crs)
        elif plot_type == "barbs":
            ax.barbs(x[40:60], y[10:30], u[10:30,40:60], v[10:30,40:60], transform=crs)
        else:
            print("Unavailable plot type, please select from 'streamplot', 'quiver', or 'barbs'")
            return None
    else:
        print("Unavailable projection method, please select from 'PlateCarree', or 'Orthographic'")
        return 
    
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND, edgecolor="black")
#     ax.set_global()
#     ax.gridlines()
    fig_name = "Wind_Vector_" + level + "_" + field_type + "_" + projection + "_" + plot_type
    if field_type == "ltm":
        fig_name = fig_name + "_" + REFERENCE_PERIOD
    else:
        fig_name = fig_name + "_" + REFERENCE_CITY
    plt.title(fig_name)
    save_fig(fig_name)
    if plt_show:
        plt.show()
    else:
        plt.close()
    return None


def contour_plot(variable_name, data=None, level="sfc", field_type="ltm", 
                 projection="PlateCarree", figsize=(13, 7), plt_show=True):
    """
    Adopted from ATMS 597 M03 N02 jupyter notebok
    The crs will be a rotated
    pole CRS, meaning that the vectors will be unevenly spaced in
    regular PlateCarree space.
    """
    plt.figure(figsize=figsize)
    
    if variable_name == "wspd":
        unit = "m/s"
    elif variable_name == "hgt":
        unit = "m"
    elif variable_name == "air":
        unit = "degreeC"
    elif variable_name == "shum":
        unit = "kg/kg"
    elif variable_name == "pr_wtr":
        unit = "kg/m^2"
    else:
        print("Unavialble variable")
        return None
    
    if not data:
        data = load_data(variable_name, level=level, field_type=field_type)
    
    if projection == "PlateCarree":
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.))
        v_max, v_min = max(map(max, data[variable_name])), min(map(min, data[variable_name]))
        data[variable_name].plot.contourf(ax=ax, # cmap="seismic",
                                          transform=ccrs.PlateCarree(),
                                          levels=10,
                                          cbar_kwargs={"label": unit}) 
                                          #, levels=np.linspace(v_min, v_max, 25))

    else:
        print("Unavailable projection method, please select from 'PlateCarree', or 'Orthographic'")
        return 

#     ax.set_global()
    ax.coastlines()
#     ax.gridlines()
    fig_name = "Contour_" + variable_name + "_" + level + "_" + field_type + "_" + projection
    if field_type == "ltm":
        fig_name = fig_name + "_" + REFERENCE_PERIOD
    else:
        fig_name = fig_name + "_" + REFERENCE_CITY
    plt.title(fig_name)
    save_fig(fig_name)
    if plt_show:
        plt.show()
    else:
        plt.close()
    return None

for field_type in field_types:
    wind_vector_plot(level="250", field_type=field_type, plt_show=False)
    contour_plot("wspd", level="250", field_type=field_type, plt_show=False)

for field_type in field_types:
    wind_vector_plot(level="500", field_type=field_type, plt_show=False)
    contour_plot("hgt", level="500", field_type=field_type, plt_show=False)

for field_type in field_types:
    wind_vector_plot(level="850", field_type=field_type, plt_show=False)
    contour_plot("wspd", level="850", field_type=field_type, plt_show=False)
    contour_plot("air", level="850", field_type=field_type, plt_show=False)
    contour_plot("shum", level="850", field_type=field_type, plt_show=False)

for field_type in field_types:
    wind_vector_plot(level="sfc", field_type=field_type, plt_show=False)
    contour_plot("air", level="sfc", field_type=field_type, plt_show=False)

for field_type in field_types:
    contour_plot("pr_wtr", level="sfc", field_type=field_type, plt_show=False)
