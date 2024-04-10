# taylor_test.py

# Usage: taylor_test.py taylored_file.nc plot_directory

# Script to test the SST and sea ice output of the bcgen code, based on the Met Office's
# tayloring_simulator code. Goal is to see how well the monthly means 
# calculated after the model's linear interpolation to daily values matches the original 
# supplied monthly means.

# Performs the following steps:
# 1. Convert SSTs from C to K
# 2. Interpolate SSTs and sea ice to daily values
# 3. Clip temp to be above 271.36 K, ice to be between 0 and 1
# 4. Calculate monthly means from the daily timeseries
# 5. Cut off first/last months - I'm unsure what the models do in the
# first and final month, and the errors here tend to dominate
# 6. Compute max monthly error between the recovered means and the means in the 
# original data
# 7. Plot the max monthly errors. 


import xarray as xr 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import sys


# Contour levels for plotting
levels_sst = [0, 5e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 1e-1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0]
levels_ice = [5e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3,0.5]

def taylor_test(taylor_output):
    """
    Simulate the linear interpolation of monthly means to daily values of SSTs and sea ice performed by 
    atmosphere models, and then recalculate monthly means from the interpolated data.
    
    Parameters:
        taylor_output (xr dataset):  the output from the bcgen Tayloring/diddling software. Contains diddled fields 
        (SST_cpl,, ice_cov), and copies of the original undiddled supplied SST and sea ice data
        (SST_cpl_prediddle, ice_cov_prediddle), with dimensions time, lat, and lon
    Returns: 
        recovered_means (xr dataset): Monthly means recalculated from the daily interpolated data. Process is applied to 
        both the diddled (i.e. Taylor processed) and raw (original supplied) data, in order to demonstrate the reduction in
        error provided by the diddling procedure.
    """
    # Step 1: Convert temperatures to K
    # Be careful here... the taylor_output ds is mutable and so this is modifying the original data.
    # That's ok for me for calculating the errors later but I should clean this up a bit.
    taylor_output["SST_cpl"] = taylor_output.SST_cpl + 273.15
    taylor_output.SST_cpl.attrs["units"] = "K"
    taylor_output["SST_cpl_prediddle"] = taylor_output.SST_cpl_prediddle + 273.15
    taylor_output.SST_cpl_prediddle.attrs["units"] = "K"

    # Step 2: Interpolate SSTs and sea ice to daily values
    # Sampling at midday appears to be important, otherwise recovered error is larger
    days = np.arange(taylor_output.time.values[0], taylor_output.time.values[-1], dtype = "datetime64[D]") +np.timedelta64(12, 'h')

    # Interpolate all variables within the dataset
    full_interp = taylor_output.copy().interp(time = days, method = "linear")
    
    # Step 3: Clip the data 
    # Met office version uses 271.36, but as far as I can tell, the UM (at least vn7.3) clips to 273.15
    full_interp["SST_cpl"] = full_interp.SST_cpl.clip(min = 271.35)
    full_interp["SST_cpl_prediddle"] = full_interp.SST_cpl_prediddle.clip(min = 271.35) 
    full_interp["ice_cov"] = full_interp.ice_cov.clip(min = 0.0, max = 1.0)
    full_interp["ice_cov_prediddle"] = full_interp.ice_cov_prediddle.clip(min = 0.0, max = 1.0)

    # Step 4: Recalculate monthly means from the interpolated and clipped data
    # First group into months, calculate mean, and recombine

    # index containing month and year to bin the data
    # Method from https://stackoverflow.com/a/58980577
    year_month_idx = pd.MultiIndex.from_arrays([full_interp['time.year'].values, full_interp['time.month'].values]).values
    full_interp.coords['year_month'] = ('time', year_month_idx)

    # Calculate mean for each month
    recovered_means = full_interp.groupby("year_month").mean(dim = "time")

    # Recovert the time data from the original dataset
    recovered_means = recovered_means.rename({"year_month":"time"})
    recovered_means["time"] = taylor_output.time

    return recovered_means 

def max_recovered_error(taylor_output, recovered_mean_ds):
    """
    Calculate the maximum absolute error in the recovered monthly means compared to the 
    original supplied data.

    Note that the first and final months are dropped, as I'm unsure how the model deals with them
    and the recovered errors there tend to dominate.

        Parameters:
            taylor_output (xr dataset):  the output from the bcgen Tayloring/diddling software. Contains diddled fields 
            (SST_cpl,, ice_cov), and copies of the original undiddled supplied SST and sea ice data
            (SST_cpl_prediddle, ice_cov_prediddle)
    
            recovered_mean_ds (xr dataset): output from taylor_test(taylor_output)

        Returns: 
            max_error_SST_taylored (xr dataarray): Maximum monthly error in the taylored SSTs compared to the original SSTs
            max_error_SST_raw (xr dataarray): Maximum monthly error in SSTs had the tayloring procedure not been applied
            max_error_ice_taylored (xr dataarray): Maximum monthly error in the taylored sea ice compared to the original SSTs
            max_error_ice_raw (xr dataarray): Maximum monthly error in the sea ice had the tayloring procedure not been applied
    """
    error_block = xr.Dataset(
        coords = dict(
            lat = (["lat"], taylor_output.lat.values),
            lon = (["lon"], taylor_output.lon.values)
        )
    )
    error_block = error_block.assign(
        dict(
            max_error_SST_taylored = np.abs((recovered_mean_ds.SST_cpl - taylor_output.SST_cpl_prediddle)[1:-1]).max(dim = "time"),
            max_error_SST_raw = np.abs((recovered_mean_ds.SST_cpl_prediddle - taylor_output.SST_cpl_prediddle)[1:-1]).max(dim = "time"),
            max_error_ice_taylored = np.abs((recovered_mean_ds.ice_cov - taylor_output.ice_cov_prediddle)[1:-1]).max(dim = "time"),
            max_error_ice_raw = np.abs((recovered_mean_ds.ice_cov_prediddle - taylor_output.ice_cov_prediddle)[1:-1]).max(dim = "time")
        )
    )


    return error_block

def plot_error(error_ds, var, levels, cbar_label, filepath):
    """
    Saves contour plot of maximum monthly error between processed and original monthly means
    """
    fig, ax = plt.subplots(figsize = (10,6),subplot_kw = dict(projection = ccrs.PlateCarree(central_longitude = 200)))
    ax.coastlines()


    # Add cyclic point to data to remove missing vertical line in contour plot
    data = error_ds[var]
    
    lon = data.coords['lon']
    lon_idx = data.dims.index('lon')
    wrap_data, wrap_lon = add_cyclic_point(data.values, coord=lon, axis=lon_idx)
    
    cyclic_data = xr.DataArray(
        wrap_data,
        coords = dict(
            lat = data.coords['lat'],
            lon = wrap_lon
        )
    )


    p = cyclic_data.plot.contourf(
        add_colorbar = False,
        levels = levels,
        transform = ccrs.PlateCarree()
    )
    plt.colorbar(p, label = cbar_label, fraction = 0.025)
    plt.title(var.replace("_"," "))
    
    plt.savefig(filepath + "/" + var + ".png")




if __name__ == "__main__":
    args = sys.argv[1:]
    
    ifile = args[0]
    plotdir = args[1]

    print("loading data")
    taylor_output = xr.open_dataset(ifile, chunks = {"lat":2,"lon":2})
    print("Interpolating to daily values and recalculating monthly means")
    recovered_means = taylor_test(taylor_output).compute()
    print("Calculating max monthly errors")
    max_errors = max_recovered_error(taylor_output.compute(), recovered_means)
    print("Saving figures")
    plot_error(max_errors, "max_error_SST_taylored", levels_sst, "Max monthly error (K)", plotdir)
    plot_error(max_errors, "max_error_SST_raw", levels_sst, "Max monthly error (K)", plotdir)
    plot_error(max_errors, "max_error_ice_taylored", levels_ice, "Max monthly error (fraction)", plotdir)
    plot_error(max_errors, "max_error_ice_raw", levels_ice, "Max monthly error (fraction)", plotdir)
    # max_errors.to_netcdf("taylor_test.nc")
    print("Done")







