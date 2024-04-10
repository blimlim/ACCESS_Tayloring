# Script to adjust format of diddled files to match that of the original ESM1.5 AMIP ancillaries.
#  Adjusts dimensions, names, units, attributes, masks data, adds a 
# singleton "surface" dimenstion and saves data using encodings from an xconved 
# netcdf of the original AMIP data.

# Reads input file and output path from command line
# Input file should have the formatting of the output from the diddling code. I.e. it should 
# have fields names "SST_cpl" and "ice_cov" for the sst's and ice concentrations respectively.

# Script should be run with two command line arguments ifile and ofile, the first pointing to the input file. 
# The script will produce two output files, "ofile_sst.nc" and "ofile_ice.nc".
#-------------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import xarray as xr
import sys
args = sys.argv[1:]
if len(args) != 2:
    raise Exception("Usage: python bcgen_to_ESM_format.py ifile ofile")

ifile = args[0]
ofile = args[1]

# Land sea mask for the ESM grid. Ocean = 1, land = 0
maskfile = "/home/565/sw6175/prepping_forcings/regrid_attempts/iris_regridding/mask_grids/Target_mask_grid.nc"

# Load the mask
lsm = xr.load_dataset(maskfile)


# Store the attributes and encodings for the original ESM AMIP ancillary files, so that we dont
# need to load in the files here

ESM_info = {
    "longitude": 
    {
        "attrs": {'long_name': 'longitude','standard_name': 'longitude','units': 'degrees_east','point_spacing': 'even','modulo': ' '},
        "encoding": {'dtype': np.dtype('float32')}
    },
    "latitude": 
    {
        "attrs": {'long_name': 'latitude','standard_name': 'latitude','units': 'degrees_north','point_spacing': 'even'},
        "encoding": {'dtype': np.dtype('float32')}
    },
    "surface": 
    {
        "attrs": {'long_name': 'Surface', 'units': 'level', 'positive': 'up'},
        "encoding": {'dtype': np.dtype('float32')}
    },
    "t": 
    {
        "attrs": {'long_name': 't', 'time_origin': '01-JAN-1870:00:00:00'},
        "encoding": {'dtype': np.dtype('float32'),'units': 'days since 1870-01-01 00:00:00', 'calendar': 'gregorian'}
    },
    "temp": 
    {
        "attrs": {'source': 'No field processing','name': 'temp','title': 'SURFACE TEMPERATURE AFTER TIMESTEP','date': '01/01/70','time': '00:00','long_name': 'SURFACE TEMPERATURE AFTER TIMESTEP','standard_name': 'surface_temperature','units': 'K'},
        "encoding": {'dtype': np.dtype('float32'),'missing_value': 2e+20,'_FillValue': 2e+20}
    },
    "iceconc": 
    {
        "attrs": {'source': 'No field processing','name': 'iceconc','title': 'FRAC OF SEA ICE IN SEA AFTER TSTEP', 'date': '01/01/70', 'time': '00:00', 'long_name': 'FRAC OF SEA ICE IN SEA AFTER TSTEP', 'standard_name': 'sea_ice_area_fraction', 'units': '0-1'},
        "encoding": {'dtype': np.dtype('float32'),'missing_value': 2e+20,'_FillValue': 2e+20}
    }   
}





# Functions to do each modification
#-------------------------------------------------------------------------------------------------------

# Apply the mask
# Assumes the mask dataset has a field called "mask" and that the longitudes and latitudes are correct
# and that 1 = ocean, 0 = land
def mask(dataset, mask):
    dataset["mask"] = mask.mask
    dataset["SST_cpl"] = dataset.SST_cpl.where(dataset.mask == 1, np.nan)
    dataset["ice_cov"] = dataset.ice_cov.where(dataset.mask == 1, np.nan)
    
    # Drop the mask from the dataset 
    dataset = dataset.drop("mask")
    return dataset   


# Convert SST's from deg C to K
# ice doesn't need a unit conversion
def convert_units(dataset):
    dataset["SST_cpl"] = dataset["SST_cpl"] + 273.15
    dataset["SST_cpl"].attrs["units"] = "K"
    return dataset

# Rename all of the dimensions and variables to match the ESM ancillaries
def rename_vars(dataset):
    dataset = dataset.rename({"time": "t", "lon": "longitude", "lat": "latitude"})
    dataset = dataset.rename({"SST_cpl": "temp"})
    dataset = dataset.rename({"ice_cov": "iceconc"})
    return dataset
    


# Function to set the encodings and attributes
def set_enc_attr(dataset, infosource):
    # Set all the attributes and encodings to match the ESM ones
    # It is assumed that dataset has variables and dimensions renamed already to match the ESM ones
    
    for coord in dataset.coords:
        dataset[coord].attrs = infosource[coord]["attrs"]
        dataset[coord].encoding = infosource[coord]["encoding"]
        # The ESM amip file has no fill value for the coordinates
        dataset[coord].encoding["_FillValue"] =  None

    for var in dataset.keys():
        dataset[var].attrs = infosource[var]["attrs"]
        dataset[var].encoding = infosource[var]["encoding"]
        
    # set t to an unlimited dimension
    dataset.encoding["unlimited_dims"] = "t"
    dataset.attrs["source"] = "Regridded and diddled SST and ice data, names and attributes modified for use as an ESM1.5 AMIP ancillary"
    return dataset




# Function to do all the formatting
def to_esm(dataset):
    # Make a deep copy of the dataset to work on
    # If we don't, changes get written to the dataset that is inputed 
    # I don't quite understand why, but that's what it is doing
    # I think it must be something to do xarray datasets being mutable.
    # This is probably memory inefficient to be doing this though (don't really know what I'm saying)
    working_set = dataset.copy(deep = True)
    # We are free to modify working_set from here
    working_set = working_set.drop_vars(["date", "datesec"])
    if "SST_cpl_prediddle" in working_set.keys():
        working_set = working_set.drop_vars("SST_cpl_prediddle")
    if "ice_cov_prediddle" in working_set.keys():
        working_set = working_set.drop_vars("ice_cov_prediddle")

    # Check that the dataset contains the correct fields, error out if not
    if not("SST_cpl" in working_set.keys()) or not("ice_cov" in working_set.keys()):
        raise Exception("Missing SST_cpl or ice_cov variables, does your variable have the correct name?")
        
    
    # Mask the data 
    print("Masking data")
    working_set = mask(working_set, lsm)
    # Rename the variables and convert units
    working_set = convert_units(working_set)
    # Get it to use "gregorian" rather than "proleptic-gregorian" calendar
    # Probably doesn't matter too much
    working_set = working_set.convert_calendar("gregorian", use_cftime = True)

    working_set = rename_vars(working_set)
    
    
    # Add a surface dimension with a single level
    working_set = working_set.expand_dims(dim = {"surface": np.array([0.0], dtype = "float32")})
    
    # Reorder the dimensions so that they match the ESM ancillary files
    working_set = working_set.transpose("t", "surface", "latitude","longitude",)
    
    # Set the attributes and encodings to match the ESM ancillary files
    working_set = set_enc_attr(working_set, ESM_info)

    return working_set


# Run and save
#==========================================================================================
unformatted_ds = xr.load_dataset(ifile) 
print("Formatting dataset")
formatted_ds = to_esm(unformatted_ds)

formatted_ds.temp.to_netcdf(ofile + "_sst.nc")
formatted_ds.iceconc.to_netcdf(ofile + "_ice.nc")
print("Done!")