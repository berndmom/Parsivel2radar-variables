# -*- coding: utf-8 -*-
# Code written by: Bernd Mom, email: bernd.mom@helsinki.fi

# Meaning of the keys (skipping not important variables)
# [date:yyyy-mm-dd HH:MM:SS, 01:Rain intensity (mm/h), 02:Rain amount accumulated (mm),
# 03:Weather code according to SYNOP wawa, 07:Radar reflectivity (dBZ), 08:MOR visibility in the precipitation (m),
# 11:Number of detected particles, 13:Serial number of instrument, 18:Sensor status,
# 34:Kinetic energy (KJ), 90:Field N (d) log(1/m3mm), 91:Field v(d) (m/s), 93:Raw data,]

import glob
import numpy as np
import xarray as xr

#from parsivel import parsivel
from DSD_parameters import DSD_parameters
from parsivel2netcdf import parsivel2netcdf
from filter_parsivel import parsivel
from apply_pytmatrix import radar_parameters

for month in range(1,13):
# 1) .txt to .nc conversion
	# Format: monthly -> ['year', 'month']; day -> ['year', 'month', 'day']; period -> ['year', 'month', 'day start', 'day end'].
	date = [2023, month]	
	
	# Example file is given in folder
	path_in = "Insert path to location of the Parsivel text file"
	path_out = "Insert folder in which to save the file"

	to_netcdf = parsivel2netcdf(date, path_in, path_out)
	data_raw = to_netcdf()

# 2) filter raw data
	path_in = path_out + "{:04d}".format(date[0]) + "{:02d}".format(date[1]) + "parsi23.nc"
	path_out = "Insert path to save file"

	run_parsivel = parsivel(path_in, path_out)
	run_parsivel.time_interval = 1
	run_parsivel.length = 5
	run_parsivel.DSD_parameters = DSD_parameters()
	data_filtered = run_parsivel.load_data()

# 3) calculate radar variables
	# Radar variables calculated for vertical geometry   (reflectivity, mean Doppler velocity, specific attenuation).
	# Radar variables calculated for horizontal geometry (reflectivity, mean Doppler velocity, specific differential phase, differential reflectivity).
	# Below an example of horizontal geometry.
	# Look-up tables added for C-band horizontal and vertical, e.g. "ltable_Choriz_8mm.txt" and "ltable_Cvert_8mm.txt" respectively.
	
	radar_scat = radar_parameters("ltable_Choriz_8mm.txt", np.arange(0, data_filtered.dims["time"], 1), data_filtered.D0.values, data_filtered.mu.values, data_filtered.Nw.values)
	radar_scat.wavelength = "C" 		# Wavelength type
	radar_scat.geometry = "horizontal" 	# horizontal or vertical
	radar_scat.limit = [.4,-1,100]		# Lower limits for the DSD-parameters [D0, mu, Nw].
	data_scatter_vert = radar_scat.prep_scatter().reset_index("index", drop=True).assign_coords(time=data_filtered.time.values).rename({"index":"time"})

# End

