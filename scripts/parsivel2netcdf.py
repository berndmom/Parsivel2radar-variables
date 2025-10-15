# -*- coding: utf-8 -*-
# Code written by: Bernd Mom, email: bernd.mom@helsinki.fi

import glob
import pandas as pd
import numpy as np
import xarray as xr
import datetime as dt

from netCDF4 import Dataset    # Note: python is case-sensitive!
from netCDF4 import date2num,num2date

class parsivel2netcdf(object):
    """parsivel2netcdf
       Reading Parsivel text file and converts it to netcdf.        
       Only the variables "synop", "number_particles", "number_density", "particle_speed" and "raw_data" are saved.
       No more variables needed to calculate radar variables. If necessary to store more variables, add own code.
    
    Attributes:
    	path_in:	Path to Parsivel folder or file (raw data).
      	path_out:	Path to save Parsivel raw data in netcdf form.
      	date:		The date for which the data needs to be read. Format: monthly -> ['year', 'month']; 
      			day -> ['year', 'month', 'day']; period -> ['year', 'month', 'day start', 'day end'].
    """
        
    def __init__(self, date, path_in=" ", path_out=" "):
        self.path = path_in
        self.save_file = path_out
        if self.path==" " or self.save_file==" ":
            print("Please insert correct path for input or output file")
        self.date = date
        
    def __call__(self):
        self.read_parsivel()
        self.write_to_netcdf()
        return self.data
    
    def read_parsivel(self):
        date_month = self.prep_data() 
        
        self.store_files = []
        
        # store the file in list (every list is a seperate time step)
        for file in self.files:
            store_file = []
            with open(file) as infile:  
                for line in infile:          
                    # get the necessary information for each time step          
                    if line[0:3] in ['[20']:
                        file_part = []
                        file_part = file_part + [line[1:-1]] 
                    elif line[0:2] in ['03','11','90','91','93']:
                        file_part = file_part + [line[:-1]]
                    elif line[-2] == ']':
                        store_file += [file_part]
            self.store_files += [store_file]
        
        count_files = 0      # Keeps track of the files present in directory)
    
        for time in self.day_range:
            print(time)
           
            if count_files < len(self.files):     # Checks whether correspoding file number is accurate
                date_check = date_month[date_month.day==time]
                #print(date_check)
            
                # Checks whether the date corresponds with the date in the file name
                if str(self.date[0]) + "{:02d}".format(self.date[1]) + "{:02d}".format(time) in self.files[count_files]:
                    self.handle_parsivel(count_files, date_check)
                    count_files +=1
        
    def prep_data(self):
        date_range = pd.date_range(start=str(self.date[0]) + "-01-01", end=str(self.date[0]+1) + "-01-01", freq="1min")[:-1]  # Creates an array of a year
        date_month = date_range[date_range.month==self.date[1]]
        
        self.path = self.path + str(self.date[0]) + "{:02d}".format(self.date[1]) + "*.txt"
        all_files = glob.glob(self.path)
        
        if len(self.date)==3:
            date_days = date_month[date_month.day==self.date[2]]
            self.day_range = list(set(date_days.day))
            dates = [dt.datetime(self.date[0],self.date[1],self.date[2])+n*dt.timedelta(minutes=1) for n in range(date_days.shape[0])]
            self.create_dict(dates)
            
            for file in all_files: 
                if str(self.date[0]) + "{:02d}".format(self.date[1]) + "{:02d}".format(self.date[2]) in file:
                    self.files = [file]     
        elif len(self.date)==4:
            date_days = date_month[(date_month.day>=self.date[2]) & (date_month.day<=self.date[3])]
            self.day_range = list(set(date_days.day))
            dates = [dt.datetime(self.date[0],self.date[1],self.date[2])+n*dt.timedelta(minutes=1) for n in range(date_days.shape[0])]
            self.create_dict(dates)
            
            self.files = []
            for i in range(self.date[2],self.date[3]+1):
                for file in all_files:
                    if str(self.date[0]) + "{:02d}".format(self.date[1]) + "{:02d}".format(i) in file:
                        self.files += [file]
            
            self.files = sorted(self.files) 
        else:
            self.day_range = list(set(date_month.day)) 
            dates = [dt.datetime(self.date[0],self.date[1],1)+n*dt.timedelta(minutes=1) for n in range(date_month.shape[0])]
            self.create_dict(dates)  
            self.files = all_files
            self.files = sorted(self.files)  
        return date_month
        
    def handle_parsivel(self, count_files, date_check_all):        
        dates = []
        count_minutes = 0
        count_time = 0

        # Retrieve data
        for date_check in date_check_all:
            count_total = (count_files) * 1440 + count_minutes
            
            if count_time < len(self.store_files[count_files]):
                date = pd.to_datetime(self.store_files[count_files][count_time][0], format='%Y-%m-%d %H:%M:%S')
            
                # Check for double times
                if date_check > date:
                    while date in dates:
                        self.data["status"]["data"][count_time] = 0
                        count_time +=1
                        date = pd.to_datetime(self.store_files[count_files][count_time][0], format='%Y-%m-%d %H:%M:%S')
            
                if date not in dates:
                    dates = dates + [date_check]

                    if date_check == date:
                        for part in self.store_files[count_files][count_time]:
                            self.data["status"]["data"][count_total] = 1
                            if part[0:2] == '03': 
                                self.data["synop"]["data"][count_total] = float(part[3:])
                            elif part[0:2] == '11':
                                self.data["number_particles"]["data"][count_total] = float(part[3:])
                            elif part[0:2] == '90':
                                temp_part = part[3:-1].split(';')
                                if len(temp_part) == 32:
                                    temp = np.float_(temp_part)
                                    self.data["number_density"]["data"][count_total,:] = np.where(temp==-9.999, 0.0, temp) 
                                else:
                                    print('Corrupted data: ', date )                                 
                            elif part[0:2] == '91':
                                temp_part = part[3:-1].split(';')
                                if len(temp_part) == 32:
                            	        self.data["particle_speed"]["data"][count_total,:] = np.float_(temp_part) 
                            elif part[0:2]=='93':
                                temp = np.float_(part[3:-1].split(';'))
                                if temp.shape[0] == 1024:
                                    # read the PSD and velocity fields
                                    for qq in range(0,32): # velocity class
                                        for nn in range(0,32): # diameter class
                                            self.data["raw_data"]["data"][count_total,qq,nn] = temp[(qq*32 + nn)]
            
                        count_time +=1   
                        count_minutes +=1
                    else:
                        self.data["status"]["data"][count_total] = 0
                        count_minutes +=1
        
    def write_to_netcdf(self):
        # Opening a file, creating a new dataset
        try: ncfile.close()  # just to be safe, make sure dataset is not already open.
        except: pass

        file_name = self.save_file + str(self.date[0]) + "{:02d}".format(self.date[1]) + 'parsi23.nc'
        ncfile = Dataset(file_name,mode='w', format='NETCDF4_CLASSIC') 
        
        # Creating dimensions
        time_dim    = ncfile.createDimension('time', None) # unlimited axis (can be appended to).
        class_D_dim = ncfile.createDimension('class_diameter', 32)
        class_V_dim = ncfile.createDimension('class_velocity', 32)
        
        # Creating attributes
        ncfile.title = 'Hyytiälä Parsivel23 data'
        
        # Creating variables
        for var_name, var_attrs in self.data.items():
            print(var_name)
            if var_name == 'time':
                var = ncfile.createVariable(var_name, "f4", ("time",)) # Set dimensions
                var[:] = var_attrs["data"] # Set data
            elif var_name == 'class_diameter' or var_name == 'class_velocity':
                var = ncfile.createVariable(var_name, "f4", (var_name))
                var[:] = var_attrs["data"]
            elif var_name == 'number_density':
                var = ncfile.createVariable(var_name, "f4", ("time", "class_diameter"))
                var[:] = var_attrs["data"]
            elif var_name == 'particle_speed':
                var = ncfile.createVariable(var_name, "f4", ("time", "class_velocity"))
                var[:] = var_attrs["data"]
            elif var_name == 'raw_data':
                var = ncfile.createVariable(var_name, "f4", ("time", "class_diameter", "class_velocity"))
                var[:] = var_attrs["data"]
            else:
                var = ncfile.createVariable(var_name, "f4", ("time",))
                var[:] = var_attrs["data"]
              
            # Setting attributes 
            for attr_name, attr_value in var_attrs.items():
                if attr_name != "data":
                    setattr(var, attr_name, attr_value)
                    
        ncfile.close()

    def create_dict(self, days):
        # Create time range for netcdf
        since = "hours since " + str(self.date[0]) + "-" + "{:02d}".format(self.date[1]) + "-01"
        times = date2num(days, since, has_year_zero=self.date[0], calendar='standard')
        
        self.data = {}
        
        variables = ["time",     "class_diameter",          "class_velocity",          "status",                     "synop",                    "number_particles",             "number_density", 
                     "particle_speed",             "raw_data"]
        data_type = [times,      np.arange(1, 32+1),        np.arange(1, 32+1),        np.empty(times.shape[0]),     np.empty(times.shape[0]),   np.empty(times.shape[0]),       np.empty((times.shape[0],32)), 
                    np.empty((times.shape[0],32)), np.empty((times.shape[0],32,32))]
        units     = [since,      '-',                       '-',                       '-',                          '-',                        "#",                            "log(1/m3mm)",
                     "m/s",                        "#"]
        long_name = ["Time UTC", "Diameter classes (1-32)", "Velocity classes (1-32)", "Status 0: no data; 1: data", "Synop code based on WaWa", "Number of detected particles", "Field N(D)",
                     "Field v(D)",                 "Number of detected particles raw data (32x32 classes)"]
        
        for i in range(len(variables)):
            if i == 0:  
                self.data[variables[i]] = {
                    "data": data_type[i],
                    "units": units[i],
                    "long_name": long_name[i],
                    "calendar": "standard"
                }
            else:
                if i > 2:
                    data_type[i][:] = np.nan
                
                self.data[variables[i]] = {
                    "data": data_type[i],
                    "units": units[i],
                    "long_name": long_name[i],
                }

#%% End
