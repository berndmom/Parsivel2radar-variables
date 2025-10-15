# -*- coding: utf-8 -*-
# Code written by: Bernd Mom, email: bernd.mom@helsinki.fi

import numpy as np
import xarray as xr

class parsivel(object):
    """Load and filter Parsivel data
    
    Class for simulating the drop size distributions of rain observed by a disdr
    
    Attributes:
        path_in:	Path to Parsivel folder or file (raw data) in netcdf.
        path_out:	Path to save file in netcdf form.
        time_interval:	Interval time (min.) to create larger interval. Parsivel interval time is 1 min.
        		This can increased by setting a value in minutes.        
        limit_synop:	Filter the data based on synop code. Default: [50, 53, 60, 63] (e.g. rain). 
        		Format: [limit 1 start, limit 1 end, limit 2 start, limit 2 end]
        limit_par_num:	Filter the date based on the number of particles observed. Default: 100 droplets.
        border:		Take into account the droplet falling at the edges of the laser beam. Default: False
        area:		Area of the laser beam. Default: [30*10**-3, 180*10**-3] (e.g. Parsivel) 
        dsd_length:	Length of consecutive bins in which droplets needs to be observed.
        gap:		The minimum gap between the consecutive bins and seperate bin or between two seperate bins,
        		to be excluded from the observations.  
        DSD_parameters: Link to file where the drop size distrubution parameters are calculated.       		
    """
    
    def __init__(self, path_in, path_out, synop=[50, 53, 60, 63], par_num=100, border=False, area=[30*10**-3, 180*10**-3], length=5, gap=3):
        self.path_in = path_in
        self.path_out = path_out
        self.time_interval = 1    # minutes    
        
        self.limit_synop = synop
        self.limit_par_num = par_num
        self.border = border
        self.area = area # 30 mm, 180 mm
        self.dsd_length = length
        self.dsd_gap = gap
        
        self.DSD_parameters = None
        
    def load_data(self):
        self.data = xr.open_dataset(self.path_in)       
        self.data["time"] = self.data["time"].dt.round("S")
        
        if self.time_interval != 1:
            self.time_integration()
        else:
            self.temp_data = self.data
   
        if self.DSD_parameters is None:
            return print('Set distributions for D0, mu and Nw')
        elif (self.DSD_parameters is not None):          
            self.filter_all()
            self.temp_data = self.DSD_parameters(self.temp_data)
        else:
            return print("Check if input variables are correct")    
        
        self.temp_data.to_netcdf(self.path_out)
        return self.temp_data  
    
    def time_integration(self): 
        data_status = self.data["status"].resample(time=str(self.time_interval)+"min").max()
        data_synop = self.data["synop"].resample(time=str(self.time_interval)+"min").max()
        data_num_par = self.data["number_particles"].resample(time=str(self.time_interval)+"min").sum()
        data_num_den = self.data["number_density"].resample(time=str(self.time_interval)+"min").mean()
        data_par_spd = self.data["particle_speed"].resample(time=str(self.time_interval)+"min").mean()     # nan or 0.0
        data_raw = self.data["raw_data"].resample(time=str(self.time_interval)+"min").sum()
        self.temp_data = xr.merge([data_status, data_synop, data_num_par, data_num_den, data_par_spd, data_raw])
        return print("Finished time integration")

    def filter_all(self):
        # Filter data by synop code and number of detected particles
        self.temp_data = self.temp_data.where((self.temp_data["synop"]>=self.limit_synop[0]) & (self.temp_data["synop"]<=self.limit_synop[1]) | 
                                              (self.temp_data["synop"]>=self.limit_synop[2]) & (self.temp_data["synop"]<=self.limit_synop[3]), np.nan)
        self.temp_data = self.temp_data.where((self.temp_data["number_particles"]>=self.limit_par_num), np.nan) 
        # Filter data by velocity-diameter mask
        self.filter_VD()
        # Filter for number of consecutive DSD and gaps
        self.filter_gaps()
    
    def filter_gaps(self):
        for i in range(0, len(self.temp_data.time)):
            if np.nansum(10**self.temp_data["number_density_filtered"][i])>0:
                data_where = np.where(10**self.temp_data["number_density_filtered"][i]>0)[0]
                
                for j in range(data_where.shape[0]):   
                    if data_where.shape[0] < self.dsd_length:
                        self.temp_data["number_density_filtered"][i][:] = np.nan
                        self.temp_data["particle_speed"][i][:] = np.nan
                    else:                  
                        if j < 1 and (data_where[j+1] - data_where[j]) > self.dsd_gap:
                            self.temp_data["number_density_filtered"][i][data_where[j]] = np.nan
                            self.temp_data["particle_speed"][i][data_where[j]] = np.nan
                        elif j > data_where.shape[0]-2 and (data_where[j] - data_where[j-1]) > self.dsd_gap:
                            self.temp_data["number_density_filtered"][i][data_where[j]] = np.nan
                            self.temp_data["particle_speed"][i][data_where[j]] = np.nan   
                        elif j > 0 and j < data_where.shape[0]-1:
                            if (data_where[j+1] - data_where[j]) > self.dsd_gap and (data_where[j] - data_where[j-1]) > self.dsd_gap:
                                self.temp_data["number_density_filtered"][i][data_where[j]] = np.nan
                                self.temp_data["particle_speed"][i][data_where[j]] = np.nan                        
                                
                data_where = np.where(10**self.temp_data["number_density_filtered"][i]>0)[0]
                indx = []
        
                for j in range(data_where.shape[0]):
                    if j == 0:
                        if (data_where[j+1] - data_where[j])<2:
                            indx += [data_where[j]]
                    elif j == data_where.shape[0]-1: 
                        if (data_where[j] - data_where[j-1])<2:
                            indx += [data_where[j]]
                    elif (data_where[j+1] - data_where[j])<2 and (data_where[j] - data_where[j-1])>1:
                        indx += [data_where[j]]
                    elif (data_where[j+1] - data_where[j])>1 and (data_where[j] - data_where[j-1])<2:
                        indx += [data_where[j]]                
                
                dif = [] # start and end
                
                for j in range(0,len(indx),2):
                    dif += [(indx[j+1] - indx[j])+1]    
                      
                if np.where(np.array(dif)>=self.dsd_length)[0].shape[0] == 0:
                    self.temp_data["number_density_filtered"][i][:] = np.nan
                    self.temp_data["particle_speed"][i][:] = np.nan
                else:
                     mid = np.where(np.array(data_where)==indx[np.where(np.array(dif)>=self.dsd_length)[0][0]*2])[0]
                     
                     for j in range(1,len(data_where[:mid[0]+1])):
                         if (data_where[:mid[0]+1][-j] - data_where[:mid[0]+1][-j-1])-1 >= self.dsd_gap:
                             self.temp_data["number_density_filtered"][i][:data_where[:mid[0]+1][-j]] = np.nan
                             self.temp_data["particle_speed"][i][:data_where[:mid[0]+1][-j]] = np.nan
                         
                     for j in range(len(data_where[mid[0]:])-1):
                         if (data_where[mid[0]:][j+1] - data_where[mid[0]:][j])-1 >= self.dsd_gap:
                             self.temp_data["number_density_filtered"][i][data_where[mid[0]:][j+1]:] = np.nan
                             self.temp_data["particle_speed"][i][data_where[mid[0]:][j+1]:] = np.nan
    
    def filter_VD(self):
        D_parsi  = np.array([0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.25,3.75,4.250,4.750,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21.5,24.5])
        dD_parsi = np.array([0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.250,0.250,0.250,0.250,0.250,0.50,0.50,0.50,0.50,0.50,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0])
        V_atlas = np.abs(9.65-10.3*np.exp(-0.6*D_parsi))
          
        if self.border==True:   # Sampling volume corrected for border cases:
            Vs_calc = (self.area[0]-((10**-3*D_parsi)/2)) * self.area[1] * V_atlas * 60 * self.time_interval
        else:   # Sampling volume: 
            Vs_calc = self.area[0] * self.area[1] * V_atlas * 60 * self.time_interval
            
        # Create mask
        mask_atlas = self.create_mask_VD(D_parsi)
        
        # Vs_parsivel = np.empty((len(self.temp_data.time),32)) 
        num_par_mask = np.empty((len(self.temp_data.time),32))  
        num_den_filtered = np.empty((len(self.temp_data.time),32))
        
        for time in range(len(self.temp_data.time)):
            # Before masking: sampling volume per diameter = number of particles per diameter / dD * N(D)
            # Vs_parsivel[time] = (np.nansum(self.temp_data["raw_data"][time], axis=0)) / (dD_parsi*(10**self.temp_data["number_density"][time]))
            # Masking all data: new number of particles per diameter
            self.temp_data["raw_data"][time] = np.where(mask_atlas==1, self.temp_data["raw_data"][time], 0.0)
            num_par_mask[time] = np.nansum(self.temp_data["raw_data"][time], axis=0)
    
            # After masking: number density based on raw data
            num_den_dt = (np.nansum(self.temp_data["raw_data"][time], axis=0)) / (dD_parsi*Vs_calc)
            num_den_filtered[time] = np.where(num_den_dt!=0, num_den_dt, np.nan)
        
        # self.temp_data["sampling_volume_parsivel"]=(['time', 'class_diameter'],  Vs_parsivel)  
        self.temp_data["number_particles_filtered"]=(['time', 'class_diameter'],  num_par_mask) 
        self.temp_data["number_density_filtered"]=(['time', 'class_diameter'],  np.log10(num_den_filtered)) 
    
    def create_mask_VD(self, D):
        V_parsi = np.array([0.050,0.150,0.250,0.350,0.450,0.550,0.650,0.750,0.850,0.950,1.1,1.3,1.5,1.7,1.9,2.2,2.6,3.0,3.4,3.8,4.4,5.2,6.0,6.8,7.6,8.8,10.4,12.0,13.6,15.2,17.6,20.8])       
        V_atlas_up = 1.5*np.abs(9.65-10.3*np.exp(-0.6*D))
        V_atlas_down = 0.5*np.abs(9.65-10.3*np.exp(-0.6*D))
        
        mask = np.empty((32,32))
        
        for i in range(0,32):
            for j in range(0,32):
                if V_atlas_down[i] > V_parsi[j] or V_atlas_up[i] < V_parsi[j]:
                    mask[i][j] = False
                else:
                    mask[i][j] = True
        
        mask = mask.T
        return mask

#%% End
