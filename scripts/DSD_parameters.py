# -*- coding: utf-8 -*-
# Code written by: Bernd Mom, email: bernd.mom@helsinki.fi

import numpy as np
from scipy.integrate import trapz

class DSD_parameters(object):
    """Drop size distrubution parameters
       Drop size distrubution parameters are calculate based on a gamma function. This is done with the third and fourth moment.
       The parameters are Dm: mass-weighted median diameter, D0: median volume diameter, mu: parameter of a gamma drop size distribution 
       and Nw: "intercept" parameter of a normalized gamma drop size distribution.
    
    Attributes:
    	psd_D:		Diameter centers measured by Parsivel.
      	psd_dD:		Diameter width measured by Parsivel.
    """
    
    def __init__(self):
        self.psd_D = np.array([0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.25,3.75,4.250,4.750,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21.5,24.5])
        self.psd_dD = np.array([0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.250,0.250,0.250,0.250,0.250,0.50,0.50,0.50,0.50,0.50,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0]) 
        
    def __call__(self, data):
        self.temp_data = data
        self.calc_para()
        return self.temp_data
    
    def calc_para(self):
        Dm          = np.empty((len(self.temp_data.time)))          # mass-weighted mean diameter (mm)
        mu          = np.empty((len(self.temp_data.time)))          # width of the dsd (-)
        D0          = np.empty((len(self.temp_data.time)))          # median volume diameter based on mu (mm)
        Nw          = np.empty((len(self.temp_data.time)))          # intercept parameter (mm-1 m-3)
        
        for time in range(len(self.temp_data.time)):
            num_den_check = 10**self.temp_data["number_density_filtered"][time].values
            
            if np.nansum(num_den_check)>0:
                num_den_check = np.where(num_den_check>0, num_den_check, 0.0)
                M3 = trapz((self.psd_D**3) * num_den_check, self.psd_D, self.psd_dD)
                M4 = trapz((self.psd_D**4) * num_den_check, self.psd_D, self.psd_dD)
                # Calculate parameters
                Dm[time] = M4/M3                                                        # Dm = M4/M3
                # sigma = sqrt((D-Dm)*M3/M3)
                sigma = np.sqrt(trapz((self.psd_D-Dm[time])**2 * ((self.psd_D**3) * num_den_check), self.psd_D, self.psd_dD) / M3)# mm   
                mu[time] = (Dm[time]**2 / sigma**2) - 4                                 # mu = (Dm**2/sigma**2) - 4               
                D0[time] = ((3.67 + mu[time]) / (4 + mu[time])) * Dm[time]              # D0 = (3.67+mu) / (4+mu) * Dm based on mu                
                Nw[time] = ((3.67)**4 / 6) * (M3 / D0[time]**4)                         # Nw = (3.67)**4 / 6 * (M3 / D0**4) 
            else: 
                Dm[time], mu[time], D0[time], Nw[time] = np.nan, np.nan, np.nan, np.nan
                
        self.temp_data["Dm"]=(['time'],  Dm)
        self.temp_data["mu"]=(['time'],  mu)
        self.temp_data["D0"]=(['time'],  D0)
        self.temp_data["Nw"]=(['time'],  Nw)

#%% End

