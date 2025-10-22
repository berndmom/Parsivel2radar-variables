# -*- coding: utf-8 -*-
# Code written by: Bernd Mom, email: bernd.mom@helsinki.fi

import numpy as np
import xarray as xr
import scipy.integrate as integrate
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import GammaPSD, PSDIntegrator
from pytmatrix import orientation, radar, refractive, tmatrix_aux, scatter

class radar_parameters(object):
    """pytmatrix + mean Doppler velocity calculations
    
    Attributes:
    	file: 		File path to the stored look-up table.
    	indx:		Number of times to run the scatterer for new triplets.
    	D0:		Input D0.
    	mu:		Input mu.
    	Nw:		Input Nw.
    	limits:		Filter the D0, mu and Nw based on minimum limits.
    	wavelength:	The wavelength used in the scatterer. E.g. 'W', 'C', 'Ka', ...
    	geometry:	The geometry used in the scatterer. E.g. 'horizontal' or 'vertical'.
    			Horizontal geometry: reflectivity in dB and dBZ (Zh/Zh_dBZ), mean Doppler velocity (mdv),
    			rainfall rate (R) and specific attenuation (Ai). 
    			Vertical geometry: reflectivity in dB and dBZ (Zh/Zh_dBZ), mean Doppler velocity (mdv),
    			(Zdr) and specific reflectivity (Kdp)
    """
    
    def __init__(self, file, indx=[], D0=0, mu=0, Nw=0, limits=[.1, -1, 100]):
        self.file = file
        self.indx = indx
        self.D0 = D0
        self.mu = mu
        self.Nw = Nw
        self.limits = limits
        self.wavelength = "W"
        self.geometry = "vertical"        
                
    def prep_scatter(self):    
        self.create_dataset()
        self.filter_dsd()
        
        if self.wavelength == 'W' or self.wavelength == 'w':
            tmatrix_wl = tmatrix_aux.wl_W
        elif self.wavelength == 'C' or self.wavelength == 'c':
            tmatrix_wl = tmatrix_aux.wl_C
        else:
            return print('No wavelength given')
            
        if self.geometry == 'horiz' or self.geometry == 'horizontal':
            tmatrix_geom_back = tmatrix_aux.geom_horiz_back,
            tmatrix_geom_forw = tmatrix_aux.geom_horiz_forw,
        elif self.geometry == 'vert' or self.geometry == 'vertical':
            tmatrix_geom_back = tmatrix_aux.geom_vert_back,
            tmatrix_geom_forw = tmatrix_aux.geom_vert_forw,
        else:
            return print('No geometry given')   
            
        print('Wavelength: ' , self.wavelength, '; Geometry: ', self.geometry, '\n File: ', self.file)
        
        self.count = 0

        for i in self.indx:
            if i in self.indx_sel:
                self.calc_scatter(tmatrix_wl, tmatrix_geom_back, tmatrix_geom_forw)
            self.count+=1
            
        self.data_temp = {"index": {"dims": ("index"), "data": self.indx},}
        for param in self.parameters:     
            self.data_temp[param] = {"dims": ("index"), "data": self.data[param]}
        
        self.data_temp["Zh_dBZ"] = {"dims": ("index"), "data": 10 * np.log10(self.data["Zh"]/1)}
        if "Zdr" in list(self.data.keys()):
            self.data_temp["Zdr_dB"] = {"dims": ("index"), "data": 10 * np.log10(self.data["Zdr"])}
        
        self.data = xr.Dataset.from_dict(self.data_temp)
        print('Finished scatterer')
        return self.data
                
    def calc_scatter(self, tmatrix_wl, tmatrix_geom_back, tmatrix_geom_forw):
        """ Scatterer class in pytmatrix
	    For more information refer to pytmatrix
        """
        scatterer = Scatterer(
            wavelength=tmatrix_wl, m=refractive.m_w_10C[tmatrix_wl]
        )
        scatterer.psd_integrator = PSDIntegrator()
        scatterer.psd_integrator.axis_ratio_func = (
            lambda D: 1.0 / tmatrix_aux.dsr_thurai_2007(D)
        )
        scatterer.psd_integrator.D_max = (
            2.5 * self.D0[self.count]
        )  # maximum diameter used for DSD integration
    
        # defining raindrop orientations 0 elevation angle and 7 deg std of the canting angle distributions
        scatterer.or_pdf = orientation.gaussian_pdf(7.0)
        scatterer.orient = orientation.orient_averaged_fixed
    
        scatterer.psd_integrator.geometries = (
            tmatrix_geom_back,
            tmatrix_geom_forw,
        )
        
        # Use look up table or calculate the initial look up table
        scatterer.psd_integrator.load_scatter_table(fn=self.file)
        
        # here you input the DSD parameters
        scatterer.psd = GammaPSD(D0=self.D0[self.count], Nw=self.Nw[self.count], mu=self.mu[self.count])
    
        # defining the scattering geometry, for reflectivity computation it should be backscattering          
        if self.geometry == 'horiz' or self.geometry == 'horizontal':
            scatterer.set_geometry(tmatrix_geom_back[0])
            self.data["Zh"][self.count] = radar.refl(scatterer)
            self.data["zdr"][self.count] = radar.Zdr(scatterer)
            scatterer.set_geometry(tmatrix_geom_forw[0])
            self.data["kdp"][self.count] = radar.Kdp(scatterer)
        elif self.geometry == 'vert' or self.geometry == 'vertical':
            scatterer.set_geometry(tmatrix_geom_back[0])              
            self.data["mdv"][self.count] = self.radar_mdv(scatterer, True, 90)
            self.data["Zh"][self.count] = radar.refl(scatterer)
            scatterer.set_geometry(tmatrix_geom_forw[0])
            self.data["Ai"][self.count] = radar.Ai(scatterer)
        
        def func(d):
            v = np.abs(9.65-10.3*np.exp(-0.6*d));
            int_part = d**3 * v * scatterer.psd(d) 
            return int_part
    
        int_rain = integrate.quad(func, 0, scatterer.psd_integrator.D_max)
        self.data["R"][self.count] = 0.6 * np.pi * 10**-3 * int_rain[0]
    
    def radar_mdv(self, scatterer_mdv, h_pol, el):
        """Mean Doppler velocity
	   Added mean Doppler velocity calculations using the pytmatrix Scatterer. 
        """
        Z_mdv = {}
        Z_table = scatterer_mdv.psd_integrator.__dict__['_Z_table'] 
        psd_D = scatterer_mdv.psd_integrator.__dict__['_psd_D']
        psd_w = scatterer_mdv.psd(psd_D)
        geometry = scatterer_mdv.get_geometry()
        V_atlas = np.abs(9.65-10.3*np.exp(-0.6*psd_D))
        
        for geom in scatterer_mdv.psd_integrator.__dict__['geometries']:
            Z_geom = Z_table[geom]
            if h_pol:
                D_xsect = 2 * np.pi * (Z_geom[0,0,:] - Z_geom[0,1,:] - Z_geom[1,0,:] + Z_geom[1,1,:]) # radar xsect
            else:
                D_xsect = 2 * np.pi * (Z_geom[0,0,:] + Z_geom[0,1,:] + Z_geom[1,0,:] + Z_geom[1,1,:]) # radar xsect

            Z_mdv[geom] = -(integrate.trapz(V_atlas * np.sin(el * np.pi / 180.) * D_xsect * psd_w, psd_D) / integrate.trapz(D_xsect * psd_w, psd_D))       
        return Z_mdv[geometry]
    
    def filter_dsd(self):
        """Filter D0, mu and Nw based on minimum values.
        """
        self.mu = np.where(self.D0>self.limits[0], self.mu, np.nan)
        self.Nw = np.where(self.D0>self.limits[0], self.Nw, np.nan)
        self.D0 = np.where(self.D0>self.limits[0], self.D0, np.nan)
        self.D0 = np.where(self.mu>self.limits[1], self.D0, np.nan)
        self.Nw = np.where(self.mu>self.limits[1], self.Nw, np.nan)            
        self.mu = np.where(self.mu>self.limits[1], self.mu, np.nan)
        self.D0 = np.where(self.Nw>self.limits[2], self.D0, np.nan)
        self.mu = np.where(self.Nw>self.limits[2], self.mu, np.nan)            
        self.Nw = np.where(self.Nw>self.limits[2], self.Nw, np.nan)  
        self.indx_sel = self.indx[np.where(self.D0>0)]
    
    def create_dataset(self):   
        """Create a data dictionary to store the data.
        """
        if self.geometry == 'horiz' or self.geometry == 'horizontal':
            self.parameters = ["Zh", "R", "Zdr", "Kdp"]
        elif self.geometry == 'vert' or self.geometry == 'vertical':
            self.parameters = ["Zh", "mdv", "R", "Ai"]
        
        self.data = {}
        data_array = np.empty(self.indx.shape)
        data_array[:] = np.nan
        
        for i in range(len(self.parameters)):
        
            self.data[self.parameters[i]] = data_array.copy()

#%% End

