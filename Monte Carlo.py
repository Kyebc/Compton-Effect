# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 09:44:10 2023

@author: kyebchoo
"""

#%%
import numpy as np
import scipy as sci
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math
import pandas as pd
import pyvisa
import time
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import serial
import re
from pygame import mixer
from time import sleep
from tqdm import tqdm
# from tqdm.auto import tqdm
import keyboard
from datetime import date, datetime
from matplotlib.pyplot import figure

#%%

'''
Notes to code:
    
    1) Mean free path evaluation:
    
        - at 600 keV, for air composing 0.78 N2 + 0.21 O2
        
            \kappa_air = 8.049e-2 cm2/g (mass attenuation coefficient)
            \rho_air   = 1.225e-3 g/cm3
            \mu_air    = \kappa_air * \rho_air = 9.86e-5 /cm (linear attenuation coefficient)
            mfp_air    = 10141 cm ~ 100 meters
            
            ---> since this is far larger distance in experiment, ASSUME NO INTERACTION IN AIR
        
        - at 600keV, for aluminium
        
            \kappa_al  = 7.762e-2 cm2/g
            \rho_air   = 2.7 g/cm3
            \mu_air    = 0.4 /cm
            mfp_air    = 2.5 cm
            
            ---> comparing this to the diameter of the rod, ASSUME that the gamma rays emmited SCATTERS ONLY ONCE in the medium
            
        * data on mass attenuation coeeficient is obtained from NIST XCOM simulation:
            https://physics.nist.gov/cgi-bin/Xcom/xcom3_3
        
    2) Assume that the EMMISION IS COLLINEAR from the source
    
    3) Assume that the SOURCE IS MODELLED BY A CIRCULAR DISC, with UNIFORM PROBABILITY OF EMMISION throughout the disc
    
    4) Assume that SCATTER ANGLE can be modelled using a UNIFORM DISTRIBUTION between 0 degrees to 360 degrees
    
    5) Assume that DETECTOR IS PERFECT
    
    6) Assume that there is an EQUAL CHANCE OF SCATTERING within the rod (not a good assumption due to mfp < diameter)
    
    7) Assume that the model is a 2D model and scattering only happens on a PLANE
    
'''
#%%

def generate_circular_signal(radius = 1):
    xrandom = np.random.uniform(low = -radius, high = radius)
    yrandom = np.random.uniform(low = -radius, high = radius)
    
    # if signal is from the circular cross sectional area, x^2 + y^2 <= R^2
    sum_squares = xrandom**2 + yrandom**2
    radius_squared = radius**2
    if sum_squares <= radius_squared:
        return (xrandom, yrandom)
    else:
        return ('NaN', 'NaN')
    pass



def test_circular(scatter = 1000, radius = 1):
    
    '''
    test by counting number of signals within the circle compared to square of same dimension
    
    should be equal to (\pi * r**2) / (4 * r**2) = 0.7853981634
    '''
        
    within_circle = 0
    x_signal = []
    y_signal = []

    
    for i in range(0, scatter):
        xout, yout = generate_circular_signal(radius)
        if xout != 'NaN' and yout != 'NaN':
            x_signal.append(xout)
            y_signal.append(yout)
            pass
        pass
    
    within_circle = len(x_signal)
    return within_circle/scatter



def compton(E_0, degrees):
    e = 1.602E-19
    E_0 = E_0 * 1000 * e
    theta = (degrees/180) * (np.pi)
    m_e = 9.1093837E-31
    c_0 = 299792458
    E_theta = E_0 / (1 + (E_0 / (m_e * (c_0**2))) * (1 - np.cos(theta)))
    E_theta = E_theta / (1000 * e)
    return E_theta



def compton_full(E_0, degrees, r1 = 3.968, d1 = 150.0, r2 = 9.875, d2= 170.0):
    theta_correct = (np.arctan(r1/d1) + np.arctan(r2/d2)) * 180 / np.pi
    theta_max = degrees + theta_correct
    theta_min = degrees - theta_correct
    E_upper = compton(E_0, theta_min)
    E_middle = compton(E_0, degrees)
    E_lower = compton(E_0, theta_max)
    return (E_upper, E_middle, E_lower)



def monte_carlo(E_0, angle, trials = 1000, radius_emmiter = 3.968, radius_rod = 10.27, distance_rod_detector = 170.0, radius_detector = 9.875):
    
    # dimensions in mm
    
    energy_array = []
    
    for i in range (0, trials):
        
        x, y = generate_circular_signal(radius = radius_emmiter)
        
        if x != 'NaN' and x <= radius_rod:
            x1 = x
        
            path_within_rod = (x1**2 + radius_rod**2)**0.5
            
            x_r = np.random.uniform(low = -path_within_rod, high = path_within_rod)
            
            angle_out = np.random.uniform(low = -180, high = 180)
            
            degree_upper_limit = np.arctan( (radius_detector - x_r) / distance_rod_detector) *180 / np.pi
            # print(degree_upper_limit)
            
            degree_lower_limit = np.arctan( (- radius_detector - x_r) / distance_rod_detector) *180 / np.pi
            # print(degree_lower_limit)
            
            if angle_out >= degree_lower_limit and angle_out <= degree_upper_limit:
                total_deflection = angle - angle_out
                energy = compton(E_0, total_deflection)
                energy_array.append(energy)
                
                pass
            pass
        pass
    
    plt.hist(energy_array, bins = 20)
    
    return energy_array


                
