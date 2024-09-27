# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 09:39:08 2023

@author: 11383
"""

import numpy as np
import scipy as sci
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math
import pandas as pd
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

me = 9.11e-31
c = 299792458
e = 1.602e-19
N_15cm = 173483
N_15cm_662 = 60886*3/0.9546 #15cm_1:36402(1 sigma) #15cm_2: 43463(1 sigma); 60886(2 sigma)

r_Al = 10.09e-3 #the error is 0.03e-3/2
h_Al = 0.1        #the error is 0.01e-3
V_Al = np.pi*(r_Al**2)*h_Al
mp = 1.67262192e-27
am = 26.981539
ma = mp * am
m_Al = 86.54e-3   #the error is 0.01e-3
na = m_Al/ma
ne = na*13

d = 0.195 #source-detector distance
r = 0.01975/2 #detector radius
A_detector = np.pi*r**2
denominator = N_15cm_662/(300*A_detector)*661.7
theta = np.arctan(1.975/39) #np.arctan(1.985/38.6)- np.arctan(1.975/39)
SA = 2*np.pi*(1-np.cos(theta))

flux = 3.7e6/(4*np.pi*0.15**2)   
A = (0.15)**2*np.pi*2
    
def error_calculation(E0=661.7,file = 'N_sc.csv'):
    data = pd.read_csv(file)
    nsc_list = data['N_act(1 sigma)']
    E_list = data['Energy']
    e_E_list = data['Energy error']
    n = N_15cm_662/300
    e_SA = 2.54e-4
    e_n = np.sqrt(60886*3)/0.9546/300
    e_A = 6.87e-6
    y_err = []
    
    
    for i in range(len(nsc_list)):
        nsc = nsc_list[i]/1200/0.6826
        e_nsc = np.sqrt(nsc_list[i])/1200/0.6826
        E = E_list[i]
        e_E = e_E_list[i]
        pe_nsc = e_nsc/nsc
        pe_A = e_A/A
        pe_E = e_E/E
        pe_n = e_n/n
        pe_SA = e_SA/SA
        total_perr = np.sqrt(pe_nsc**2+pe_A**2+pe_E**2+pe_n**2+pe_SA**2)
        y_err.append(total_perr)
    return y_err

    
def gammaE(E0, theta=0.0):
    theta_0 = theta*np.pi/180
    E_0 = E0 * e *1000
    deno = 1 + (E_0/(me*c**2))*(1-np.cos(theta_0))
    E_1 = E_0/deno
    E = E_1/(e*1000)
    return E

def fitting_curve(angle,A,b):
    theta = angle*np.pi/180
    E0 = 661.7
    E = gammaE(E0, theta = angle)
    ra = E/E0
    c1 = 7.94e-30
    DCS = A*ne*0.5*c1*(ra**2)*(ra+1/ra-np.sin(theta)**2)+b
    return DCS

    

def DCS(filename = 'N_sc.csv',t = 1200):
    data = pd.read_csv(filename)
    E1 = data['Energy']
    N_fitr = np.array(data['N_fit (range)'])
    N_fit1 = np.array(data['N_fit(1 sigma)'])
    N_act1 = np.array(data['N_act(1 sigma)'])
    N_act2 = np.array(data['N_act(2 sigma)'])
    DCS_fr0 = []
    DCS_f10 = []
    DCS_a10 = []
    DCS_a20 = []
    for i in range(len(N_fitr)): 
        
        DCS_fr = N_fitr[i]*E1[i]/(t*SA*denominator)
        DCS_f1 = N_fit1[i]*E1[i]/(t*SA*denominator)/0.6826
        DCS_a1 = N_act1[i]*E1[i]/(t*SA*denominator)/0.6826
        DCS_a2 = N_act2[i]*E1[i]/(t*SA*denominator)/0.9546
        
        DCS_fr0.append(DCS_fr)
        DCS_f10.append(DCS_f1)
        DCS_a10.append(DCS_a1)
        DCS_a20.append(DCS_a2)
    
    return[DCS_f10, DCS_a10, DCS_a20]
    
def K_N_formula(E,angle):
    theta = angle*np.pi/180
    E0 = 661.7
    c1 = 7.94e-30
    ra = E/E0
    dcs = 0.5*c1*(ra**2)*(ra+1/ra-np.sin(theta)**2)
    DCS = dcs*ne
    return [dcs,DCS]


def DCS_plot(filename1 = 'N_sc.csv',filename2 = 'DCS results.csv'):
    data = pd.read_csv(filename1)
    theta = data['theta']
    data = DCS(filename = filename1, t = 1200)
    DCS_act1 = data[1]
    err_DCS_act1 = []
    yerr = error_calculation()
    for i in range(len(DCS_act1)):
        err_DCS = DCS_act1[i]*yerr[i]
        err_DCS_act1.append(err_DCS)
        
    print(err_DCS_act1)
    DCS_act2 = data[2]
    DCS_fit = data[0]
    DCS_t1 = pd.read_csv(filename2)['K-N DCS']
    theta_t = np.linspace(0,140,140)
    DCS_t = []
    for i in range(len(theta_t)):
        Energy = gammaE(E0 = 661.7, theta = theta_t[i])
        DCS_T = K_N_formula(E=Energy,angle = theta_t[i])[1]
        DCS_t.append(DCS_T)
        
    
    theta_fit = theta
    DCS_fit = DCS_act1
    p0 = [0.5,0.00002]
    popt, pcov = curve_fit(fitting_curve, theta_fit, DCS_fit, p0 = p0)
    print(popt,pcov)
    theta_fitted = np.linspace(0,140,140)
    DCS_fitted = fitting_curve(angle = theta_fitted, A = popt[0], b = popt[1])
    DCS_expected = fitting_curve(angle = theta, A = popt[0], b = popt[1])
    print(sci.stats.chisquare(DCS_act1,DCS_expected,ddof = 2))
    yerr = error_calculation()
    
    
    plt.plot(theta_fitted, DCS_fitted, label = 'fit')
    
    
    plt.plot(theta_t,DCS_t,label = 'theoretical', color = 'r')
    plt.scatter(theta,DCS_act1, label = '1 sigma')
    plt.errorbar(theta, DCS_act1, yerr = err_DCS_act1,fmt = 'none')
    # plt.scatter(theta,DCS_act2, label = '2 sigma')
    # plt.scatter(theta,DCS_fit, label = 'fit')
    
    plt.scatter(theta,DCS_t1, label = 'theoretical DCS with measured Energy')
    plt.legend()
    plt.show()
    
    


    
