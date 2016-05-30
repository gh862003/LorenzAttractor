# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:14:18 2016

@author: christiana
"""
import numpy as np

#constant time step  
dt = 0.001
#variable time step making sure duration and number of timesteps is the same
#as the constant timestep case
variable_dt =np.linspace(0.0001,0.0139,100000)


nt = 100000

x_o = -10.
y_o = -10.
z_o =  10.




beta  = 8./3.
rho   = 28.
sigma = 10.