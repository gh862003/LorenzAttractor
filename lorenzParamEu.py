# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:00:42 2016

@author: christiana
Contains the patrameters to run the lorenz attractor using a euler foreward scheme as it requires
differenct parameters to that of RK, kept in separate 
"""

import numpy as np

#constant time step  
dt = 0.002
nt = 100

#Number from solving via wolfram, added 5 on end to make both durations match
variable_dt = [(0.0001)*1.0464753417137**i for i in range(nt)]
t = np.cumsum(variable_dt)
tf = [dt]*nt
tconst= np.cumsum(tf)

beta  = 8./3.
rho   = 28.
sigma = 10.