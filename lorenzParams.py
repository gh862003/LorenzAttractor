# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:14:18 2016

@author: christiana
"""
import numpy as np

#constant time step  
dt = 0.002
nt = 7000

#variable time step making sure duration and number of timesteps is the same
#as the constant timestep case
#variable_dt = [(4*dt/10)+n*2*dt/100000 for n in range(0,nt)], need nt=40001
'''
variable_dt1 = [0.0002]*(nt/2) 
variable_dt2 = [0.002]*(nt/2)
variable_dt = variable_dt1+variable_dt2



t = np.cumsum(variable_dt)
th = [dt/10]*(nt*10)
'''
#Number from solving via wolfram, added 5 on end to make both durations match
variable_dt = [(0.0001)*1.00064511098971**i for i in range(nt)]
t = np.cumsum(variable_dt)
tf = [dt]*nt
tconst= np.cumsum(tf)

beta  = 8./3.
rho   = 28.
sigma = 10.