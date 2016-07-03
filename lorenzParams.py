# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:14:18 2016

@author: christiana
"""
import numpy as np

#Lorenz Parameters for a 'chaotic solution'
beta  = 8./3.
rho   = 28.
sigma = 10.
#constant time step  
dt = 0.002
#the number of timesteps if this length required to access the region of exponential error growth
#is different for Euler and RK4 and BF so the require different nts
nt_RK = 15000 #7000
nt_FB = 5000 #3500
nt_Eu = 75000 #2500

#variable time step making sure duration and number of timesteps is the same
#as the constant timestep case

############################## RK 4 #######################################
#variables for the *a^n case
#Number from solving via wolfram, added 5 on end to make both durations match
variable_dtRK = [(0.0002)*1.000516604005029**i for i in range(nt_RK)]
variable_timeRK = np.cumsum(variable_dtRK)
const_timeRK = np.cumsum([dt]*nt_RK)
'''

variable_dtRK = [(0.00100029)+i*(2.8567223889e-7)for i in range(nt_RK)]
variable_timeRK = np.cumsum(variable_dtRK)
const_timeRK = np.cumsum([dt]*nt_RK)
'''

########################## Euler Forward ###################################

#variables for the *a^n case
variable_dtEu = [(0.0008)*1.000647930908985**i for i in range(nt_Eu)]
variable_timeEu = np.cumsum(variable_dtEu)
const_timeEu = np.cumsum([dt]*nt_Eu)

'''
#variables for the +n*a case
variable_dtEu = [(0.0010029)+i*(7.9799919968e-7) for i in range(nt_Eu)]
variable_timeEu = np.cumsum(variable_dtEu)
const_timeEu = np.cumsum([dt]*nt_Eu)
print (0.0010029)+(nt_Eu-1)*(7.9799919968e-7)
'''


######################### Backwards forwards ###############################

#variables for the *a^n case
variable_dtFB = [(0.0008)*1.0004627229274844**i for i in range(nt_FB)]
variable_timeFB = np.cumsum(variable_dtFB)
const_timeFB = np.cumsum([dt]*nt_FB)
'''

#variables for the +n*a case
variable_dtFB = [(0.0010029)+i*(5.69934266933e-7) for i in range(nt_FB)]
variable_timeFB = np.cumsum(variable_dtFB)
const_timeFB = np.cumsum([dt]*nt_FB)
'''


#variable_dt = [(4*dt/10)+n*2*dt/100000 for n in range(0,nt)], need nt=40001
'''
variable_dt1 = [0.0002]*(nt/2) 
variable_dt2 = [0.002]*(nt/2)
variable_dt = variable_dt1+variable_dt2



t = np.cumsum(variable_dt)
th = [dt/10]*(nt*10)
'''