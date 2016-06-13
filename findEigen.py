# -*- coding: utf-8 -*-
"""
This module contains exclusively designed to solve and plot the 
eigen values and numerical stability amplitude of the Runge-Kutta
Scheme when used to solve the recharge oscillator problem. It conatins
3 methods: numericalAmp, eigenSolver and eigenPlots.
"""
from __future__ import absolute_import, division
import numpy as np
import matplotlib.pyplot as plt







def eigenSolver(lorenzParams,dt):
    '''takes 6 inputs: b_o, mu, a,c,y,r that are all single numbers that
    represent constants in the recharge oscillator problem, a full description
    of their meaning can be found in the documentation of method 'Osc' in 
    project1.py. 
    This function sets up and solves for the eigen values of the
    matrix equation that defines the
    coupled ODEs of the recharge oscillator system. It returns the two eigen values.
    '''
    sigma = lorenzParams['sigma']
    rho  = lorenzParams['rho']
    beta = lorenzParams['beta']
    '''
    #Define Matrix
    a_11 = 1+ sum([(-1)**i *sigma**i *dt**i *1/math.factorial(i) for i in range(1,5)])
    a_12 = sum([(-1)**(i+1) *sigma**i *dt**i *1/math.factorial(i) for i in range(1,5)]) 
    a_13 = 0   
    a_21 = rho*sum([(-1)**(i+1) *dt**i *1/math.factorial(i) for i in range(1,5)])
    a_22 = 1+sum([(-1)**(i) *dt**i *1/math.factorial(i) for i in range(1,5)])
    a_23 = 0
    a_31 = 0
    a_32 = 0
    a_33 = 1+sum([(-1)**(i) *beta**i*dt**i *1/math.factorial(i) for i in range(1,5)])
    '''
    
    A =np.matrix([[-sigma, sigma,0],[rho, -1,0],[0,0,-beta]])
    e =np.linalg.eigvals(A) 
    return e
    

def numericalAmp(e,dt):
    '''Takes two inputs, e, a set of eigen values and dt a time step.
    It returns the amplification factor corresponding to these eign values.
    '''
    u = dt*e
    #calculate amplification factor
    A_euler = u
    musRK = 0.5*u**2 +(1/6)*u**3+(1/24)*u**4
    return musRK, A_euler 

def eigenPlots(lorenzParams):  
    '''
    Going to try and find the stability of the lorenz attractor that I think
    depends on the initial conditions
    '''
    B = lorenzParams['beta']
    dt_range = np.linspace(0.000001,0.05,80)
    plt.figure(figsize=(5,4))
    #calculate and plot all background values
    e =eigenSolver(lorenzParams,0)
    print e
    allmu =[]
    A_E = []
    for dt in dt_range:
        e_abs = [abs(i) for i in e]
        k = max(e_abs)
        plt.plot(dt,k*dt, '.')
        '''
        k = max(e)
        
        #k =e[ e_abs.index(max(e_abs))]
        mu = abs(e*dt)
        mus,A_euler= numericalAmp(e,dt)
        mus = abs(mus)
        A_euler = abs(A_euler)
        
        if mu[1]>mus[1]:
            print dt
        
        #A_VN = abs(np.exp(B*dt)*(dt-B*dt)*(1-0.5*B*dt+(1/6)*B**2*dt**2-(1/24)*B**3*dt**3))
        plt.plot(mu,mus,'.g')
        
        
        allmu.append(mu)
        A_E.append(max(A_euler))
        '''

        
        
        
   #add horizontal line at y=1 to show limit of numerical stability 
    xhigh = 3.5
    
    
    
    #plt.plot([0,xhigh],[0,xhigh],'-k')
    #plt.xlabel('$\lambda \Delta$ t')
    #plt.ylabel('1/2 $ \lambda \Delta$ t$^{2}$ +1/6 $ \lambda \Delta$ t$^{3}$ +1/24 $ \lambda \Delta$ t$^{4}$')
    #plt.show()    
    #plt.clf()
    #plt.plot(dt_range,A_E)
    plt.xlabel('dt')
    plt.ylabel('max(|lambda|) dt')
    plt.axhline(y=1)
    plt.show()

    
def main():

   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   eigenPlots(lorenzParam)
    
main()
    
    

