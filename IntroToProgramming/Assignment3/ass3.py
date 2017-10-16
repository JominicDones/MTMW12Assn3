# -*- coding: utf-8 -*-
"""

File that runs code necessary for MTMW12 Assignment 3. I wasn't able to get 
nearly as much done as I'd hoped as I'm still getting to grips with Python and 
Git and I also missed the class and wasn't able to catch up. 

Hopefully the code is styled correctly and commented okay at least.

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def pressureCalc(_P_a, _P_b, _y, _L):
    """
    Calculates the pressure as a function of y
    
    Input: _P_a : pressure at ground
    Input: _P_B : pressure difference at maximum height
    Input: _y : height
    Input: _L : characteristic length scale
    Variable: pressure : the pressure as a function of height(y)
    """
    pressure = np.zeros(len(_y))
    for x in xrange(0, len(_y)):
        pressure[x] = _P_a + _P_b * np.cos(_y[x] * np.pi / _L)
        
    #Plot the pressure as a function of height (to check)
    plt.figure(1, figsize=(10,4))
    plt.plot(_y,pressure)
    plt.ylabel("Pressure (Pa)")
    plt.xlabel("Height (m)")
    plt.title("Variation of Pressure with Height")
    plt.show()
    
    return pressure





def pressureGradient(_y, _P, N):
    """
    Calculates the gradient of the pressure, dP/dy using  centred, 2nd order 
    finite difference formula. The boundary points are calculated using
    forward and backward differences.
    
    Input: _y : height
    Input: _P : pressure
    Variable: delta_y : step/distance between points
    Variable: P_prime : grad of pressure
    """
    #Calculate the non-boundary values of pressure gradient
    P_prime = np.zeros(len(_P))
    delta_y = (_y[len(_y)-1] - _y[0])/N
    for x in xrange(1, len(_y)-1):
        P_prime[x] = (_P[x+1] - _P[x-1])/(2 * delta_y)
        
    #Calculate the boundary values
    P_prime[0] = (_P[1] - _P[0])/delta_y
    P_prime[len(P_prime)-1] = (_P[len(P_prime)-1] - _P[len(P_prime)-2])/delta_y
    
    '''    
    #Plot the pressure gradient as a function of height (to check)
    plt.figure(2, figsize=(10,4))
    plt.plot(_y,P_prime)
    plt.ylabel("Pressure gradient (Pa/m)")
    plt.xlabel("Height (m)")
    plt.title("Variation of Pressure gradient with Height")
    plt.show()
    '''
    
    return P_prime

def errPlot(_y, _P_b, _L, _P_prime, _u, _rho, _f):
    """
    Plots the analytical and numerical solutions of P', and plots the errors.

    Input: _y : height
    Input: _P_b : pressure difference at the top
    Input: _L : length scale
    Input: _P_prime : numerical values of P' 
    Input: _u : numerically calculated geostrophic wind
    Param: anal_Pprime : the analytical value of P'
    Variable: errors : the error between the analytical and numerical
    values of P'
    Variable: Relerrors : the relative error between the analytical and numerical
    values of P'
    Variable: windErr : the error between the analytical and numerical
    values of u
    Variable: relWindErr : the relative error between the analytical and numerical
    values of u
    
    """
    #Calculate the analytical solution for P'
    anal_Pprime = np.zeros(len(_P_prime))
    for x in xrange(0, len(anal_Pprime)):
        anal_Pprime[x] = - _P_b * (np.pi / _L) * np.sin(_y[x] * np.pi / _L)
        
    #Plot the pressure gradient as a function of height,\
    #analytical and numerical
        
    plt.figure(3, figsize=(10,4))
    plt.plot(_y,_P_prime)
    plt.plot(_y,anal_Pprime)
    plt.ylabel("Pressure gradient (Pa/m)")
    plt.xlabel("Height (m)")
    plt.title("Analytical vs. Numerical pressure gradient solutions")
    plt.show()
    
    '''
    #Calculate the relative error 
    errors = np.zeros(len(_P_prime))
    Relerrors = np.zeros(len(_P_prime))
    for x in xrange(0, len(errors)):
        errors[x] = abs(anal_Pprime[x] - _P_prime[x])
        Relerrors[x] = 100 * errors[x] / anal_Pprime[x]
            
    plt.figure(4, figsize=(10,4))
    plt.plot(_y,errors)
    plt.ylabel("relative error %")
    plt.xlabel("Height (m)")
    plt.title("Error of numerical solution")
    plt.show()
    
    plt.figure(5, figsize=(10,4))
    plt.plot(_y,Relerrors)
    plt.ylabel("relative error")
    plt.xlabel("Height (m)")
    plt.title("Relative error of numerical solution")
    plt.show()
    '''
    #Calculate the analytic value for the wind
    anal_u = np.zeros(len(_u))
    coeff = -1 * 1/_rho * 1/_f
    for x in xrange(0, len(_y)):
        anal_u[x] = coeff * anal_Pprime[x]
    
    #plot the numerical and analytical solutions for wind
    plt.figure(7, figsize=(10,4))
    plt.plot(_y,_u)
    plt.plot(_y,anal_u)
    plt.ylabel("Wind Speed")
    plt.xlabel("Height (m)")
    plt.title("Analytical vs. Numerical Wind")
    plt.show()
    
    #Calculate the errors on u
    #Calculate the relative error 
    windErr = np.zeros(len(_u))
    relWindErr = np.zeros(len(_u))
    for x in xrange(0, len(windErr)):
        windErr[x] = abs(anal_u[x] - _u[x])
        relWindErr[x] = 100 * windErr[x] / anal_u[x]
            
    plt.figure(8, figsize=(10,4))
    plt.plot(_y,windErr)
    plt.ylabel("absolute error")
    plt.xlabel("Height (m)")
    plt.title("Error of numerical solution of u")
    plt.show()
    
    plt.figure(9, figsize=(10,4))
    plt.plot(_y,relWindErr)
    plt.ylabel("relative error %")
    plt.xlabel("Height (m)")
    plt.title("Relative error of numerical solution of u")
    plt.show()
    
    
    
def windCalc(_y, _P_prime, _f, _rho):
    """
    Calculates the windspeed.
    
    Input: _y : height
    Input: _P_prime : pressure gradient
    Input: _f : coriolis parameter
    Input: _rho : density
    Variable: u : Geostropic wind 
    """    
    u = np.zeros(len(_y))
    coeff = -1 * 1/_rho * 1/_f
    for x in xrange(0, len(_y)):
        u[x] = coeff * _P_prime[x]
        
    plt.figure(6, figsize=(10,4))
    plt.plot(_y,u)
    plt.ylabel("Wind Speed")
    plt.xlabel("Height (m)")
    plt.title("Geostrophic Wind Relation")
    plt.show()
    
    return u
    
def main():
    
    """
    Main function for the programme. 
    
    Param: P_a : pressure at ground
    Param: P_b : pressure at maximum height
    Param: f : coriolis parameter
    Param: rho : density
    Param: L : characteristic length scale (??)
    Param: N : maximum height (meters)
    Param: y: height
    Variable: P : pressure
    Variable: P_prime : pressure gradient
    Variable: u : Geostrophic wind
    
    """
    #Define variables
    P_a = 1e5
    P_b = 200
    f = 1e-4
    rho = 1
    L = 2.4e6
    N = 10
    ymin = 0
    ymax = 1000000
    y = range(ymin, ymax + int((ymax - ymin) / N), int((ymax - ymin) / N))
    
    #Calculate the pressure as a function of height
    P = pressureCalc(P_a, P_b, y, L)
    
    #Calculate the gradient of the pressure
    P_prime = pressureGradient(y,P,N)
    
    #Calculate and plot the Geostrophic Wind Relation
    u = windCalc(y, P_prime, f, rho)    
    
    #Plot the differnce between the analytical and numerical solutions for P'
    errPlot(y,P_b,L,P_prime, u, rho, f)    
      
    '''
    The error is bigger at the boundaries, where the finite difference formula 
    is only 1st order accurate, as opposed to 2nd. 
    
    An experiment to test this would be to measure the error relative to the 
    step size and see if it scales as the square of the order of accuracy.
    
    (2) A more accurate differentiation scheme would be a 4th order accurate 
    central difference scheme. 
    '''
    
main()
