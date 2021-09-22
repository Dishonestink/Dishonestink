# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 14:28:08 2020

@author: Admin
"""
import numpy as np
from scipy.optimize import minimize

def objetivo(x):
    #x --> [d,h,e]
    d = x[0]
    h = x[1]
    e = x[2]
    A1 = (np.pi/4)*d**2
    A2 = np.pi*d*h
    return (100.0 + (0.1*(1000.0*e - 1.0)*100.0))*(2.0*A1+A2)

def geometria(x):
    d = x[0]
    h = x[1]
    return h - 2.0*d

def volumen(x):
    d = x[0]
    h = x[1]
    return h*(np.pi/4)*d**2 - 0.5e-3

x0 = [0.09,0.05,0.009]
b1 = [0.06,0.09] # rango de la variable d
b2 = [0.05,0.2] # rango de la variable h
b3 = [0.001,0.01] # rango de la variable e

fronteras = (b1,b2,b3)

r1 = {'type': 'ineq', 'fun': geometria}
r2 = {'type': 'eq', 'fun': volumen}

res = [r1,r2]

sol = minimize(objetivo,x0,method='SLSQP',bounds=fronteras,constraints=res)














