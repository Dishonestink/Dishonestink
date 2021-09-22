# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 08:56:05 2020

@author: Admin
"""

from scipy.optimize import minimize

def objetivo(x):
    f = x[0]*x[3]*(x[0]+x[1]+x[2])+x[2]
    return f

def restriccion1(x):
    r1 = x[0]*x[1]*x[2]*x[3]-25
    return r1

def restriccion2(x):
    suma = 0.0
    for i in range(4):
        suma = suma + x[i]**2
    r2 = suma - 40.0
    return r2

x0 = [1,5,5,1]
b = (1.0,5.0)
fronteras = (b,b,b,b)
r1 = {'type': 'ineq', 'fun': restriccion1}
r2 = {'type': 'eq', 'fun': restriccion2}
res = [r1,r2]

sol = minimize(objetivo,x0,method='SLSQP',bounds=fronteras,constraints=res)




