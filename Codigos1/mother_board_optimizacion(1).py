# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 08:56:05 2020

@author: Admin
"""


from scipy.optimize import minimize

def objetivo(x):
    xa = x[0]
    xb = x[1]
    xc = x[2]
    f = -(55*18*xa + 50*18*xb + 50*21*xc)
    return f

def r_costo(x):
    xa = x[0]
    xb = x[1]
    xc = x[2]
    r1 = 12e6 - 1e5*(4*xa + 6*xb + 7*xc)
    return r1

def r_operarios(x):
    xa = x[0]
    xb = x[1]
    xc = x[2]
    r2 = 100 - (3*xa + 6*xb + 6*xc)
    return r2

def r_espacio(x):
    xa = x[0]
    xb = x[1]
    xc = x[2]
    r3 = 30 - (3*xa - xb + xc)
    return r3

x0 = [10,5,5]
rango1 = [0,10]
rango2 = [0,20]
rango3 = [0,18]
fronteras = (rango1,rango2,rango3)

r1 = {'type': 'ineq', 'fun': r_costo}
r2 = {'type': 'ineq', 'fun': r_operarios}
r3 = {'type': 'ineq', 'fun': r_espacio}

r4 = {'type':'eq','fun': lambda x : max([x[i]-int(x[i]) for i in range(len(x))])}

res1 = [r1,r2,r3]
res2 = [r1,r2,r3,r4]

sol1 = minimize(objetivo,x0,method='SLSQP',bounds=fronteras,constraints=res1)
sol2 = minimize(objetivo,x0,method='SLSQP',bounds=fronteras,constraints=res2)

print(sol1)
print("----------\n----------")
print(sol2)

