# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 21:08:29 2020

@author: Admin
"""

import numpy as np

def fx(x):
    f = ((x**3)-(9*x**2)+25*x*(1+(((np.sin(x))**2)/25))+(x*(np.cos(x))**2)-24)
    return f    

def dfx(x,dx):
    df = (fx(x + 0.5*dx) - fx(x - 0.5*dx))/dx
    return df

dx = 1e-4 
x_old = 1.0   #valor semilla  1.0
fx_new = 1.0   
diff = 1.0  #diferencia entre el x_new y el x_old
eps_y = 1e-4
eps_x = 1e-3
i = 0

while((abs(fx_new)>eps_y)and(abs(diff)>eps_x)):
    i = i + 1
    x_new = x_old - (fx(x_old)/dfx(x_old,dx))
    fx_new = fx(x_new)
    diff = (x_new - x_old)
    print(x_old,x_new,fx_new,diff)
    x_old = x_new

print("Convergió a ",i,"iteraciones")
print("Una raíz es: ",x_new)

