# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 14:35:03 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

I0 = 0.0
k1 = 0.6931
k2 = 0.0231

def funciones(t,fi,i):
    x = fi[0]
    y = fi[1]
    if(i==0):
        f = I0 - k1*x
    else:
        f = k1*x - k2*y
    return f

def euler(fi_in,n,dt,t):
    #fi_in = [x,y] del tiempo conocido
    #fi_out = [x,y] del tiempo futuro
    #n es el número de variables
    #t es el instante de tiempo
    fi_out = np.zeros(n)
    for i in range(n):
        fi_out[i] = fi_in[i] + funciones(t,fi_in,i)*dt
    return fi_out

x0 = 200.0
y0 = 0.0
dt = 0.1    #paso de tiempo
tf = 231.0  #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros(it+1)
y = np.zeros(it+1)

x[0] = x0
y[0] = y0

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    [x[i],y[i]] = euler([x[i-1],y[i-1]],2,dt,t[i-1])


c1 = x0 - I0/k1 
c2 = ((1.0/(k2-k1))-(1.0/k2))*I0 + y0 - (k1*x0/(k2-k1))
x_analitica = np.exp(-k1*t)*c1 + I0/k1
y_analitica = (k1*np.exp(-k1*t)*c1/(k2-k1)) + np.exp(-k2*t)*c2 + I0/k2

err_x = abs(x-x_analitica)
err_y = abs(y-y_analitica)

plt.plot(t,x,"--b",label="x_euler")
plt.plot(t,y,"--g",label="y_euler")
plt.plot(t,x_analitica,"b-",label="x_analitica")
plt.plot(t,y_analitica,"g-",label="y_analitica")
plt.legend(loc="upper right")
plt.xlabel('tiempo(min)')
plt.ylabel('masa de medicamento (g)')
plt.grid()


