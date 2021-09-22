# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 08:37:55 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

def funciones(t,fi,i):
    x = fi[0]
    y = fi[1]
    if(i==0):
        f = -(x/3.0) + (y/12.0)
    else:
        f = (x/3.0) - (y/3.0)
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

x0 = 20.0   #en gramos
y0 = 5.0   #en gramos
dt = 0.5    #paso de tiempo en minutos
tf = 30.0   #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros(it+1)
y = np.zeros(it+1)

x[0] = x0
y[0] = y0

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    [x[i],y[i]] = euler([x[i-1],y[i-1]],2,dt,t[i-1])

x_analitica = -((y0-2.0*x0)/4.0)*np.exp(-0.5*t) + ((y0+2.0*x0)/4.0)*np.exp(-t/6.0)
y_analitica = ((y0-2.0*x0)/2.0)*np.exp(-0.5*t) + ((y0+2.0*x0)/2.0)*np.exp(-t/6.0)

err_x = abs(x_analitica-x)
err_y = abs(y_analitica-y)

plt.plot(t,x,"--b",label="x_euler")
plt.plot(t,x_analitica,"-b",label="x_analitica")
plt.plot(t,y,"--g",label="y_euler")
plt.plot(t,y_analitica,"-g",label="y_analitica")
plt.legend(loc="upper right")
plt.xlabel('tiempo(min)')
plt.ylabel('Masa de la sustancia (g)')
plt.grid()
plt.show()

# plt.plot(t,err_x)
# plt.plot(t,err_y)
    

        
    
    