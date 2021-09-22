# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 18:21:22 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

g = 9.81 #m/s^2 
A = 1.0 #m^2
a = 1e-2 #m^2

def funciones(t,fi,i): 
    #x --> [h,v]
    h = fi[0]
    v = fi[1]
     
    if(i==0):
        f = v
    elif(i==1):
        f = (1.0/(2.0*h))*v**2*((A/a)**2-1) - g
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

dt = 0.0005
tf = 40.0 #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros(it+1)
y = np.zeros(it+1)

x[0] = 1.0
y[0] = -0.043 #m/s

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    [x[i],y[i]] = euler([x[i-1],y[i-1]],2,dt,t[i-1])

plt.plot(t,x,"-g",label="Altura")
#plt.plot(t,y,"-b",label="Velocidad")

plt.legend(loc="upper right")
plt.xlabel('Tiempo(s)')
plt.ylabel('Nivel de agua (m)')
plt.grid()
