# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:36:12 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

def funciones(t,fi,i):
    x = fi[0]
    y = fi[1]
    if(i==0):
        f = -(1.0/2.0)*x - (1.0/8.0)*y + 0.5*np.exp(-0.5*t)
    else:
        f = 2.0*x - (1.0/2.0)*y
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

x0 = 0.0
y0 = 0.0
dt = 0.01    #paso de tiempo
tf = 20.0    #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros(it+1)
y = np.zeros(it+1)

x[0] = x0
y[0] = y0

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    [x[i],y[i]] = euler([x[i-1],y[i-1]],2,dt,t[i-1])

c1 = y0 - 4.0
c2 = 4.0*x0

x_analitica = -0.25*c1*np.exp(-0.5*t)*np.sin(0.5*t)+0.25*c2*np.exp(-0.5*t)*np.cos(0.5*t)
y_analitica = c1*np.exp(-0.5*t)*np.cos(0.5*t)+c2*np.exp(-0.5*t)*np.sin(0.5*t)+4.0*np.exp(-0.5*t)

plt.plot(t,x,"-b",label="x_euler")
plt.plot(t,y,"-g",label="y_euler")
plt.plot(t,x_analitica,"b--",label="x_analitica")
plt.plot(t,y_analitica,"g--",label="y_analitica")
plt.legend(loc="upper right")
plt.xlabel('tiempo(s)')
plt.grid()

