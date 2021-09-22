# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 14:35:03 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

I0 = 5.0
k1 = 0.6931
k2 = 0.0231
t1 = 5.0
t2 = 10.0
t3 = 15.0

def dosis(t):
    I = I0*(escalon(t)-escalon(t-t1)+escalon(t-t2)-escalon(t-t3))
    return I

def escalon(t):
    if(t>=0.0):
        u = 1.0
    else:
        u = 0.0
    return u

def funciones(t,fi,i):
    x = fi[0]
    y = fi[1]
    if(i==0):
        f = dosis(t) - k1*x
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

x0 = 0.0
y0 = 0.0
dt = 0.1  #paso de tiempo
tf = 60.0 #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros(it+1)
y = np.zeros(it+1)
x_analitica = np.zeros(it+1)
y_analitica = np.zeros(it+1)

x[0] = x0
y[0] = y0

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    [x[i],y[i]] = euler([x[i-1],y[i-1]],2,dt,t[i-1])

c1 = x0 + (I0/k1)*(np.exp(k1*t1)-np.exp(k1*t2)+np.exp(k1*t3)-1.0)
c2 = y0 - (k1*c1/(k2-k1)) + (I0/(k2-k1))*(np.exp(k1*t1)-np.exp(k1*t2)+np.exp(k1*t3)-1.0+(k1/k2)*(1.0-np.exp(k2*t1)+np.exp(k2*t2)-np.exp(k2*t3)))

for i in range(it+1): 
    x_analitica[i] = c1*np.exp(-k1*t[i]) - (I0/k1)*((1.0-(escalon(t[i]-t1)))*(np.exp(k1*(t1-t[i]))-1.0) + (1.0-escalon(t[i]-t2))*(1.0-np.exp(k1*(t2-t[i]))) + (1.0-escalon(t[i]-t3))*(np.exp(k1*(t3-t[i]))-1.0))
    y_analitica[i] = (k1*c1*np.exp(-k1*t[i])/(k2-k1)) + c2*np.exp(-k2*t[i]) - (I0/(k2-k1))*((1.0-escalon(t[i]-t1))*(np.exp(k1*(t1-t[i]))-1.0+(k1/k2)*(1.0-np.exp(k2*(t1-t[i])))) + (1.0-escalon(t[i]-t2))*(1.0-np.exp(k1*(t2-t[i]))+(k1/k2)*(np.exp(k2*(t2-t[i]))-1.0)) + (1.0-escalon(t[i]-t3))*(np.exp(k1*(t3-t[i]))-1.0+(k1/k2)*(1.0-np.exp(k2*(t3-t[i])))))
    
plt.plot(t,x,"--b",label="x_euler")
plt.plot(t,y,"--g",label="y_euler")
plt.plot(t,x_analitica,"b-",label="x_analitica")
plt.plot(t,y_analitica,"g-",label="y_analitica")
plt.legend(loc="upper right")
plt.xlabel('tiempo(min)')
plt.ylabel('masa de medicamento (mg)')
plt.grid()


