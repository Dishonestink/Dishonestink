# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 20:42:59 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

def fy(y,v,t):
    g = 9.81
    y = y + v*t - 0.5*g*t**2.0
    return y

def fv(y1,y2,v):
    g = 9.81
    v = (v**2.0 -2.0*g*(y2-y1))**0.5
    if(y2>y1):  
        v = v
    else:
        v = -v
    return v

def arrastre(y1,y2):
    import random
    val = 0
    i = 0    
   
    while(val==0):
        x = random.random() #Le asigna a x un numero entre 0 y 1
        if(x<=0.15):
            val = 1
        else:
            i = i + 1
    if(y2>y1):  
        y2 = y2*(1.0-0.1*x)
        #y2 = y2
    else:
        y2 = y2*(1.0+0.1*x)
        #y2 = y2
    return y2,i

v0 = 10.0
dt = 0.01
n = 300
y = np.zeros(n)
t = np.zeros(n)
v = np.zeros(n)

y[0] = 0.0
t[0] = 0.0
v[0] = v0

for i in range(1,n):
    t[i] = t[i-1] + dt
    y[i] = fy(y[i-1],v[i-1],dt)
    [y[i],cont] = arrastre(y[i-1],y[i])
    v[i] = fv(y[i-1],y[i],v[i-1])
    print("Iteración:",i)
    print("posicion calculada:",y[i])
    print("velocidad calculada:",v[i])
    print("iteraciones para calcular x:" , cont)
    if(y[i]<=0.0):
        import sys
        plt.plot(t[0:i],y[0:i])
        sys.exit("La partícula ha regresado a su posición inicial")
        
#pueden usar break en vez de sys.exit
