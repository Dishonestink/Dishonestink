# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 14:39:29 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

g = 9.81
m1 = 1.0
m2 = 1.0
L1 = 1.0
L2 = 1.0

def dw1dt(teta1,teta2,w1,w2,T1,T2,m1,m2,L1,L2):
    dw1 = ((T2*np.sin(teta2))-(T1*np.sin(teta1)))/(m1*L1*np.cos(teta1))+(w1**2*np.tan(teta1))
    return dw1

def dw2dt(teta1,teta2,w1,w2,T1,T2,m1,m2,L1,L2):
    teta1_2puntos = dw1dt(teta1,teta2,w1,w2,T1,T2,m1,m2,L1,L2)
    dw2 = (m1*L1*(w1**2*np.sin(teta1)-teta1_2puntos*np.cos(teta1))-T2*np.sin(teta2)+m2*L2*(w2**2*np.sin(teta2)))/(m2*L2*np.cos(teta2))
    return dw2

def funciones(t,x,i): 
    #x --> [teta1,teta2,w1,w2]
    teta1 = x[0]
    teta2 = x[1]
    w1 = x[2]
    w2 = x[3]
       
    T2 = m2*g*np.cos(teta2)
    T1 = T2*np.cos(teta2-teta1) + m1*g*np.cos(teta1)
    
    if(i==0):
        f = w1
    elif(i==1):
        f = w2
    elif(i==2):
        f = dw1dt(teta1,teta2,w1,w2,T1,T2,m1,m2,L1,L2)
    elif(i==3):
        f = dw2dt(teta1,teta2,w1,w2,T1,T2,m1,m2,L1,L2)
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

dt = 0.001
tf = 20.0 #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1)
x = np.zeros((it+1,4)) 
#x --> [teta1,teta2,w1,w2] en cualquier fila

x[0,0] = 0.1*np.pi

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    x[i,:] = euler(x[i-1,:],4,dt,t[i-1])

plt.plot(t,x[:,0],"-g",label="teta1")
plt.plot(t,x[:,1],"-b",label="teta2") 

#plt.plot(t,x_rk4[:,2],"-g",label="w1")
#plt.plot(t,x_rk4[:,3],"-b",label="w2")

plt.legend(loc="upper right")
plt.xlabel('Tiempo(s)')
plt.ylabel('Ángulo (rad)')
plt.grid()


