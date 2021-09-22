# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 08:35:45 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

T = 1500.0 #K

def velocidad_reaccion():
    Ru = 8.315   #J/molK
    kf = np.zeros(2)
    kr = np.zeros(2)
    with open('kf.dat', 'r') as file1:
        for line in range(2):
            A,b,Ea = [float(x) for x in next(file1).split()] # read first line
            Ea = 4.184*Ea   #J/mol
            kf[line] = A*T**b*np.exp(-Ea/(Ru*T))
    with open('kr.dat', 'r') as file2:
        for line in range(2):
            A,b,Ea = [float(x) for x in next(file2).split()] # read first line
            Ea = 4.184*Ea   #J/mol
            kr[line] = A*T**b*np.exp(-Ea/(Ru*T))
    return kf,kr
           
def funciones(t,x,i): 
    #x --> [x0,x1,x2,x3,x4]
    #x --> [h ,h2,o ,o2,oh]
    [kf,kr] = velocidad_reaccion()
    if(i==0):
        f = -kf[0]*x[0]*x[3] + kr[0]*x[2]*x[4] + kf[1]*x[1]*x[2] - kr[1]*x[0]*x[4]
    elif(i==1):
        f = -kf[1]*x[1]*x[2] + kr[1]*x[0]*x[4]
    elif(i==2):
        f = kf[0]*x[0]*x[3] - kr[0]*x[2]*x[4] - kf[1]*x[1]*x[2] + kr[1]*x[0]*x[4]
    elif(i==3):
        f = -kf[0]*x[0]*x[3] + kr[0]*x[2]*x[4]
    else:
        f = kf[0]*x[0]*x[3] - kr[0]*x[2]*x[4] + kf[1]*x[1]*x[2] - kr[1]*x[0]*x[4]
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

def rk4(fi_in,n,dt,t):
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)
    fi_out = np.zeros(n)
    for i in range(n):
        k1[i] = funciones(t,fi_in,i)
    for i in range(n):
        k2[i] = funciones(t+0.5*dt,fi_in+0.5*k1*dt,i)
    for i in range(n):
        k3[i] = funciones(t+0.5*dt,fi_in+0.5*k2*dt,i)
    for i in range(n):
        k4[i] = funciones(t+dt,fi_in+k3*dt,i)
    for i in range(n):
        fi_out[i] = fi_in[i] + (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*dt/6.0
    return fi_out

dt = 1.0e-8    #paso de tiempo en minutos
tf = 1.0e-6   #Tiempo de simulación 
it = int(tf/dt)

t = np.zeros(it+1) 
x_ef = np.zeros((it+1,5))
x_rk4 = np.zeros((it+1,5))

x_ef[0,0] = 1.0e-5 #mol/cm^3
x_ef[0,3] = 1.0e-6 #mol/cm^3

x_rk4[0,0] = 1.0e-5 #mol/cm^3
x_rk4[0,3] = 1.0e-6 #mol/cm^3

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    x_ef[i,:] = euler(x_ef[i-1,:],5,dt,t[i-1])
    x_rk4[i,:] = rk4(x_rk4[i-1,:],5,dt,t[i-1])
    
plt.plot(t,x_ef[:,0],"-k",label="h")
plt.plot(t,x_ef[:,1],"-b",label="h2")
plt.plot(t,x_ef[:,2],"-r",label="o")
plt.plot(t,x_ef[:,3],"-g",label="o2")
plt.plot(t,x_ef[:,4],"-y",label="oh")

plt.plot(t,x_rk4[:,0],"--k",label="h")
plt.plot(t,x_rk4[:,1],"--b",label="h2")
plt.plot(t,x_rk4[:,2],"--r",label="o")
plt.plot(t,x_rk4[:,3],"--g",label="o2")
plt.plot(t,x_rk4[:,4],"--y",label="oh")

plt.legend(loc="upper right")
plt.xlabel('Tiempo(s)')
plt.ylabel('Concentración (mol/cm^3)')
plt.grid()
    




