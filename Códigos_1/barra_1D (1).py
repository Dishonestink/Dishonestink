# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:22:49 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

#parametros
L = 1.0 #m
alpha = 1e-4 #m^2/s
teta1 = 300.0 #K --> Temperatura inicial
teta0 = 500.0 #K --> Temperatura ambiente
nodos = 100
dt = 1.0 #s
tf = 1000.0

dx = L/(nodos - 1)
dX = dx/L
dtau = alpha*dt/L**2
gamma = dtau/dX**2
it = int(tf/dt)

#vetores y matrices
T = np.zeros(nodos)
vec_F = np.zeros(nodos)
mat_futura = np.zeros((nodos,nodos))
mat_actual = np.zeros((nodos,nodos))
minv = np.zeros((nodos,nodos))
mat_T = np.zeros((nodos,it+1))

eq = 0
for i in range(nodos):
        print("Ecuación del nodo: ",eq)
        nc = eq #nodo central
        nl = nc - 1 #nodo izquierdo
        nr = nc + 1 #nodo derecho

        if((i==0)): #Nodos tipo 10
            mat_futura[eq,nc] = 1.0
            
            print("Nodo en x=0")
            
        elif((i<nodos-1)):
            mat_futura[eq,nl] = -gamma
            mat_futura[eq,nc] = 2.0*gamma + 1.0
            mat_futura[eq,nr] = -gamma
            
            mat_actual[eq,nc] = 1.0
            
            print("Nodo interno")
            
        else:
            mat_futura[eq,nl] = -2.0*gamma
            mat_futura[eq,nc] = 2.0*gamma + 1.0
            
            mat_actual[eq,nc] = 1.0
            
            print("Nodo en x=L")
            
        eq = eq + 1
        
        
minv = np.linalg.inv(mat_futura)

T[:] = 1.0 #K condicion inicial
mat_T[:,0] = 1.0

for j in range(1,it+1):
    print("iteracion: ",j)
    T = minv@((mat_actual@T)+vec_F)
    mat_T[:,j] = (teta1 - teta0)*T + teta0    

    
x = np.zeros(nodos)
for i in range(nodos):
    x[i] = dx*(i)    
  
plt.plot(x,mat_T[:,1],"-k",label="it = 1")
plt.plot(x,mat_T[:,100],"-g",label="it = 100")
plt.plot(x,mat_T[:,200],"-r",label="it = 200")
plt.plot(x,mat_T[:,500],"-b",label="it = 300")
plt.plot(x,mat_T[:,1000],"-y",label="it = 1000")

plt.legend(loc="upper right")
plt.xlabel('Distancia (m)')
plt.ylabel('Temperatura (K)')
plt.grid()





    
            
            
            






