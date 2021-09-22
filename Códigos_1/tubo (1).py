# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 18:47:39 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

#parametros
miu = 0.001 #kg/ms
L = 0.05 #m
R = 0.01 #m
rho = 1000.0 #kg/m3
P_grad = -1.4 #Pa/m
Re = 800.0
U0 = Re*miu/(rho*2*R)

#malla
n = 100 #nodos lado horizontal
m = 20 #nodos lado vertical
nodos = n*m
dx = L/(n-1)
dr = R/(m-1)
mat = np.zeros((nodos,nodos))
U = np.zeros(nodos)
vec_F = np.zeros(nodos) 

eq = 0
for j in range(m):
    for i in range(n):
        print("Ecuación del nodo: ",eq)
        nc = eq #nodo central
        nl = nc - 1 #nodo izquierdo
        nr = nc + 1 #nodo derecho
        nu = nc + n #nodo superior
        nd = nc - n #nodo inferior
        rij = dr*j
        if(i==0): #Nodos tipo 1
            mat[eq,nc] = 1.0
            vec_F[eq] = U0
            print("Nodo de entrada")
        elif(j==0 and i<n-1): #Nodos tipo 2
            mat[eq,nc] = -((2.0/dr**2)+(2.0/dx**2))
            mat[eq,nl] = 1.0/dx**2
            mat[eq,nr] = 1.0/dx**2
            mat[eq,nu] = 2.0/dr**2
            vec_F[eq] = P_grad/miu
            print("Nodo para r = 0")
        elif(j==0 and i==n-1): #Nodos tipo 3
            mat[eq,nc] = -((2.0/dr**2)+(2.0/dx**2))
            mat[eq,nl] = 2.0/dx**2
            mat[eq,nu] = 2.0/dr**2
            vec_F[eq] = P_grad/miu  
            print("Nodo para r = 0, x = L")
        elif(j<m-1 and i<n-1): #Nodos tipo 4 (internos)
            mat[eq,nc] = -((2.0/dr**2)+(2.0/dx**2))
            mat[eq,nl] = 1.0/dx**2
            mat[eq,nr] = 1.0/dx**2
            mat[eq,nu] = (1.0/dr**2) + (1.0/(2.0*rij*dr))
            mat[eq,nd] = (1.0/dr**2) - (1.0/(2.0*rij*dr))
            vec_F[eq] = P_grad/miu
            print("Nodo interno")
        elif(i==n-1 and j<m-1): #Nodos tipo 5
            mat[eq,nc] = -((2.0/dr**2)+(2.0/dx**2))
            mat[eq,nl] = 2.0/dx**2
            mat[eq,nu] = (1.0/dr**2) + (1.0/(2.0*rij*dr))
            mat[eq,nd] = (1.0/dr**2) - (1.0/(2.0*rij*dr))
            vec_F[eq] = P_grad/miu
            print("Nodo de salida")
        else: #Nodos tipo6
            mat[eq,nc] = 1.0
            vec_F[eq] = 0.0  
            print("Nodo en la pared del tubo")
        eq = eq + 1
        
minv = np.linalg.inv(mat)
U = minv@vec_F

r = np.zeros(m)
for i in range(m):
    r[i] = dr*(i)

perfil = np.zeros(m)
distancia = L/3.0
piv = int(distancia/dx)
cont = 0

for i in range(piv,nodos,n):
    perfil[cont] = U[i]
    cont = cont + 1
    
plt.plot(perfil,r,"-k")
plt.ylabel('Distancia radial (m)')
plt.xlabel('Velocidad (m/s)')
plt.grid()    

    
    



