# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 20:44:51 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

#paramtetros
T1 = 400.0 #K
T2 = 350.0 #K
T_inf = 300 #k
h = 500.0 #W/m2K
k = 205 #W/mK
L = 1.0 #m

#malla
n = 10 #nodos lado horizontal
m = 10 #nodos lado vertical
nodos = n*m
dx = L/(n-1)
dy = L/(m-1)

#memoria
mat = np.zeros((nodos,nodos))
T = np.zeros(nodos)
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
        if(i==0): #Nodos tipo 1
            mat[eq,nc] = 1.0
            vec_F[eq] = T1
            print("Nodo pared lateral izquierda")
        elif(j==0 and i<n-1): #Nodos tipo 2
            mat[eq,nc] = -((2.0/dx**2)+(2.0/dy**2))
            mat[eq,nl] = 1.0/dx**2
            mat[eq,nr] = 1.0/dx**2
            mat[eq,nu] = 2.0/dy**2
            print("Nodo pared inferior")
        elif(i==n-1): #Nodos tipo 3
            mat[eq,nc] = 1.0
            vec_F[eq] = T2
            print("Nodo pared lateral derecha")
        elif(j<m-1): #Nodos tipo 4 (internos)
            mat[eq,nc] = -((2.0/dx**2)+(2.0/dy**2))
            mat[eq,nl] = 1.0/dx**2
            mat[eq,nr] = 1.0/dx**2
            mat[eq,nu] = 1.0/dy**2
            mat[eq,nd] = 1.0/dy**2
            print("Nodo interno")
        else: #Nodos tipo 5
            mat[eq,nc] = -((2.0/dx**2)+(2.0/dy**2)+(2.0*h/(k*dy)))
            mat[eq,nl] = 1.0/dx**2
            mat[eq,nr] = 1.0/dx**2
            mat[eq,nd] = 2.0/dy**2
            vec_F[eq] = -2.0*h*T_inf/(k*dy)
            print("Nodo pared superior")
        eq = eq + 1
        
minv = np.linalg.inv(mat)
T = minv@vec_F

x = np.zeros(n)
for i in range(n):
    x[i] = dx*(i)

plt.plot(x,T[0:n],"-k",label="j=0")
plt.plot(x,T[n:2*n],"-b",label="j=1")
plt.plot(x,T[2*n:3*n],"-g",label="j=2")
plt.plot(x,T[3*n:4*n],"-y",label="j=3")
plt.plot(x,T[4*n:5*n],"-r",label="j=4")
plt.plot(x,T[5*n:6*n],"-k",label="j=5")
plt.plot(x,T[6*n:7*n],"-b",label="j=6")
plt.plot(x,T[7*n:8*n],"-g",label="j=7")
plt.plot(x,T[8*n:9*n],"-y",label="j=8")
plt.plot(x,T[9*n:10*n],"-r",label="j=9")
plt.legend(loc="upper right")
plt.xlabel('Distancia (m)')
plt.ylabel('Temperautura (K)')
plt.grid()



    


