# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:11:53 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D




#parametros
h1 = 13000 #W/m2K
T_inf1 = 273.15 + 50.0 #k
h2 = 1000.0 #W/m2K
T_inf2 = 273.15 + 20.0 #k
k = 70.0 #W/mK
Cp = 300 #J/kgK
rho = 14950.0 #kg/m³
q = 720000.0 #W/m²
x1 = 12.0e-3 #m 
x2 = 40.0e-3 #m 
x_total = 50.0e-3 #m 
y_total = 20.0e-3 #m 


dx = 1.0e-3
dy = dx
alpha = k/(rho*Cp)
Bi1 = h1*dy/k
Bi2 = h2*dy/k
shi = 0.5
dt = 0.2
Fo = alpha*dt/dx**2
#Fo = 0.24
#dt = Fo*dx**2/alpha
tf = 100.0
it = int(tf/dt)


#nodos
nx = int(x_total/dx) + 1
ny = int(y_total/dy) + 1
nx1 = int(x1/dx) + 1
nx2 = int(x2/dx) + 1
nodos = nx*ny


T = np.zeros(nodos)
vec_F = np.zeros(nodos)
mat_futura = np.zeros((nodos,nodos))
mat_actual = np.zeros((nodos,nodos))
minv = np.zeros((nodos,nodos))
mat_T = np.zeros((nodos,it+1))

eq = 0
for j in range(ny):
    for i in range(nx):
        print("Ecuación del nodo: ",eq)
        nc = eq #nodo central
        nl = nc - 1 #nodo izquierdo
        nr = nc + 1 #nodo derecho
        nu = nc + nx #nodo superior
        nd = nc - nx #nodo inferior

        if((i==0) and (j==0)): #Nodos tipo 10
            mat_futura[eq,nc] = 1.0+(4.0-2.0*Bi2)*shi*Fo
            mat_futura[eq,nr] = -2.0*shi*Fo
            mat_futura[eq,nu] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-(4.0-2.0*Bi2)*(1.0-shi)*Fo
            mat_actual[eq,nr] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nu] = 2.0*(1.0-shi)*Fo

            vec_F[eq] = -2.0*Bi2*Fo*T_inf2
            print("Nodos tipo 10")
        elif((j==0)and(i<=nx2-1)): #Nodos tipo 9
            mat_futura[eq,nc] = 1.0+(4.0-2.0*Bi2)*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nu] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-(4.0-2.0*Bi2)*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nu] = 2.0*(1.0-shi)*Fo
            
            vec_F[eq] = -2.0*Bi2*Fo*T_inf2
            print("Nodos tipo 9")
        elif((j==0)and(i<nx-1)): #Nodos tipo 5
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nu] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nu] = 2.0*(1.0-shi)*Fo
            print("Nodos tipo 5")
        elif((j==0)and(i==nx-1)): #Nodos tipo 7
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -2.0*shi*Fo
            mat_futura[eq,nu] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nu] = 2.0*(1.0-shi)*Fo
            print("Nodos tipo 7")
        elif((i==0)and(j<ny-1)): #Nodos tipo 2
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nr] = -2.0*shi*Fo
            mat_futura[eq,nu] = -shi*Fo
            mat_futura[eq,nd] = -shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nr] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nu] = (1.0-shi)*Fo
            mat_actual[eq,nd] = (1.0-shi)*Fo
            print("Nodos tipo 2")
        elif((i<nx-1)and(j<ny-1)): #Nodos tipo 1
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nu] = -shi*Fo
            mat_futura[eq,nd] = -shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nu] = (1.0-shi)*Fo
            mat_actual[eq,nd] = (1.0-shi)*Fo
            print("Nodos tipo 1")
        elif((i==nx-1)and(j<ny-1)): #Nodos tipo 3
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -2.0*shi*Fo
            mat_futura[eq,nu] = -shi*Fo
            mat_futura[eq,nd] = -shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nu] = (1.0-shi)*Fo
            mat_actual[eq,nd] = (1.0-shi)*Fo
            print("Nodos tipo 3")
        elif((i==0)and(j==ny-1)): #Nodos tipo 12
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nr] = -2.0*shi*Fo
            mat_futura[eq,nd] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nr] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nd] = 2.0*(1.0-shi)*Fo
            
            vec_F[eq] = 2.0*dy*Fo*q/k
            print("Nodos tipo 12")
        elif((i<=nx1-1)and(j==ny-1)): #Nodos tipo 11
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nd] = -2.0*shi*Fo
            
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nd] = 2.0*(1.0-shi)*Fo
            
            vec_F[eq] = 2.0*dy*Fo*q/k
            print("Nodos tipo 11")
        elif((i<=nx2-1)and(j==ny-1)): #Nodos tipo 8
            mat_futura[eq,nc] = 1.0+(4.0+2.0*Bi1)*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nd] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-(4.0+2.0*Bi1)*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nd] = 2.0*(1.0-shi)*Fo
            
            vec_F[eq] = 2.0*Bi1*Fo*T_inf1
            print("Nodos tipo 8")
        elif((i<nx-1)and(j==ny-1)): #Nodos tipo 4
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -shi*Fo
            mat_futura[eq,nr] = -shi*Fo
            mat_futura[eq,nd] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = (1.0-shi)*Fo
            mat_actual[eq,nr] = (1.0-shi)*Fo
            mat_actual[eq,nd] = 2.0*(1.0-shi)*Fo
            print("Nodos tipo 4")
        else: #Nodos tipo 6
            mat_futura[eq,nc] = 1.0+4.0*shi*Fo
            mat_futura[eq,nl] = -2.0*shi*Fo
            mat_futura[eq,nd] = -2.0*shi*Fo
            
            mat_actual[eq,nc] = 1.0-4.0*(1.0-shi)*Fo
            mat_actual[eq,nl] = 2.0*(1.0-shi)*Fo
            mat_actual[eq,nd] = 2.0*(1.0-shi)*Fo
            print("Nodos tipo 6")
        eq = eq + 1

        
minv = np.linalg.inv(mat_futura)

T[:] = 290.0 #K condicion inicial
mat_T[:,0] = 290.0

for j in range(1,it+1):
    print("iteracion: ",j)
    T = minv@((mat_actual@T)+vec_F)
    mat_T[:,j] = T
    
x = np.zeros(nx)
for i in range(nx):
    x[i] = dx*(i)
    
    
    

    

plt.plot(x,mat_T[nodos-nx:nodos,int(0.1*it)],"-k",label="0.1tf")    
plt.plot(x,mat_T[nodos-nx:nodos,int(0.3*it)],"-g",label="0.3tf")    
plt.plot(x,mat_T[nodos-nx:nodos,int(0.5*it)],"-r",label="0.5tf")    
plt.plot(x,mat_T[nodos-nx:nodos,int(0.8*it)],"-b",label="0.8tf")
plt.plot(x,mat_T[nodos-nx:nodos,it],"-y",label="tf")
# print(mat_T[nodos-nx:nodos,it])
# print(nodos-nx,nodos)

# tiempo = np.zeros(it+1)
# for i in range(it+1):
#     tiempo[i] = i*dt

# tie_m = np.meshgrid(tiempo)

# fig = plt.figure()
# ax = Axes3D(fig)
# xl = np.linspace(0,x_total,50)
# yl = np.linspace(0,tf,50)
# xm,ym = np.meshgrid(x,tiempo)
    

# ax.contourf(tiempo,x,mat_T[nodos-nx:nodos,:])
# plt.show()
plt.legend(loc="upper right")
plt.xlabel('Distancia (m)')
plt.ylabel('Temperatura (K)')
plt.grid()










