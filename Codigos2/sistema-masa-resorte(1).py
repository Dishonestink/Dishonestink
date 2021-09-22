# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 18:15:30 2020

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

m1 = 1.0
m2 = 0.5
m3 = 2.0
k1 = 500.0
k2 = 100.0
k3 = 200.0
k4 = 1000.0
b1 = 1.0
b2 = 3.0

def funciones(t,x,i):
    #x --> [x1,v1,x2,v2,x3,v3]
    x1 = x[0]
    v1 = x[1]
    x2 = x[2]
    v2 = x[3]
    x3 = x[4]
    v3 = x[5]
       
    if(i==0):
        f = v1
    elif(i==1):
        f = (-k1*x1+k2*(x2-x1)+b1*(v2-v1))/m1
    elif(i==2):
        f = v2
    elif(i==3):
        f = (-k2*(x2-x1)-b1*(v2-v1)+k3*(x3-x2)+b2*(v3-v2))/m2
    elif(i==4):
        f = v3
    elif(i==5):
        f = (-k3*(x3-x2)-b2*(v3-v2)-k4*x3)/m3
    return f

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

A = np.array([[0,1,0,0,0,0],
              [-(k1+k2)/m1,-b1/m1,k2/m1,b1/m1,0,0],
              [0,0,0,1,0,0],
              [k2/m2,b1/m2,-(k2+k3)/m2,-(b1+b2)/m2,k3/m2,b2/m2],
              [0,0,0,0,0,1],
              [0,0,k3/m3,b2/m3,-(k3+k4)/m3,-b2/m3]])

eigen_val,eigen_vec = np.linalg.eig(A)

#condiciones iniciales
x1_0 = 0.1
v1_0 = 0.0
x2_0 = 0.0
v2_0 = 0.0
x3_0 = 0.0
v3_0 = 0.0

B = [x1_0,v1_0,x2_0,v2_0,x3_0,v3_0]

minv = np.linalg.inv(eigen_vec)
C = minv@B

dt = 0.05
tf = 5.0 #Tiempo de simulación 
it = int(tf/dt)
t = np.zeros(it+1)
x = np.zeros((it+1,6))
x_rk4 = np.zeros((it+1,6))
err = np.zeros((it+1,6))
#x --> [x1,v1,x2,v2,x3,v3]

x[0,:] = [x1_0,v1_0,x2_0,v2_0,x3_0,v3_0]
x_rk4[0,:] = [x1_0,v1_0,x2_0,v2_0,x3_0,v3_0]

for i in range(1,it+1):
    t[i] = t[i-1] + dt
    for k in range(6):
        for j in range(6):
            alpha = np.real(eigen_val[j])
            beta = np.imag(eigen_val[j])
            a = np.real(C[j]*eigen_vec[k,j])
            b = np.imag(C[j]*eigen_vec[k,j])
            x[i,k] = x[i,k] + np.exp(alpha*t[i])*(a*np.cos(beta*t[i])-b*np.sin(beta*t[i]))
    
    x_rk4[i,:] = rk4(x_rk4[i-1,:],6,dt,t[i-1])
    
ac = np.zeros(it+1)
for i in range(1,it+1):
    ac[i] = x[i,1]/t[i]
    
print(ac)  
print(len(ac))
err = abs(x-x_rk4)

# plt.plot(t,x[:,0],"-b",label="x1_analítica")
# plt.plot(t,x[:,2],"-g",label="x2_analítica")
# plt.plot(t,x[:,4],"-y",label="x3_analítica")
# plt.legend(loc="upper right")
# plt.xlabel('Tiempo(s)')
# plt.ylabel('Posición de las masas (m)')
# plt.grid()

plt.plot(t,x[:,0],"-k",label="x1_analítica")
plt.plot(t,x_rk4[:,0],"--k",label="x1_rk4")
plt.legend(loc="upper right")
plt.xlabel('Tiempo(s)')
plt.ylabel('Posición de las masas (m)')
plt.grid()

# plt.plot(t,err[:,0],"-k",label="Error en x1")
# plt.legend(loc="upper right")
# plt.xlabel('Tiempo(s)')
# plt.ylabel('Error absoluto (m)')
# plt.grid()



            
    

