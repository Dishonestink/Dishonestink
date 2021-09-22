import numpy as np

P = 1.0 
T = 3000.0
Ru = 8.315 #J/molK
fi = 1.5
n = 8

def equilibrio(x,i):
    #i representa el numero de la reaccion
    Gx = np.array([46182.0,0.0,54554.0,0.0,-4245.0,-77326.0,479933.0,381915.0])
    
    vr = np.array([[0,1,0,0,0,0,0,0],
	       	       [0,0,0,1,0,0,0,0],
	               [0,0,0,0,0,1,0,0],
                   [1,0,0,1,0,0,0,0],
	               [0,1,0,1,0,0,0,0],
	               [0,0,1,0,0,1,0,0]])
    
    vp = np.array([[2,0,0,0,0,0,0,0],
	               [0,0,2,0,0,0,0,0],
	               [0,0.5,0,0,1,0,0,0],
	               [0,0,0,0,0,0,1,0],
                   [0,0,0,0,0,0,0,1],
	               [0,0,0,0,2,0,0,0]])

    a = 1.0
    b = 1.0
    dG_r = 0.0
    dG_p = 0.0
    for j in range(n):
        a = a*x[j]**vp[i,j]
        b = b*x[j]**vr[i,j]
        dG_r = dG_r + vr[i,j]*Gx[j]
        dG_p = dG_p + vp[i,j]*Gx[j]        
        
    dG = dG_p - dG_r
    f_reac = ((a/b)*P**(np.sum(vp[i,:])/np.sum(vr[i,:])) 
             - np.exp(-dG/(Ru*T)))
    return f_reac

def funciones(x,i):
    #i representa el numero de la funcion
    if(i<6):
        f = equilibrio(x,i)
    elif(i==6):
        f = (x[0] + 2*x[1] -2*fi*x[2] - 4*fi*x[3] + (1 - 2*fi)*x[4]
            +(2 - 2*fi)*x[5] + (1 - 4*fi)*x[6] + (2 - 4*fi)*x[7])
    else:
        f = -1.0 + np.sum(x)
    return f

def derivada(x,i,j):
#i es el numero de la funcion a derivar numericamente
#j la especie con respecto a la cual derivar
    if(abs(x[j])>1.0):
        eps = 1.0e-5*abs(x[j])
    else:
        eps = 1.0e-5
    
    x_eps = np.zeros(8) 
    x_eps[:] = x
    x_eps[j] = x_eps[j] + eps
    d_num = (funciones(x_eps,i)-funciones(x,i))/eps
    return d_num

x = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

k = 0
cont = 0
Jac = np.zeros((n,n))
f = np.zeros(n)
d = np.zeros(n)

while(cont<n):
    k = k + 1
    for i in range(n):
        for j in range(n):
            Jac[i,j] = derivada(x,i,j)
            
        f[i] = -funciones(x,i)
    
    minv = np.linalg.inv(Jac)
    d = minv@f
    
    xf = x + d
    norm_old = 0.0
    norm_new = 0.0
    cont = 0
    
    for i in range(n):
        norm_old = norm_old + abs(funciones(x,i))
        norm_new = norm_new + abs(funciones(xf,i)) 
        val = abs(d[i]/x[i])
        if(val<1.0e-7):
            cont = cont + 1
            
    if(norm_new>norm_old):
        print("La nueva norma resultó ser mayor en iteracion:",k)
        xf = x + d/5.0  
        
    x[:] = xf
    
print("El método convergio a ",k,"iteraciones")
print("x_h:",x[0])
print("x_h2:",x[1])
print("x_o:",x[2])
print("x_o2:",x[3])
print("x_oh:",x[4])
print("x_h2o:",x[5])
print("x_ho2:",x[6])
print("x_h2o2:",x[7])




