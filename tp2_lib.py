import numpy as np
import scipy as sp
import scipy.sparse as spsp    ## Sparse matrix 
import matplotlib as mpl        
import matplotlib.pyplot as plt
from math import * 

def indk(i,j, K):
    return i+j*(K+1)

def vit(x,y,V0,L):
    V1 = -V0*np.sin(np.pi *x/L)*np.cos(np.pi*y/L)
    V2 = V0*np.sin(np.pi *y/L)*np.cos(np.pi*x/L)
    return np.array([V1,V2])

def source(x,y,L):
    xl = x/L 
    yl = y/L 
    return 256*(xl**2)*(yl**2)*((1-xl)**2)*((1-yl)**2)

def MatA(L, V0, D, h, dt, K, indk, vit):

    K2 = (K+1)**2
    A = spsp.lil_matrix((K2, K2))   ## Declaration of A as a sparse matrix 
## Loop on the internal nodes
    for i in range(1,K):
        for j in  range(1,K):
            v1, v2 = vit(L*i/K,L*j/K,V0,L)
            k = indk(i,j,K)
            ke = indk(i+1,j,K)
            ko = indk(i-1,j,K)
            kn = indk(i,j+1,K)
            ks = indk(i,j-1,K)
            ##to be completed 
            A[k,k] = 1 - 4*dt*D/(h**2)
            A[k,ke] = -dt*v1/(2*h)+dt*D/(h**2)       #i+1,j 
            A[k,ko] = dt*v1/(2*h)+dt*D/(h**2)       #i-1,j
            A[k,kn] = -dt*v2/(2*h)+dt*D/(h**2)
            A[k,ks] = dt*v2/(2*h)+dt*D/(h**2)
            
            
 ## Loop on the nodes located in x = 0 (Neumann Boundary Condition) 
    for j in range(1,K):
        i = 0 
        k = indk(i,j,K)
        ke = indk(i+1,j,K)
        kn = indk(i,j+1,K)
        ks = indk(i,j-1,K)
        v2 = vit(L*i/K,L*j/K,V0,L)[1]
        A[k,k] = 1 - 3*dt*D/(h**2)
        A[k,ke] = dt*D/(h**2)
        A[k,kn] = -dt*v2/(2*h)+dt*D/(h**2)
        A[k,ks] = dt*v2/(2*h)+dt*D/(h**2)
 ## Loop on the nodes located in x = L (Neumann Boundary Condition) 
    for j in range(1,K):
        i = K
        k = indk(i,j,K)
        ko = indk(i-1,j,K)
        kn = indk(i,j+1,K)
        ks = indk(i,j-1,K)
        v2 = vit(L*i/K,L*j/K,V0,L)[1]
        A[k,k] = 1 - 3*dt*D/(h**2)
        A[k,ko] = dt*D/(h**2)
        A[k,kn] = -dt*v2/(2*h)+dt*D/(h**2)
        A[k,ks] = dt*v2/(2*h)+dt*D/(h**2)
    
    return A


def secmb(K, dt, h, L, Tbas, Thaut, indk, source):
    S = np.zeros((K+1)**2)
    for i in range(K+1):
        j = 0 
        k = indk(i,j,K)
        S[k] = Tbas 
    for i in range(K+1):
        j = K 
        k = indk(i,j,K)
        S[k] = Thaut 
    for i in range(1,K):
        for j in  range(1,K):
            k = indk(i,j,K)
            S[k] = dt * source(L*i/K, L*j/K,L)
    return S 
           
    
def matT(T, K, indk):
    T2D = np.zeros((K+1,K+1))
    for i in range(K+1):
        for j in range(K+1):
            T2D[i,j] = T[indk(i,j,K)]
    return T2D


# def evolution(X, Y, Theta, K, dt, L, D, V0, Tbas, Thaut, source, vit):

#     Vit = np.zeros([K,2])
#     S = np.zeros(K)
#     alpha = np.random.randn(K,2)
    
#     for k in range(0, K):
#         V = vit(X[k],Y[k],V0,L)
#         Vit[k,0]= V[0]
#         Vit[k,1]= V[1]
#         S[k] = source(X[k],Y[k],L)

#     ## to be completed
    
    
#     return X,Y,Theta


# def tmoy(X, Y, Theta, K, M, L):

#     Tm = np.zeros([M,M])    ## mean temperature in a cell 
#     nbp = np.zeros([M,M])   ## number of particles in a cell 
    
#     for k in range(0, K):
#         ## to be completed 

#     return Tm
    

# def posinit(K, L):

#     ## to be completed
    
#     return X,Y
    
    


    
    
