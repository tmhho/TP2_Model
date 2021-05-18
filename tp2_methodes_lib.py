import numpy as np
import scipy as sp
import scipy.sparse as spsp  
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import * 
import tp2_lib as tp


def methodeDF(K,tmax,L):
#######################################################################
##                 Finite Difference Method                          ##
#######################################################################

    ## Continuous Problem Data
    V0 = 1.0
    D = 0.2
    Tbas = 0.   ## (BC en y = 0)
    Thaut = 1.  ## (BC en y = L)


    ## Numerical parameters
    K2 = (K+1)**2 
    h = L/K
    dt = 0.8 * h*h / (4*D)
    
    # Initialisation of T
    T = np.zeros(K2)
    print("T shape ", T.shape)

    # Assembling of the Right Hand Side
    S = tp.secmb(K, dt, h, L, Tbas, Thaut, tp.indk, tp.source)
    print("S shape ", T.shape)
    # Assembling of the Scheme Matrix 
    A = tp.MatA(L, V0, D, h, dt, K, tp.indk, tp.vit)
    print("A shape ", A.shape)
    # Time Loop 
    time = 0.0
    while time < tmax :
        ## to be completed
        # print(T.shape)
        # print(A.shape)
        # print(S.shape)
        T = A.dot(T) + S 
        # print(T)
        time += dt 
    return T



# def methodeMC(K,M,tmax,L):
# #######################################################################
# ##                       Monte Carlo Method                          ##
# #######################################################################

#     ## Continuous Problem Data
#     V0 = 1.0
#     D = 0.2
#     Tbas = 0.   ## (BC en y = 0)
#     Thaut = 1.  ## (BC en y = L)


#     ## Numerical parameters
#     eps = L/M
#     dt = 0.25 * eps * eps / D 

#     # Initialisation of Theta and X,Y
#     X,Y = tp.posinit(K,L)
#     Theta = np.zeros(K)

#     # Time Loop 
#     time = 0.0
#     while time < tmax :
#         ## to be completed

        
#     # Compute Temperaure field
#     T = tp.tmoy(X, Y, Theta, K, M, L)
    
#     return T
