# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:03:16 2018

@author: niemel34
"""
import scipy as sp
from matplotlib import pyplot as plt
from scipy import optimize

#-------------------------------------------------------------------

def Re_func(nu, V, d):
    
    return V*d/nu

#-------------------------------------------------------------------

def f_func(Re, Pr, eps, d):
    
    if Re <= 4000:
        
        raise ValueError('f_func: Virtaus laminaari, korrelaatio ei pade')
    
    elif Re > 4000 and Re < 5e6:
        
        def to_zero(f):
            return 1
        
        return optimize.fsolve(to_zero, 1e-2)[0]

    else:
        
        raise ValueError('f_func: Re-luku liian suuri, korrelaatio ei pade')

#-------------------------------------------------------------------
        
def Nu_func(Re, Pr, f):
        
    pass
        
#-------------------------------------------------------------------
        
def h_func(Nu, d, k):
        
    pass

#-------------------------------------------------------------------

def T_out_func(T_w, T_in, h, d, m, L, cp):
        
    pass

#-------------------------------------------------------------------

def dp_func(f, L, d, rho, V):
        
    pass

#-------------------------------------------------------------------

def main(): 

    # Aineominaisuudet, Vesi 41 Celciusta
    cp = 0 # Ominaislämpökapasiteetti, J/kgK
    rho = 0 # Tiheys, kg/m3
    Pr = 0 # Prandtlin luku, -
    k = 0 # Lämmönjohtavuus, W/mK
    nu = 1 # Kinemaattinen viskositeetti, m2/s
    
    # Vakiot
    L = 0 # putken pituus metreinä
    eps = 0 # putken karheus metreinä
    T_in = 0 # Sisääntulolämpötila, K
    T_w = 0 # Putken sisäpinnan lämpötila, K
    m = 0 # Massavirta, kg/s
    d = sp.array([0, 0]) # putken halkaisija metreinä
    V = sp.array([0, 0])
    
    #Laskenta
    
    n_d = len(d) # Eri halkaisijoiden lukumäärä
    Re = sp.zeros(n_d)
    f = sp.zeros(n_d)
    Nu = sp.zeros(n_d)
    h = sp.zeros(n_d)
    T_out = sp.zeros(n_d)
    dp = sp.zeros(n_d)
    
    for i in range(n_d):
        Re[i] = Re_func(nu, V[i], d[i])



    print('Re = ', dp)

    
#    plt.plot(d*1000, T_out-273)
#    plt.xlabel('d (mm)')
#    plt.ylabel('Ulostulolampotila (C)')
    
#    plt.plot(d*1000, Re)
#    plt.xlabel('d (mm)')
#    plt.ylabel('Re (-)')

    
if __name__ == "__main__":
    main()