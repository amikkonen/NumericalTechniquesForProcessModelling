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
        
    elif Pr > 2000 or Pr < 0.5:
        
        raise ValueError('f_func: Prandtlin luku ei patevyysalueella (0.5-2000)')
    
    elif Re > 4000 and Re < 5e6:
        
        def to_zero(f):
            return f**(-0.5) + 2*sp.log10( eps/(3.7*d) + 2.51/(Re*f**0.5) )
        
        return optimize.fsolve(to_zero, 1e-2)[0]

    else:
        
        raise ValueError('f_func: Re-luku liian suuri, korrelaatio ei pade')

#-------------------------------------------------------------------
        
def Nu_func(Re, Pr, f):
        
    if Re <= 3000:
        
        raise ValueError('Nu_func: Virtaus laminaari, korrelaatio ei pade')
    
    elif Re > 3000 and Re <= 5e6:
        
        return ((f/8)*(Re-1000)*Pr) / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
    
    else:
        
        raise ValueError('Nu_func: Re-luku liian suuri, korrelaatio ei pade')
        
#-------------------------------------------------------------------
        
def h_func(Nu, d, k):
        
    return Nu*k/d

#-------------------------------------------------------------------

def T_out_func(T_w, T_in, h, d, m, L, cp):
        
    return T_w - (T_w-T_in) * sp.exp( -h*2*sp.pi*(d/2)*L / (m*cp) )

#-------------------------------------------------------------------

def dp_func(f, L, d, rho, V):
        
    return f*(L/d)*0.5*rho*V**2

#-------------------------------------------------------------------

def main(): 

    # Aineominaisuudet, Vesi 40 Celciusta
    cp  = 4180 # Ominaislämpökapasiteetti, J/kgK
    rho = 992 # Tiheys, kg/m3
    Pr  = 4.16 # Prandtlin luku, -
    k   = 0.634 # Lämmönjohtavuus, W/mK
    nu  = 6.58e-7 # Kinemaattinen viskositeetti, m2/s
    
    # Vakiot
    L    = 6.3 # putken pituus metreinä
    eps  = 1.3e-6 # putken karheus metreinä
    T_in = 30 + 273 # Sisääntulolämpötila, K
    T_w  = 48 + 273 # Putken sisäpinnan lämpötila, K
    m    = 0.11 # Massavirta, kg/s
    d    = sp.array([0.012, 0.013, 0.014, 0.016,
                     0.020, 0.023, 0.026, 0.033]) # putken halkaisijat metreinä
    V    = m/(rho*sp.pi*(d/2)**2)
#    print('V = ', V)
    
    #LASKENTA
    n_d = len(d) # Eri halkaisijoiden lukumäärä
    
    # Alustetaan arrayt ensin nollilla, ja lasketaan sitten for-silmukassa
    # muuttujien arvot kaikille eri putkien halkaisjoille
    Re    = sp.zeros(n_d) 
    f     = sp.zeros(n_d)
    Nu    = sp.zeros(n_d)
    h     = sp.zeros(n_d)
    T_out = sp.zeros(n_d)
    dp    = sp.zeros(n_d)
    
    for i in range(n_d):
        Re[i]    = Re_func(nu, V[i], d[i])
        f[i]     = f_func(Re[i], Pr, eps, d[i])
        dp[i]    = dp_func(f[i], L, d[i], rho, V[i])
        Nu[i]    = Nu_func(Re[i], Pr, f[i])
        h[i]     = h_func(Nu[i], d[i], k)
        T_out[i] = T_out_func(T_w, T_in, h[i], d[i], m, L, cp)

    # Tarkistetaan tulosten järkevyys tulostamalla arvot
    print('V =     ', V)
    print('Re =    ', Re)
    print('f =     ', f)
    print('dp =    ', dp)
    print('Nu =    ', Nu)
    print('h =     ', h)
    print('T_out = ', T_out)
    
    # Piirretään ulostulolämpötila ja painehäviö halkaisijan funktiona ja
    # tallennetaan kuvat PDF-tiedostoiksi (tallennus kommentoituna)
    plt.figure(0)
    plt.plot(d*1000, T_out-273, 'k-d')
    plt.xlabel('d (mm)')
    plt.ylabel('Ulostulolampotila (C)')
#    plt.savefig('ulostulolampotila.pdf')
    
    plt.figure(1)
    plt.plot(d*1000, dp/1000, 'k-d')
    plt.xlabel('d (mm)')
    plt.ylabel('Painehavio (kPa)')
#    plt.savefig('painehavio.pdf')

    
if __name__ == "__main__":
    main()