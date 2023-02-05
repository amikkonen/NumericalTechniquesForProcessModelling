# -*- coding: utf-8 -*-
"""
Very simple real time plotter for Grove - Digital Infrared Temperature Sensor
    http://wiki.seeedstudio.com/Grove-Digital_Infrared_Temperature_Sensor/
linked to Arduino. 

Read data from serial and plots temperature real time.

Easy to modify for other similar data plotting from serial.

Antti Mikkonen, a.mikkonen@iki.fi, winter 2019.

"""

import matplotlib.pyplot as plt
import scipy as sp
from CoolProp.CoolProp import PropsSI

C2K = 273.15
# Graaviational acceleration
g = 9.81
# Stefan-Boltman constant
sigma = 5.67036713e-8
    
def read(path, first_index=0, last_index=1e99999):
    with open(path, "r") as ifile:
        lines = ifile.readlines()[first_index:last_index]
        
    T_red = []
    T_sensor = []
    t = []
    
    for line in lines:
        parts = line.split(",")
        t.append(float(parts[0])*60)
        T_red.append(float(parts[1]))
        T_sensor.append(float(parts[2]))
    
    return sp.array(t)-t[0], sp.array(T_red)+C2K, sp.array(T_sensor)+C2K

def hVertical(T, T_inf, L):
    """
    Natural convection for vertical plate.
    
    Correlation grapped from:
        https://en.wikipedia.org/wiki/Natural_convection#Natural_convection_from_a_vertical_plate
        https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Natural%20Convection.pdf
    Applicapility or primary sources have not been checked! 
    Do not use in real applications without chekking! 
    
    For convinience, we use CoolProp here for fluid properties. In real 
    applications it would be adviceble to find same faster solution. CoolProps 
    is slow and deep inside iterative loops performace matters.
    """
    # Reference pressure
    p = 1e5
    fluid = "air"
    T_ref = (T+T_inf)/2
    
    # density
    rho = PropsSI("D", "T", T_ref, "P", p, fluid)
    # dynamic viscosity, directly available in CoolProps
    mu  = PropsSI("V", "T", T_ref, "P", p, fluid)
    # Kinematic viscosity
    nu  = mu/rho
    # Prandtl number
    Pr  = PropsSI("Prandtl", "T", T_ref, "P", p, fluid)
    # Heat conductivity
    kf  = PropsSI("conductivity", "T", T_ref, "P", p, fluid)
    # Coefficient of thermal expansion (equal to approximately 1/T, for ideal gases)
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_ref, "P", p, fluid)

    # Grashof number, a typical dimensionless number for natural convection
    Gr = (beta*g*L**3*(T-T_inf)) / (nu**2)
    
    # Nusslet number
#    Nu = 0.478*Gr**0.25
    
    # https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Natural%20Convection.pdf
    Ra = Gr*Pr
    if Ra < 1e4 or Ra > 1e13:
        raise ValueError("Ra < 1e4 or Ra > 1e13", "Ra=", Ra)
    elif Ra < 1e9:
#        print("small")
        Nu = 0.59*Ra**0.25
    else:
#        print("big")
        Nu = 0.1*Ra**(1/3)
    
    # Heat transfer coefficient
    h = Nu * kf / L
    
    return h

def calculate(t, T_init, T_inf, L_vertical, L_horizontal, L_depth,
              m, cp, eps, k_obj):

    # No up or down
    A = (
             2*L_vertical*L_horizontal 
           + 2*L_vertical*L_depth
        )

    Ts =  sp.zeros(len(t))
    Ts[0] = T_init

    hs  =  sp.zeros_like(t)
    h_cs =  sp.zeros_like(t)
    h_rs =  sp.zeros_like(t)
    
    # Explicit time stepping
    for k in range(len(t)-1):
        
        T   = Ts[k]
        dt  = t[k+1] - t[k]
        h_c   = hVertical(T, T_inf, L_vertical)
        Q_c = h_c*A*(T - T_inf)*dt
        Q_r = A*eps*sigma*(T**4-T_inf**4)*dt
        Q   = Q_c + Q_r
        dT  = Q/(m*cp) 
        Ts[k+1] = T-dT
        
        # Convections
        h_r = Q_r / (A*(T - T_inf)*dt)
        h = h_c + h_r
        
        # Add to storage
        hs[k]   = h
        h_cs[k] = h_c
        h_rs[k] = h_r
        
        # Biot number, if Bi << 1 -> isothermal
        Bi_depth      = L_depth*h/k_obj
        Bi_vertical   = L_vertical*h/k_obj
        Bi_horizontal = L_horizontal*h/k_obj
        assert (Bi_depth<1e-3)
        assert (Bi_vertical<1e-2)
        assert (Bi_horizontal<1e-2)        
    
    # Touch
    hs[-1] = hs[-2]
    h_rs[-1] = h_rs[-2]
    h_cs[-1] = h_cs[-2]
    
    return Ts, hs, h_cs, h_rs


    
def plot(t, T_red, T_sensor, T_calc, 
         T_calc_low_eps, T_calc_high_eps,
         T_inf, 
         h, h_c, h_r, 
         eps, eps_low, eps_high):
    fig, axes = plt.subplots(2, figsize=(4,8))
    tmin = t/60
    
    # Temperature
    ax = axes[0]
    ax.plot(tmin, T_red, "r", label="IR")
#    ax.plot(t_RKmin, T_RK, "y--",  label="RK45")
    ax.fill_between(tmin, T_calc_low_eps, T_calc_high_eps, color="g",label="$\mathrm{Nu_c} = 0.59\mathrm{Ra}^{0.25}$\n$q_r=\epsilon \sigma (T^4-T_\infty^4)$"+ ", $\epsilon="+str(eps_low)+"-"+str(eps_high)+"$")
    ax.plot(tmin, T_calc, "k--",  label="$\mathrm{Nu_c} = 0.59\mathrm{Ra}^{0.25}$\n$q_r=\epsilon \sigma (T^4-T_\infty^4)$, $\epsilon="+str(eps)+"$")
    ax.plot(tmin, T_sensor, "b",  label="sensor")
    ax.axhline(T_inf, color="k", label="enviroment")
    ax.set_ylabel("T (C)")
    
    # htc
    ax = axes[1]
    ax.plot(tmin, h, "r",  label="total")
    ax.plot(tmin, h_c, "b",  label="convection")
    ax.plot(tmin, h_r, "y",  label="radiation")
    ax.set_ylabel("h (W/m2K)")
    
    # Prettyfy
    for ax in axes.ravel():
        ax.legend(frameon=False)    
        ax.set_xlim(tmin.min(), tmin.max())
        ax.set_xlabel("t (min)")
        
    fig.savefig("verticalPlate.pdf")
    
def main():
    # Log parameters
    path = "temperature.log"
    first_index = 1*60
    last_index  = 41*60
    
    # Plate parameters
    L_vertical = 12.4e-2
    L_horizontal = 2.9e-2
    L_depth = 5e-3
    m_obj = 48e-3
    cp_obj = 0.91e3 #J/kgK # https://www.engineeringtoolbox.com/specific-heat-metals-d_152.html
    k_obj = 237
    
    # emissivity
    # https://www.flir.com/discover/rd-science/
    # use-low-cost-materials-to-increase-target-emissivity/
    eps = 0.96 
    eps_low = 0.8
    eps_high = 1
    
    #################################################################
    # Read
    t, T_red, T_sensor = read(path, first_index, last_index)    

    # Calculate
    T_init = T_red[0]
    T_inf  = 23.5+C2K
    
    # Main
    T_calc_low_eps,h, h_c, h_r  = calculate(t, T_init, T_inf, 
                       L_vertical, L_horizontal, L_depth,
                       m_obj, cp_obj, eps_low, k_obj)

    # High eps
    T_calc_high_eps ,h, h_c, h_r  = calculate(t, T_init, T_inf, 
                       L_vertical, L_horizontal, L_depth,
                       m_obj, cp_obj, eps_high, k_obj)

    # Low eps
    T_calc,h, h_c, h_r  = calculate(t, T_init, T_inf, 
                       L_vertical, L_horizontal, L_depth,
                       m_obj, cp_obj, eps, k_obj)
    
    # Plot
    plot(t, T_red-C2K, T_sensor-C2K, T_calc-C2K, 
         T_calc_low_eps-C2K, T_calc_high_eps-C2K,
         T_inf-C2K,
         h, h_c, h_r, eps, eps_low, eps_high)
    
if __name__ == "__main__":
    main()
