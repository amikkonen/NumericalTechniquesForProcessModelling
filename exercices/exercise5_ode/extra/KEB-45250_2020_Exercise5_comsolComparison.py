# -*- coding: utf-8 -*-
import scipy as sp
from matplotlib import pyplot as plt
from scipy import optimize
from CoolProp.CoolProp import PropsSI

C2K = 273.15
# Graaviational acceleration
g = 9.81
# Stefan-Boltman constant
sigma = 5.67036713e-8

def read(path, first_index=0, last_index=1e99999):
    t, T_red, T_sensor = sp.loadtxt(path, delimiter=',',unpack=True)
    t *= 60 # to seconds 
    t = t[first_index:last_index] - t[first_index]
    T_red = T_red[first_index:last_index] + C2K
    T_sensor = T_sensor[first_index:last_index] + C2K

    return t, T_red, T_sensor
    
    return sp.array(t)-t[0], sp.array(T_red)+C2K, sp.array(T_sensor)+C2K

def read_comsol(path):
    with open(path, "r") as ifile:
        lines = ifile.readlines()
    t = []
    T = []
    for line in lines[8:]:
        print(line)
        parts = line.split()
        t.append(float(parts[0]))
        T.append(float(parts[1]))
    
    t = sp.array(t)*60 # to seconds 
    T = sp.array(T) + C2K
    return t, T




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
    
    # https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Natural%20Convection.pdf
    Ra = Gr*Pr
    if Ra < 1e4 or Ra > 1e13:
        raise ValueError("Ra < 1e4 or Ra > 1e13", "Ra=", Ra)
    elif Ra < 1e9:
        Nu = 0.59*Ra**0.25
    else:
        Nu = 0.1*Ra**(1/3)
    
    # Heat transfer coefficient
    h = Nu * kf / L
    
    return h


def analytic(T_i, T_inf, t, a):
    return T_inf + (T_i - T_inf)*sp.exp(-a*t)
    
def explicit(T_i, T_inf, dt, L, A, m_obj, cp_obj):
    #h = 5
    h = hVertical(T_i, T_inf, L)
    a = h*A/(m_obj*cp_obj)
    b = eps*sigma*A/(m_obj*cp_obj)
    T_new = T_i - a*(T_i-T_inf)*dt - b*(T_i**4-T_inf**4)*dt
    return T_new
    
def implicit(T_i, T_inf, dt, L, A, m_obj, cp_obj):

    # Compare to analytical
#    h = 5
#    a = h*A/(m_obj*cp_obj)
#    T_end = a*dt/(1+a*dt)*T_inf + T_i / (1+a*dt)

    def to_zero(T_end_guess):
        h = hVertical(T_end_guess, T_inf, L)
        a = h*A/(m_obj*cp_obj)
        b = eps*sigma*A/(m_obj*cp_obj)
        T_new = T_i - a*(T_end_guess-T_inf)*dt - b*(T_end_guess**4-T_inf**4)*dt 
        return T_end_guess - T_new
    res = optimize.root(to_zero, T_i)
    assert(res.success)

    return res.x[0]
    
def customised(T_i, T_inf, dt, L, A, m_obj, cp_obj, iteraitions=2):
    
    T_new = T_i
    for k in range(iteraitions):
        T_mean = (T_new+T_i) / 2
        h_rad = eps*sigma*(T_mean**4-T_inf**4) / (T_mean-T_inf)
        h_conv = hVertical(T_mean, T_inf, L)
        h = h_rad + h_conv
        a = h*A/(m_obj*cp_obj)
        T_new = analytic(T_i, T_inf, dt, a)

    return T_new



##########################################################################
# MAIN
##########################################################################

#--------------Measurement data------------------------------------
# Read, 
path = "temperature.log"
first_index = 1*60
last_index  = 41*60
T_inf = 23.5 + C2K

t_measurement, T_measured, T_sensor = read(path, first_index, last_index) 

#--------------Plate parameters------------------------------------
L_vertical   = 12.4e-2
L_horizontal = 2.9e-2
L_depth      = 5e-3
m_obj        = 48e-3
cp_obj       = 0.91e3 
k_obj        = 237
eps          = 0.96 
T_0   = T_measured[0]
#print(T_0)
 
#--------------Solution times------------------------------------
N = 7 # The customised method does not need this many time steps.

#--------------Help variables------------------------------------
A = 2*L_vertical*L_horizontal + 2*L_vertical*L_depth
t_end = t_measurement[-1]
t = sp.linspace(0, t_end, N)
dt = t_end/(N-1)

#--------------Analytic solution------------------------------------
h_analytic = 10
a_analytic = h_analytic*A/(m_obj*cp_obj) # hA/mc_p
T_analytic = analytic(T_0, T_inf, t, a_analytic)

#-------------Finite-difference-------------------------------------
T_explicit, T_implicit = sp.zeros_like(t)+T_0, sp.zeros_like(t)+T_0
T_customised = sp.zeros_like(t)+T_0
for i in range(len(t)-1):
    
    T_explicit[i+1] = explicit(T_explicit[i], T_inf, dt, 
                               L_vertical, A, m_obj, cp_obj)
    T_implicit[i+1] = implicit(T_implicit[i], T_inf, dt, 
                               L_vertical, A, m_obj, cp_obj)
    T_customised[i+1] = customised(T_customised[i], T_inf, dt, 
                               L_vertical, A, m_obj, cp_obj)

#--------------Plotting---------------------------------------------
#Results
plt.figure(0)
plt.plot(t_measurement/60, T_measured-C2K, 'k', label='Measured')
#plt.plot(t/60, T_analytic-C2K, 'b--', label='Analytic')
plt.plot(t/60, T_explicit-C2K, 'r*', label='Explicit')
plt.plot(t/60, T_implicit-C2K, 'g+', label='Implicit')
plt.plot(t/60, T_customised-C2K, 'yd', label='Customised')

t_comsol, T_comsol = read_comsol("comsolCorrelation.txt")
plt.plot(t_comsol/60, T_comsol-C2K, '-', label='Comsol')



plt.xlabel('Time (min)'), plt.ylabel('Temperature (C)')
plt.xlim(0, t_end/60)
plt.legend()

