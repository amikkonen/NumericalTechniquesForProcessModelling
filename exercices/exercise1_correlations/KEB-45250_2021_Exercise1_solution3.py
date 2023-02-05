# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import json
from CoolProp.CoolProp import PropsSI

def read_data():
    """This is a comment descriving this function.
    
    Reads data from stagnationPointNu.json and returns it.
    """
    with open("stagnationPointNu.json", "r") as ifile:
        data = json.load(ifile)
    return data

def get_data(d, data):
    """Get the cases for which d == d in data and returns Re and Nu0.
    
    """
    Re = []
    Nu0 = []
    for item in data:
        if item["d"] == d:
            Re.append(item["Re"])
            Nu0.append(item["Nu0"])
    return np.array(Re), np.array(Nu0)
            

def Nu0_correlation_for_large_d(Re, Pr, Hd):
    Nu0 = 0.56*Re**0.53*Pr**0.4*Hd**0.064
    return Nu0

def Nu0_correlation_for_small_d(Re, Pr, Hd):
    Nu0 = 0.5*Re**0.3*Pr**0.4*Hd**0.064
    return Nu0


def make_correlation():
    # Prandtl 
    Pr = 0.7
    # Dimendionless distance from the wall
    Hd = 2
    # Get data from file
    data = read_data()
    # Re range for plotting correlations
    Res = np.linspace(0,100000)
    
    ########################################################
    # Large diameter nozzles
    ########################################################
    
    # d=8mm
    Re, Nu0 = get_data(8e-3, data)
    plt.plot(Re,Nu0, "d", label="d=8mm")
    # d=6mm
    Re, Nu0 = get_data(6e-3, data)
    plt.plot(Re,Nu0, "d", label="d=6mm")
    
    # Correlation for large
    Nu0_c = Nu0_correlation_for_large_d(Res, Pr, Hd)
    plt.plot(Res,Nu0_c, "-", label=r"$\mathrm{Nu}_{0}=0.56\mathrm{Re}^{0.53}\mathrm{Pr^{0.4}}\left(H/d\right)^{0.064}$")
    
    
    ########################################################
    # Small diameter nozzle
    ########################################################
    
    # d=1.1mm
    Re, Nu0 = get_data(1.1e-3, data)
    plt.plot(Re,Nu0, "d", label="d=1.1mm")
    
    # Correlation for small
    Nu0_c = Nu0_correlation_for_small_d(Res, Pr, Hd)
    plt.plot(Res,Nu0_c, "-", label="correlation for small pipe")
    
    
    
    ########################################################
    # Make a prettyer plot
    ########################################################
    plt.legend()
    plt.title("Stagnation point heat transfer")
    plt.xlabel("Re")
    plt.ylabel("Nu0")
    
def solve_heat_flux():
    d = 6e-3
    Hd = 2
    V = 50
    T_w = 50+273
    T_n = 20+273
    mu  = PropsSI("V", "T", T_n, "P", 1e5, "air")
    rho = PropsSI("D", "T", T_n, "P", 1e5, "air")
    nu  = mu/rho
    Pr = PropsSI("Prandtl", "T", T_n, "P", 1e5, "air")
    kf = PropsSI("conductivity", "T", T_n, "P", 1e5, "air")
    
    Re = V*d/nu
    Nu0 = Nu0_correlation_for_large_d(Re, Pr, Hd)
    
    print("Re", Re)
    print("Nu0", Nu0)
    
    h = Nu0*kf/d
    print("h", h)
    dT = T_w-T_n
    q=h*dT
    print("q",q,"W/m2K")

    
if __name__ == '__main__':
    make_correlation()
    solve_heat_flux()