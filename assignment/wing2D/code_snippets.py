#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:05:27 2018

@author: ami
"""
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib import transforms


def read_ladson(plot=False):
    
    with open("LadsonForceData.txt") as ifile:
        lines = ifile.readlines()
    
    # 80 grit    
    n = 17
    alpha80 = sp.zeros(n)
    cl80 = sp.zeros(n)
    cd80 = sp.zeros(n)
    for k, line in enumerate(lines[5:5+n]):
        parts = line.split()
        alpha80[k] = float(parts[0])
        cl80[k] = float(parts[1])
        cd80[k] = float(parts[2])
        
    # 120 grit    
    n = 18
    alpha120 = sp.zeros(n)
    cl120 = sp.zeros(n)
    cd120 = sp.zeros(n)
    for k, line in enumerate(lines[23:23+n]):
        parts = line.split()
        alpha120[k] = float(parts[0])
        cl120[k] = float(parts[1])
        cd120[k] = float(parts[2])    
    
    # 180 grit    
    n = 18
    alpha180 = sp.zeros(n)
    cl180 = sp.zeros(n)
    cd180 = sp.zeros(n)
    for k, line in enumerate(lines[42:42+n]):
        parts = line.split()
        alpha180[k] = float(parts[0])
        cl180[k] = float(parts[1])
        cd180[k] = float(parts[2])       
        

    if plot:
        fig, axs = plt.subplots(2)
        ax = axs[0]
        ax.plot(alpha80, cl80, '-d')
        ax.plot(alpha120, cl120, '-d')
        ax.plot(alpha180, cl180, '-d')
        ax = axs[1]
        ax.plot(alpha80, cd80, '-d')
        ax.plot(alpha120, cd120, '-d')
        ax.plot(alpha180, cd180, '-d')
        
        axs[1].set_ylim(0.005, 0.03)
        axs[0].set_ylabel("$C_L$")
        axs[1].set_ylabel("$C_D$")
        axs[1].set_xlabel(r"$\alpha$ ($^\circ$)")
        
        axs[0].axhline(0, color='k')
        
        
def wing_shape():
    x = sp.linspace(0,1, 500)
    yp = (0.594689181*(0.298222773*sp.sqrt(x) - 0.127125232*x - 0.357907906*x**2 + 0.291984971*x**3 - 0.105174606*x**4) )
    yn= -(0.594689181*(0.298222773*sp.sqrt(x) - 0.127125232*x - 0.357907906*x**2 + 0.291984971*x**3 - 0.105174606*x**4) )
    
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(-15)
    
    plt.plot(x,yp,'k', transform= rot + base)
    plt.plot(x,yn,'k', transform= rot + base)
    plt.axis("off")
    plt.axis('equal')
    plt.savefig("wing_raw.png")


if __name__ == "__main__":
    print("START")
    wing_shape() 
    read_ladson(plot=True)
    print("END")

