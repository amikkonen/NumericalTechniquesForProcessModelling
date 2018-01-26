#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on

@author: a.mikkonen@iki.fi
"""
import time
import scipy as sp
from matplotlib import pyplot as plt

def power_law_profile(nn):
    n = 7
    yperR = sp.linspace(0,1, 1000)
    u_rat = (yperR)**(1/n)    
#    plt.plot(yperR, u_rat, 'k:', label=r"$\frac{\overline{u}}{u_{max}}=\frac{y}{R}^{1/n}$")
    plt.plot(yperR, u_rat, 'k:', label="pipe velocity profile")
    
    
    fig = plt.gcf()
    fig.set_size_inches(3.5,2.5)
    
    yperR_disc = sp.linspace(0,1, nn+1)
    yperR_f = (yperR_disc[:-1] + yperR_disc[1:])/2
    u_rat_disc = (yperR_f)**(1/n)    
    
    for k in range(len(yperR_disc)-1):
        plt.plot(yperR_disc[k:k+2], [u_rat_disc[k],u_rat_disc[k]], 'k')
        
    plt.plot(yperR_disc[k:k+1], [u_rat_disc[k]], 'k', label="discretization")    
    
    ax = plt.gca()
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlim(-0.01, 1)
    plt.ylim(0, 1.1)
    
    plt.xlabel("$y/R$")
    plt.ylabel("$u/u_{max}$")
    
    plt.tight_layout() 
    plt.legend(frameon=False, loc='lower right')    
    plt.savefig("wallProfile"+str(nn)+".pdf")

def main():
    pass
    
if __name__ == "__main__":
    start = time.time()
    print("START")
#    main()
    power_law_profile(10)
#    power_law_profile(50)
    print("END %.4f s" % (time.time()-start))


