#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on

@author: a.mikkonen@iki.fi
"""
import time
import scipy as sp
from matplotlib import pyplot as plt

def main():
    n = 100
    x = sp.linspace(0,1,n)
    y = sp.random.randn(n)/10 + 1
    
    plt.figure(figsize=(3.5,2.5))
    plt.plot(x,y,'k')
    m = y.mean()
    plt.axhline(m, color='k')
    plt.xlim(0,1)
    plt.ylim(0,2)
    
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.text(1.05, m, "$\overline{u}$" )
    ax.text(0.3, m+0.2, "$u$" )
    
    
    plt.xlabel("t")
    plt.ylabel("u")
    plt.xticks([])
    plt.yticks([])
    
    plt.savefig("turbulence.pdf")
    
if __name__ == "__main__":
    start = time.time()
    print("START")
    main()
    print("END %.4f s" % (time.time()-start))



