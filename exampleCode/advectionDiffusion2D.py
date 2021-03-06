#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of steady 2D advection diffusion heat problem.


NOTE: USE OF SPARSE MATRIXES IN sparseMatrix.py

Created on Feb 12 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from matplotlib import pyplot as plt
import sparseMatrix
##############################################################################



def solver(nx, ny, kt, Lx, Ly, T_x0, T_xL, T_y0, T_yL, crho, ux, uy, 
           advection_scheme):
    
    n_tot = nx*ny
    
    # Cell sizes
    dx = Lx / nx
    dy = Ly / ny
    
    # Coefficients
    aCx  = kt*dy/dx
    aCy  = kt*dx/dy

    # Source vector
    b = sp.zeros(n_tot, dtype=sp.float_)
    # Coefficient matrix
    A = sparseMatrix.SparseMatrixA()


    # List of y boundary index
    yMin = sp.arange(nx)
    yMax = sp.arange(n_tot-nx, n_tot)

    #########################################################################
    # Build diffusion
    #########################################################################
    
    for k in range(n_tot):

        # Add aW to matrix A
        if k % nx != 0:
            A.add(k,k, aCx)
            A.add(k,k-1,-aCx)
        else:
            # NaN for zerogradient
            if not sp.isnan(T_x0):
                # Constant boundary
                A.add(k,k,2*aCx)
                b[k] += 2*aCx*T_x0    

        # Add aE to matrix A
        if (k + 1) % nx != 0:
            A.add(k,k, aCx)
            A.add(k,k+1,-aCx)
        else:
            # NaN for zerogradient
            if not sp.isnan(T_xL):
                # Constant boundary
                A.add(k,k,2*aCx)
                b[k] += 2*aCx*T_xL

        # Add aS to matrix A
        if k not in yMin:
            A.add(k,k, aCy)
            A.add(k,k-nx,-aCy)
        else:
            # NaN for zerogradient
            if not sp.isnan(T_y0):
                # Constant boundary
                A.add(k,k,2*aCy)
                b[k] += 2*aCy*T_y0

        # Add aN to matrix A
        if k not in yMax:
            A.add(k,k, aCy)
            A.add(k,k+nx,-aCy)
        else:
            # NaN for zerogradient
            if not sp.isnan(T_yL):
                # Constant boundary
                A.add(k,k,2*aCy)
                b[k] += 2*aCy*T_yL    
        
    #########################################################################
    # Build advection
    #########################################################################
    
    #########################################################################
    # Central differencing
    #########################################################################
    
    if advection_scheme == "central_differencing":
        aFx = crho*ux*dy/2
        aFy = crho*uy*dx/2
        
        for k in range(n_tot):        
            if ux != 0:
                # West
                if k % nx != 0:
                    A.add(k,k, -aFx)
                    A.add(k,k-1,-aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_x0):
                        # Constant boundary
                        b[k] += 2*aFx*T_x0
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
                
                # East
                if (k + 1) % nx != 0:
                    A.add(k,k, aFx)
                    A.add(k,k+1,aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_xL):
                        # Constant boundary
                        b[k] -= 2*aFx*T_xL
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
            if uy != 0:                
                # South
                if k not in yMin:
                    A.add(k,k, -aFy)
                    A.add(k,k-nx,-aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_y0):
                        # Constant boundary
                        b[k] += 2*aFy*T_y0
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
                # North
                if k not in yMax:
                    A.add(k,k, aFy)
                    A.add(k,k+nx,aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_yL):
                        # Constant boundary
                        b[k] -= 2*aFy*T_yL
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
                    
    #########################################################################
    # Upwind 
    #########################################################################                
    
    elif advection_scheme == "upwind":
        aFx = crho*ux*dy
        aFy = crho*uy*dx
        
        for k in range(n_tot):  
            # Positive velocity
            if ux > 0:
                # West
                if k % nx != 0:
                    A.add(k,k-1,-aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_x0):
                        # Constant boundary
                        b[k] += aFx*T_x0
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
                # East
                if (k + 1) % nx != 0:
                    A.add(k,k, aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_xL):
                        # Constant boundary
                        A.add(k,k, aFx)
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")

            if ux < 0:            
                # West
                if k % nx != 0:
                    A.add(k,k,-aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_x0):
                        # Constant boundary
                        A.add(k,k,-aFx)
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
                # East
                if (k + 1) % nx != 0:
                    A.add(k,k+1, aFx)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_xL):
                        # Constant boundary
                        b[k]   -= aFx*T_xL
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
            
            if uy > 0:
                # South
                if k not in yMin:
                    A.add(k,k-nx,-aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_y0):
                        # Constant boundary
                        b[k] += aFy*T_y0
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
        
                # North
                if k not in yMax:
                    A.add(k,k, aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_yL):
                        # Constant boundary
                        A.add(k,k, aFy)
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
            if uy < 0: 
                # South
                if k not in yMin:
                    A.add(k,k,-aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_y0):
                        # Constant boundary
                        A.add(k,k,-aFy)
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")
        
                # North
                if k not in yMax:
                    A.add(k,k+nx, aFy)
                else:
                    # NaN for zerogradient
                    if not sp.isnan(T_yL):
                        # Constant boundary
                        b[k]   -= aFy*T_yL
                    else:
                        raise NotImplementedError("NOT IMPLEMENTED!")




    # Solver and post
    T = A.solve(b)
    # [iy,ix]
    Tar = T.reshape((ny,nx))
    return Tar
    



##############################################################################
def analytical1D(n, kt, L, T_A, T_B, x, crho, u):
    if u != 0:
        return T_A + (T_B-T_A) * (sp.exp(crho*u*x/kt)-1)/(sp.exp(crho*u*L/kt)-1)
    else:
        return T_A + (T_B-T_A) * x /L

def test1Dx():
    # Thermal conductivity (W/mK)
    kt = 0.1
    # Lenght (m)
    Lx = 1
    Ly = 1
    # Boundary temperatures (C)
    T_x0 = 1
    T_xL = 0.5
    
    T_y0 = sp.nan
    T_yL = sp.nan

    advection_scheme = "central_differencing"
#    advection_scheme = "upwind"
    
    # Product of heat capacity and density
    crho = 1
    
    # Velocity and Number of control volumes
#    # Case a
    nx = 5;  ux = 0.1
    ny = 5;  uy = 0

    Tar = solver(nx, ny, kt, Lx, Ly, T_x0, T_xL, T_y0, T_yL, 
                   crho, ux, uy, advection_scheme)
    
    dx = Lx / nx
    xCFD = sp.linspace(dx/2,Lx-dx/2,nx)  
    xana = sp.linspace(0,Lx)
    Tana1D = analytical1D(nx, kt, Lx, T_x0, T_xL, xana, crho, ux)
    
    plt.plot(xana, Tana1D)
    plt.plot(xCFD, Tar[0],'-dk')
    print(Tar[0])
#    print("[ 0.93373341  0.7879469   0.6130031   0.40307053  0.15115145]")
    print("[ 0.84884855  0.59692947  0.3869969   0.2120531   0.06626659]")


def test1Dy():
    # Thermal conductivity (W/mK)
    kt = 0.1
    # Lenght (m)
    Lx = 1
    Ly = 1
    # Boundary temperatures (C)
    T_x0 = sp.nan
    T_xL = sp.nan
    
    T_y0 = 1
    T_yL = 0.5

    advection_scheme = "central_differencing"
#    advection_scheme = "upwind"
    
    # Product of heat capacity and density
    crho = 1
    
    # Velocity and Number of control volumes
#    # Case a
    nx = 5; ux = 0
    ny = 5;  uy = 0.1

    Tar = solver(nx, ny, kt, Lx, Ly, T_x0, T_xL, T_y0, T_yL, 
                   crho, ux, uy, advection_scheme)
    
    dy = Ly / ny
    yCFD = sp.linspace(dy/2,Ly-dy/2,ny)  
    yana = sp.linspace(0,Ly)
    Tana1D = analytical1D(ny, kt, Ly, T_y0, T_yL, yana, crho, uy)
    
    plt.plot(yana, Tana1D)
    plt.plot(yCFD, Tar[:,0],'-dk')
    print(Tar[:,0])
#    print("[ 0.93373341  0.7879469   0.6130031   0.40307053  0.15115145]")
    print("[ 0.84884855  0.59692947  0.3869969   0.2120531   0.06626659]")


def main():
    
    # Thermal conductivity (W/mK)
    kt = 0.
    # Lenght (m)
    Lx = 1
    Ly = 1
    # Boundary temperatures (C)
    T_x0 = 1
    T_xL = 0
    
    T_y0 = T_xL
    T_yL = T_x0

#    advection_scheme = "central_differencing"
    advection_scheme = "upwind"
    
    # Product of heat capacity and density
    crho = 1
    
    # Velocity and Number of control volumes
#    # Case a
    nx = 11;  ux = 1
    ny = nx;  uy = ux

    Tar = solver(nx, ny, kt, Lx, Ly, T_x0, T_xL, T_y0, T_yL, 
                   crho, ux, uy, advection_scheme)
    
    plt.matshow(Tar,aspect='equal',
#                cmap='gray'
                cmap='jet'
                )
    ax = plt.gca()
    ax.set_xticks(sp.arange(0, nx, 1));
    ax.set_yticks(sp.arange(0, ny, 1));
    plt.colorbar(label="T (C)")
    plt.xlabel("x")
    plt.ylabel("y")
    ax.invert_yaxis()
    for k in range(nx):
        plt.axvline(k+0.5, color='k')
    for k in range(ny):
        plt.axhline(k+0.5, color='k')        
    plt.axhline(int(ny/2), color="r")
    plt.axvline(int(nx/2), color="b")
    plt.plot([0,10], [10,0], "--k")
    
    plt.annotate('', xy=(1,1), xytext=(-1, -1),
            arrowprops=dict(facecolor='black', shrink=0.01),
            )
    plt.text(-1,-2,"u", fontsize=15)
    
    plt.savefig("numericalDiffusion2D.pdf", figsize=(3.5,2.5))

    fig, axs = plt.subplots(3)
    dx = Lx / nx
    dy = Ly / ny
    x = sp.linspace(dx/2,Lx-dx/2,nx)
    y = sp.linspace(dy/2,Ly-dy/2,ny)
    ax = axs[0]
    ax.plot(x, Tar[int(ny/2)], 'r')
    ax.set_xlim(x.min(), x.max())
    ax = axs[1]
    ax.plot(y, Tar[:,int(nx/2)], "b")
    ax.set_xlim(y.min(), y.max())
    ax = axs[2]
    
    TnumDiff = []
    for k in range(nx):
        TnumDiff.append(Tar[ny-k-1, k])    
    ax.plot(sp.arange(nx), TnumDiff, "--k", label="numerical")
    ax.plot([0,int(nx/2),int(nx/2),nx], [T_x0,T_x0,T_xL,T_xL], "-k", label="exact")
    ax.set_xlim(0, nx)
    ax.legend(frameon=False)
    for ax in axs:
        ax.set_ylim(T_xL-0.1, T_x0+0.1)
        ax.set_ylabel("T (C)")
    
    fig.tight_layout()
    fig.savefig("numericalDiffusion2Dcurves.pdf",figsize=(3.5,2.5))
    

def false_diffusion():
    
    # Thermal conductivity (W/mK)
    kt = 0.
    # Lenght (m)
    Lx = 1
    Ly = 1
    # Boundary temperatures (C)
    T_x0 = 1
    T_xL = 0
    
    T_y0 = T_xL
    T_yL = T_x0

#    advection_scheme = "central_differencing"
    advection_scheme = "upwind"
    
    # Product of heat capacity and density
    crho = 1
    
    # Velocity and Number of control volumes
#    # Case a
    nx = 11;  ux = 1
    ny = nx;  uy = ux

    Tar = solver(nx, ny, kt, Lx, Ly, T_x0, T_xL, T_y0, T_yL, 
                   crho, ux, uy, advection_scheme)
    
    plt.matshow(Tar,aspect='equal',
#                cmap='gray'
                cmap='jet'
                )
    ax = plt.gca()
    ax.set_xticks(sp.arange(0, nx, 1));
    ax.set_yticks(sp.arange(0, ny, 1));
    plt.colorbar(label="T (C)",anchor=(1.0,0.0))
#    plt.xlabel("x")
#    plt.ylabel("y")
    plt.axis('off')


    ax.invert_yaxis()
    for k in range(nx+1):
        plt.axvline(k-0.5, color='k')
    for k in range(ny+1):
        plt.axhline(k-0.5, color='k')        
    plt.axhline(int(ny/2), color="r")
    plt.axvline(int(nx/2), color="b")
    plt.plot([0,10], [10,0], "--k")
    
    plt.annotate('', xy=(1,1), xytext=(-1, -1),
            arrowprops=dict(facecolor='black', shrink=0.01),
            )
    plt.text(-1,-2,"u", fontsize=15)
    plt.text(-1.5,5.5,"T=1 C", fontsize=15, rotation=90)
    plt.text(11,5.5,"T=0 C", fontsize=15, rotation=90)
    plt.text(4,-1.5,"T=0 C", fontsize=15, )
    plt.text(4,11,"T=1 C", fontsize=15, )
    
    plt.savefig("numericalDiffusion2D.pdf", figsize=(3.5,2.5))

    fig, axs = plt.subplots(3)
    dx = Lx / nx
    dy = Ly / ny
    x = sp.linspace(dx/2,Lx-dx/2,nx)
    y = sp.linspace(dy/2,Ly-dy/2,ny)
    ax = axs[0]
    ax.plot(x, Tar[int(ny/2)], 'r')
    ax.set_xlim(x.min(), x.max())
    ax = axs[1]
    ax.plot(y, Tar[:,int(nx/2)], "b")
    ax.set_xlim(y.min(), y.max())
    ax = axs[2]
    
    TnumDiff = []
    for k in range(nx):
        TnumDiff.append(Tar[ny-k-1, k])    
    ax.plot(sp.arange(nx), TnumDiff, "--k", label="numerical")
    ax.plot([0,int(nx/2),int(nx/2),nx], [T_x0,T_x0,T_xL,T_xL], "-k", label="exact")
    ax.set_xlim(0, nx)
    ax.legend(frameon=False)
    for ax in axs:
        ax.set_ylim(T_xL-0.1, T_x0+0.1)
        ax.set_ylabel("T (C)")
    
    fig.tight_layout()
    fig.savefig("numericalDiffusion2Dcurves.pdf",figsize=(3.5,2.5))
if __name__ == "__main__":
    print("START")
#    test1Dx()
    test1Dy()
#    main()
#    false_diffusion()
    print("END")