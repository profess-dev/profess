# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
from scipy.optimize import curve_fit
import sys

def get_function_on_grid(function, shape, vectors, r0=None):
    """ (1) function must take 3d numpy arrays as input
        (2) minimum image convention brings each component to within (0.5,0.5)
        (3) expects box vectors in rows of 'vectors'
    """
    vectors = np.asarray(vectors)
    if r0 is None:
        r0 = np.zeros(3)
    else:
        r0 = np.linalg.inv(vectors.T).dot(r0) # convert to scaled coords
    # compute possible elements of dr
    x = np.arange(shape[0], dtype='float')/shape[0] - r0[0]
    y = np.arange(shape[1], dtype='float')/shape[1] - r0[1]
    z = np.arange(shape[2], dtype='float')/shape[2] - r0[2]
    # bring components within [0.5,0.5]
    x = x - np.rint(x)
    y = y - np.rint(y)
    z = z - np.rint(z)
    # create grids for vectorized calls
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    # convert to cartesian coordinates
    x = vectors[0,0]*xx + vectors[1,0]*yy + vectors[2,0]*zz
    y = vectors[0,1]*xx + vectors[1,1]*yy + vectors[2,1]*zz
    z = vectors[0,2]*xx + vectors[1,2]*yy + vectors[2,2]*zz
    # evaluate function
    return function(x,y,z)

def approximate_potential(get_energy, density, volume):
    d = 1e-6
    potential = np.empty(density.shape())
    energy_0 = get_energy(density)
    for i in range(density.shape()[0]):
        for j in range(density.shape()[1]):
            for k in range(density.shape()[2]):
                # nudge density up
                density[i,j,k] += d
                energy_p = get_energy(density)
                # approximate potential
                derivative = (energy_p-energy_0) / d
                potential[i,j,k] = derivative / (volume/density.size())
                # restore density
                density[i,j,k] -= d
    return potential

def approximate_stress(get_energy, box_vectors):
    d = 1e-4
    A = np.array(box_vectors).T
    stress = np.empty([3,3])
    for i in range(6):
        eps = np.zeros([3,3])
        # nudge up by tiny amount
        if i==0 or i==1 or i==2:
            eps[i,i] = d
        elif i==3:
            eps[1,2] = d
        elif i==4:
            eps[0,2] = d
        elif i==5:
            eps[0,1] = d
        A_eps = (np.eye(3) + eps).dot(A)
        ene_p = get_energy(A_eps.T)
        # nudge down by tiny amount
        A_eps = (np.eye(3) - eps).dot(A)
        ene_m = get_energy(A_eps.T)
        # compute approximate stress
        s = 1.0/np.linalg.det(A)*(ene_p-ene_m)/(2*d)
        if i==0 or i==1 or i==2:
            stress[i,i] = s
        elif i==3:
            stress[1,2] = stress[2,1] = s
        elif i==4:
            stress[0,2] = stress[2,0] = s
        elif i==5:
            stress[0,1] = stress[1,0] = s
    return stress

def approximate_forces(get_energy, xyz_coords, box_vectors):
    d = 1e-4
    a = np.asarray(box_vectors).T
    ai = np.linalg.inv(a)
    forces = np.empty(xyz_coords.shape)
    for i in range(xyz_coords.shape[0]): # loop over atoms
        xyz_i = np.copy(xyz_coords[i,:])
        for j in range(3):               # loop over x/y/z
            # add tiny amount to coordinate
            xyz_coords[i,j] += d
            x = ai.dot(xyz_coords[i,:]) # get scaled position
            x = x - np.floor(x)         # ensure position is within box
            xyz_coords[i,:] = a.dot(x)  # recover cartesian position
            ene_p = get_energy(xyz_coords)                
            xyz_coords[i,:] = xyz_i     # reset coordinates
            # subtract tiny amount from coordinate 
            xyz_coords[i,j] -= d
            x = ai.dot(xyz_coords[i,:]) # get scaled position
            x = x - np.floor(x)         # ensure position is within box
            xyz_coords[i,:] = a.dot(x)  # recover cartesian position
            ene_m = get_energy(xyz_coords)
            xyz_coords[i,:] = xyz_i     # reset coordinates
            # compute approximate force
            forces[i,j] = -(ene_p-ene_m)/(2*d)
    return forces

def murnaghan(vol, ene, plot=False):
    vol = np.asarray(vol)
    ene = np.asarray(ene)
    # initial guess is harmonic solid (E = E0 + 0.5*B0*(V-V0)^2/V0)
    apar, bpar, cpar= np.polyfit(vol, ene, 2)
    B0_g= -bpar
    V0_g= B0_g/(2*apar)
    E0_g= cpar-0.5*B0_g*V0_g
    B0prime_g= 3.5 
    # fit to murnaghan equation of state
    def murn(v, B0, B0prime, E0, V0):
        return  E0 + (B0*v/B0prime)*((((V0/v)**B0prime)/ \
                     (B0prime-1))+1)-B0*V0/(B0prime-1) 
    params, pcov = curve_fit(
        murn, vol, ene, p0=(B0_g,B0prime_g,E0_g,V0_g), maxfev=1000)
    # optional plotting
    if plot:
        import matplotlib
        matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
        plt.plot(vol,ene,'bo')
        vfit = np.linspace(0.99*vol[0],1.01*vol[-1]+0.25)
        efit = murn(vfit, params[0],params[1],params[2],params[3])
        plt.plot(vfit,efit,'b-')
        plt.legend(['data','fit'],loc='best')
        plt.show()
    return params
