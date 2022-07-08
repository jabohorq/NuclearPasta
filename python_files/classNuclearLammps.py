#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 10:56:38 2022
NuclearLammps class:
    Set of funcions and solutionsto simulate nuclear pasta using 
    the results of Molecular Dynamics simulations obtained using LAMMPS.
@author: jbohorquez
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from itertools import product as iterate
#from quantimpy import minkowski as mk

class NuclearLammps(object):
    
    def __init__(self, lammps_file_name):
        '''
        Constructor of the NuclearLammps class.
        Take as argument the name of the dump file with the coordinates 
        of the nucleons. In the dump file we can have more that one set
        coordintate (configurations for different times).
        After the dump file is read, the properties:
        lammps_file_name (name  of the dump file), 
        N_par (number of nucleons), 
        and L_sim_box (lenght of the simulation box)
        are created
        '''
        self.lammps_file_name = lammps_file_name
        d_file = open(self.lammps_file_name, 'r')
        # Timestep
        d_file.readline()
        d_file.readline()
        # read number of particles
        d_file.readline()
        N_par = int(d_file.readline())
        self.N_par = N_par
        d_file.readline()
        L0, L_sim_box = d_file.readline().split()
        L_sim_box = float(L_sim_box)
        self.L_sim_box = L_sim_box
        num_lines = sum(1 for line in open(self.lammps_file_name))
        num_conf = num_lines /(9 + N_par)
        self.num_conf = int(num_conf)
##################################################################
    def print_info(self):
        """
        print_info: method to print class variable or properties
        Returns
        -------
        None.
        """
        print('***********************************')
        print('Class variables')
        print('***********************************')
        print('Dump file name: {}'.format(self.lammps_file_name))
        print('Number of particles: {0:d}'.format(self.N_par))
        print('Simulation box length: {0:1.2e} A'.format(self.L_sim_box))
        print('Number of configurations in this file: {0:d}'.format(self.num_conf))
        print('***********************************')
##################################################################
    def read_conf(self, n_conf=0):
        """
        read_conf method: reads lammps dump file and extracts the n_conf 
        configuration.

        Parameters
        ----------
        n_conf : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        Yp: proton fraction
        x_proton, y_proton, z_proton: numpy arrays with proton coord
        x_neutron, y_neutron, z_neutron: numpy arrays with neutron coord
        """
        # n_conf: number of configuration to read
        d_file = open(self.lammps_file_name, 'r')
        for i in range(9+(self.N_par + 9)*int(n_conf)):
            d_file.readline()
        x_data = np.zeros(self.N_par)
        y_data = np.zeros(self.N_par)
        z_data = np.zeros(self.N_par)
        type_data = np.zeros(self.N_par)
        #read coordinates and type
        for i in range(self.N_par):
            idp, type_p, xx, yy, zz = d_file.readline().split()
            xx = float(xx)
            yy = float(yy)
            zz = float(zz)
            type_p = int(type_p)
            x_data[i] = xx
            y_data[i] = yy
            z_data[i] = zz
            type_data[i] = type_p
        d_file.close()
        ind_neutron = [i for i, n in enumerate(type_data) if n == 1]
        ind_proton = [i for i, n in enumerate(type_data) if n == 2]
        # proton positions
        x_proton = np.array([x_data[i] for i in ind_proton])
        y_proton = np.array([y_data[i] for i in ind_proton])
        z_proton = np.array([z_data[i] for i in ind_proton])
        # neutron positions
        x_neutron = np.array([x_data[i] for i in ind_neutron])
        y_neutron = np.array([y_data[i] for i in ind_neutron])
        z_neutron = np.array([z_data[i] for i in ind_neutron])
        # number of protons
        n_proton = len(ind_proton)
        # Proton fraction
        Yp = n_proton / self.N_par
        return Yp, x_proton, y_proton, z_proton,\
    x_neutron, y_neutron, z_neutron
##################################################################
## END OF CLASS
##################################################################
def plot_conf(x,y,z,outname):
    """
    plot_conf: quick 3d plot of a set op particles with 
    coordinates stored in numpy arrays x, y, z.
    Saves plot in file with name outname 

    Parameters
    ----------
    x : numpy array
        x coordinates
    y : numpy array
        y coordinates.
    z : numpy array
        z coordinates.
    outname : string
        plot out file name.

    Returns
    -------
    None.
    """
    # 3d plots
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z, c='r', marker='o')
    plt.savefig(outname, dpi=300)
    plt.show()
##################################################################
# Integer coordinates transformation
def integer_coord(x, div):
    """
    integer_coord: when transforming coordinates to bins, this function
    takes care of the periodic boundary conditions

    Parameters
    ----------
    x : float
        coordinate.
    div : float
        number of divisions or bins.

    Returns
    -------
    xn : integer
        transformed coordinate.

    """
    if x == 0:
        xn = int(0)
    else:
        xn = int(x/div+0.5)- 1
    return xn
##################################################################
# wrapped integer coordinates
def wrapped_coord(x, n, bins):
    """
    wrapped_coord: periodic boundary conditions for bins coordinates
    when finding closets neighbors

    Parameters
    ----------
    x : float
        DESCRIPTION: bin coordinate
    n : integer
        DESCRIPTION: neighbor cell
    bins : integer
        DESCRIPTION. total number of bins

    Returns
    -------
    xn : float
        DESCRIPTION. coordinate of n neighbor

    """
    if (x+n) < 0:
        xn = x + bins + n
    else:
        xn = x + n
    return xn
##################################################################
# distance for periodic boundary conditions
def dist_pbc(r, bins):
    """
    dist_pbc: calculate distance when using PBC (to be use in 
    potential calculation)

    Parameters
    ----------
    r : real
        DESCRIPTION. distance
    bins : integer
        DESCRIPTION. Total number of bins
        
    Returns
    -------
    new_r : ral
        DESCRIPTION. transformed distace

    """
    new_r = r - bins * int(r/bins)
    return new_r
##################################################################
def voxel(xp,yp,zp,L_sim_box, bins, n_th, name_seq):
    """
    voxel: takes discrete distribution of protons and create a discrete map 
    of charge assuming the charge has a gaussian distribution.
    Also creates a gray-scale map of the charge distribution to calculate
    the Minkowski functionals.
    Finally, creates text files with slices of the gray-scale map

    Parameters
    ----------
    xp, yp, zp : float arrays
        DESCRIPTION. coordinate array of proton coordinates
    L_sim_box : float
        DESCRIPTION. simulation box lenght 
    bins : integer
        DESCRIPTION. total number of bins
    n_th : real
        DESCRIPTION. threshold charge distribution
    name_seq : string
        DESCRIPTION. name for the slice files 

    Returns
    -------
    div : real
        DESCRIPTION. simulation box length / total number of bins (same as
    resolution lenght/bin)
    greytones : array
        DESCRIPTION. 3d array with the gray-scale map
    voxels : array
        DESCRIPTION. 3d array with the binary map

    """
    div = L_sim_box / bins
    print('div: {0}'.format(div))
    delta_v = div*div*div
    print('delta_v: {0}'.format(delta_v))
    sigma = 1.5e-5
    charges = np.zeros([bins,bins,bins]) #x,y,z
    n_protons = len(xp) # number of protons
    print('n_protons: {0}'.format(n_protons))
    for m in range(n_protons):
        xt = integer_coord(xp[m], div)
        yt = integer_coord(yp[m], div)
        zt = integer_coord(zp[m], div)
    
        for i,j,k in iterate(range(-2,3),repeat=3):
            xtt = wrapped_coord(xt, i, bins)
            ytt =  wrapped_coord(yt, j, bins)
            ztt =  wrapped_coord(zt, k, bins)

            delta_x = xp[m]/div - xtt -1
            r_x = dist_pbc(delta_x, bins)
            delta_y = yp[m]/div - ytt -1
            r_y = dist_pbc(delta_y, bins)
            delta_z = zp[m]/div - ztt -1
            r_z = dist_pbc(delta_z, bins)
            r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
            # 3D normal dist
            sigma_n = sigma/div
            konst = 2*np.pi*sigma_n**2
            
            norm3d = norm.pdf(r,loc=0,scale=sigma_n)/konst
            
            charges[xtt%bins, ytt%bins, ztt%bins] += norm3d
            
    voxels = np.zeros([bins, bins, bins])

    # Threshold density = 0.03 fm-3 = 0.03e15 A-3
    # threshold in bin units: 0.03e15 * div^3
    #n_th = 0.03e15
    charge_th = n_th * delta_v
    print('threshold charge: {}'.format(charge_th))
    for i,j,k in iterate(range(bins),repeat=3):
         if charges[i,j,k] > charge_th:
             voxels[i,j,k] = 1
            
    voxels = np.array(voxels, dtype=bool)
    greytones = np.zeros([bins, bins, bins]) #one greytone at each voxel
    for i,j,k in iterate(range(bins),repeat=3):
    #at each grey tone, go through the 7 values within the adjacent cell,
    #and count their voxel values with the following equation:
        for l,m,n in iterate([0,1],repeat=3):
            greytones[i,j,k]+=2**(l+2*m+4*n)*voxels[(i+l)%bins,\
                                (j+m)%bins,(k+n)%bins]
    seq_slice(name_seq, greytones, bins)
    # testvoxels is merely a list of the voxel coordinates with values=1
    # testvoxels = np.array([[i,j,k] \
    #         for i,j,k in iterate(range(bins),repeat=3) if voxels[i,j,k]])
    return div, greytones, voxels
##################################################################
def seq_slice(name_seq, arr3d, bins):
    '''
    Save slices of charge distribution into text files
    Parameters
    ----------
    name_seq : TYPE String
        DESCRIPTION: root name for the sequence of slices
    arr3d : TYPE 3d array
        DESCRIPTION. 3d array with active voxels
    bins : TYPE integer
        DESCRIPTION. Number of bins (per unit length)

    Returns
    -------
    None.
    '''
    for i in range(bins):
        file_seq_name = name_seq + '%03d.txt'% i
        print(file_seq_name)
        #print(arr3d[:,:,i])
        np.savetxt(file_seq_name, arr3d[:,:,i], delimiter=' ', fmt='%d')
##################################################################
def plot_voxel(grey_arr, bins, outname):
    """
    plot_voxel: makes a 3d plot of the gray-scale 
    Parameters
    ----------
    grey_arr : 3d array
        DESCRIPTION.
    bins : integer
        DESCRIPTION. number of bins
    outname : string
        DESCRIPTION. plot file name

    Returns
    -------
    None.
    """
    gray_arr = grey_arr / 255
    #gray = np.ones((bins,bins,bins))-gray
    filled = np.ones((bins, bins, bins), dtype=bool)
    colors = np.repeat(gray_arr[:, :, :, np.newaxis], 3, axis=3)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.voxels(filled, facecolors=colors)
    plt.savefig(outname, dpi=300)
    plt.show()
 #############################################
# if __name__ == '__main__':
#     dump_file_arr = ['dump01.pos', 'dump02.pos','dump03.pos','dump04.pos','dump05.pos','dump06.pos']
#     name_seq_arr = ['rho01res15seq', 'rho02res15seq', 'rho03res15seq','rho04res15seq','rho05res15seq','rho06res15seq']
#     bin_arr = [52,41,36,33,30,28]
#     #################################
#     for i in range(6):
#         sim = NuclearLammps(dump_file_arr[i])
#         #sim.print_info()
#         Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(5)
#     #print('Proton fraction: {0:f}'.format(Yproton))
#     #plot_conf(xn, yn, zn, conf_name)
#         bins = bin_arr[i]
#         div, gray, vox = voxel(xp,yp,zp,sim.L_sim_box, bins, 1e13, name_seq_arr[i])
#         print('Done config {}'.format(dump_file_arr[i]))
#         print('***********************************')
    #print(vox)
    # kv, mv = minkowski(gray, sim.L_sim_box, bins)
    # print('***********************************')
    # print('Configuration: {}'.format(5))
    # print('Euler characteristic: {0:f}'.format(kv))
    # print('Average curvature: {0:f}'.format(mv))
    # print("You have made {}".format(type_of_pasta(kv, mv)))
    #res = np.array([div,div,div])
    #mink = mk.functionals(vox, res)
    # print('minkowski func: {}'.format(mink))
    #print('bin: {0}'.format(bin))
    #print('res: {0}'.format(div))
    #print(mink)
    #print('Volume:{0}'.format(mink[0]))
    #S = mink[1]*8
    #print('Surface: {0}'.format(S))
    #H = 2 * np.pi* np.pi * mink[2]
    #print('Integral mean curvature: {0}'.format(H))
    #chi = 4 * np.pi *mink[3] /3
    #print('Euler: {0}'.format(chi))
    #plot_voxel(gray, bins, gray_map)
    ########################
    #out_file = 'rho01-mink-Smallbin.dat'
    #ofile = open(out_file, 'w')
    #sim = NuclearLammps(dump_file_name)
    #sim.print_info()
    #Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(5)
    #bins = [31,39,53,79,158]
    #n_ths = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    #print('# n_th volume surface curvature euler')
    #bin = 300
    #n_th = 1e13
    #ofile.write('# bin \t vol \t surf \t curv euler \n ')
    #ofile.write('# n_th \t vol \t surf \t curv euler \n ')
    #div, testvoxels, charg = simple_voxel(xp,yp,zp,sim.L_sim_box, bins, n_th)
     #plot these as we did the points.
    #res = np.array([div,div,div])
    #mink = mk.functionals(charg, res)
    # #print('bin: {0}'.format(bin))
    # #print('res: {0}'.format(div))
    # #print(mink)
    # #print('Volume:{0}'.format(mink[0]))
    #S = mink[1]*8
    # #print('Surface: {0}'.format(S))
    #H = 2 * np.pi* np.pi * mink[2]
    # #print('Integral mean curvature: {0}'.format(H))
    #chi = 4 * np.pi *mink[3] /3
    # #print('Euler: {0}'.format(chi))
    # print('{0}\t  {1}\t  {2}\t {3}\t {4}'.format(n_th, mink[0], S, H, chi))

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(testvoxels[:][:,0],testvoxels[:][:,1],testvoxels[:][:,2])
    # ax.set_xlim3d(0, bins)
    # ax.set_ylim3d(0, bins)
    # ax.set_zlim3d(0, bins)
    # plt.show()
    # for n_th in n_ths:
    #for bin in bins:
    #     #bin = 200
    #    div, charg = simple_voxel(xp,yp,zp,sim.L_sim_box, bin, n_th)
    #    res = np.array([div,div,div])
    #    mink = mk.functionals(charg, res)
    #     #print('bin: {0}'.format(bin))
    #     #print('res: {0}'.format(div))
    #     #print(mink)
    #     #print('Volume:{0}'.format(mink[0]))
    #    S = mink[1]*8
    #     #print('Surface: {0}'.format(S))
    #    H = 2 * np.pi* np.pi * mink[2]
    #     #print('Integral mean curvature: {0}'.format(H))
    #    chi = 4 * np.pi *mink[3] /3
    #     #print('Euler: {0}'.format(chi))
        #print('{0}\t  {1}\t  {2}\t {3}\t {4}'.format(n_th, mink[0], S, H, chi))
    #    ofile.write('{0}\t  {1}\t  {2}\t {3}\t {4}\n'.format(bin, mink[0], S, H, chi))
    #ofile.close()
    ##########################
    ##########################
        #print('***********************')
    # for con in range(5):
    #     Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(con)
    #     #print('Proton fraction: {0:f}'.format(Yproton))
    #     #plot_conf(xn, yn, zn)
    #     bins = 50
    #     charg, gray = voxel(xp,yp,zp,sim.L_sim_box, bins)
    #     kv, mv = minkowski(gray, sim.L_sim_box, bins)
    #     print('***********************************')
    #     print('Configuration: {}'.format(con))
    #     print('Euler characteristic: {0:f}'.format(kv))
    #     print('Average curvature: {0:f}'.format(mv))
    #     print("You have made {}".format(type_of_pasta(kv, mv)))
    #     plot_voxel(gray, bins)
    ##########################
    ##########################
    # xp = np.array([10e-5, 40e-5])
    # yp = np.array([25e-5, 25e-5])
    # zp = np.array([25e-5, 25e-5])
    # Lb = 50e-5
    # bins = 100
    # n_th = 1e10
    # div, testvoxels = simple_voxel(xp,yp,zp,Lb, bins, n_th)
    # #plot these as we did the points.
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(testvoxels[:][:,0],testvoxels[:][:,1],testvoxels[:][:,2])
    # ax.set_xlim3d(0, bins)
    # ax.set_ylim3d(0, bins)
    # ax.set_zlim3d(0, bins)
    # plt.show()
