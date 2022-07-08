#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:24:13 2022
Program to generate slices of gray-scale map for different densities
@author: jbohorquez
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from itertools import product as iterate
from classNuclearLammps import *

# dump_file_name = 'dumpFile2sample/dump.pos-rho06'
# sim = NuclearLammps(dump_file_name)
# sim.print_info()

dump_file_arr = ['dumpFile2sample/dump.pos-rho01', \
                 'dumpFile2sample/dump.pos-rho02', \
                 'dumpFile2sample/dump.pos-rho03', \
                 'dumpFile2sample/dump.pos-rho04', \
                 'dumpFile2sample/dump.pos-rho05', \
                 'dumpFile2sample/dump.pos-rho06']
name_seq_arr = ['imagRes/rho01res15seq', 'imagRes/rho02res15seq',\
                'imagRes/rho03res15seq','imagRes/rho04res15seq',\
                'imagRes/rho05res15seq','imagRes/rho06res15seq']
bin_arr = [52,41,36,33,30,28]
conf_number = [60, 50, 25, 25, 24, 24]
#################################
for i in range(6):
    sim = NuclearLammps(dump_file_arr[i])
    #sim.print_info()
    Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(5)
#print('Proton fraction: {0:f}'.format(Yproton))
#plot_conf(xn, yn, zn, conf_name)
    bins = bin_arr[i]
    div, gray, vox = voxel(xp,yp,zp,sim.L_sim_box, bins, 1e13, name_seq_arr[i])
    print('Done config {}'.format(dump_file_arr[i]))
    print('***********************************')

