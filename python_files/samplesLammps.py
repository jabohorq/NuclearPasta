#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:54:30 2022

@author: jbohorquez
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from itertools import product as iterate
from classNuclearLammps import *

files_names= ['dumpFile1sample/dump.pos-rho01',\
              'dumpFile1sample/dump.pos-rho02',\
              'dumpFile1sample/dump.pos-rho03',\
                  'dumpFile1sample/dump.pos-rho04',\
                      'dumpFile1sample/dump.pos-rho05',\
                          'dumpFile1sample/dump.pos-rho06']

dens = ['01','02','03','04','05','06']
bins = [52,41,36,33,30,28]
#dump_file_name = 'dumpFile2sample/dump.pos-rho01'

for m in range(6):
    sim = NuclearLammps(files_names[m])
    init_idx = sim.num_conf - 10
    idxs = range(init_idx, sim.num_conf)
    name_seq = 'ML-data/sample1/rho'+dens[m]+'/rho'+dens[m]+'sam'
    # #################################
    for i in idxs:
        file_seq_name = name_seq + '%02d-slc'% i 
        #sim = NuclearLammps(dump_file_name)
        Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(i)
        #bins = 52
        div, gray, vox = voxel(xp,yp,zp,sim.L_sim_box, bins[m],\
                           1e13, file_seq_name)
        print('Done config {}'.format(i))
        print('***********************************')

