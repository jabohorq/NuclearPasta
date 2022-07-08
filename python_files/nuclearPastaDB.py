#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:21:54 2022

@author: jbohorquez
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from itertools import product as iterate
from classNuclearLammps import *

dump_file_name = 'dumpFile1sample/dump.pos-rho01'
out_db_name = 'nuclear_db_test.csv'

out_db = open(out_db_name, 'w')

sim = NuclearLammps(dump_file_name)
#sim.print_info()


configs = range(2)
for m in configs:
    Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(m)
    #print(m)
    #sim.print_info()
    isproton = 1
    isneutron = 0
    #print(len(xp))
    for i in range(len(xp)):
        #print(xp[i]/sim.L_sim_box, yp[i]/sim.L_sim_box, zp[i]/sim.L_sim_box, isproton, isneutron)
        out_db.write('{0},{1},{2},{3},{4},'.format(xp[i]/sim.L_sim_box, yp[i]/sim.L_sim_box, zp[i]/sim.L_sim_box, isproton, isneutron))
    isproton = 0
    isneutron = 1    
    for i in range(len(xn)):
        #print(xn[i]/sim.L_sim_box, yn[i]/sim.L_sim_box, zn[i]/sim.L_sim_box, isproton, isneutron)
        if i == len(xn)-1:
            out_db.write('{0},{1},{2},{3},{4}\n'.format(xn[i]/sim.L_sim_box, yn[i]/sim.L_sim_box, zn[i]/sim.L_sim_box, isproton, isneutron))
        else:
            out_db.write('{0},{1},{2},{3},{4},'.format(xn[i]/sim.L_sim_box, yn[i]/sim.L_sim_box, zn[i]/sim.L_sim_box, isproton, isneutron))
out_db.close()