#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:29:42 2022

@author: jbohorquez
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from itertools import product as iterate
from classNuclearLammps import *

def db_file(dump_file,db_name,n_samples):
    sim = NuclearLammps(dump_file)
    init_idx = sim.num_conf - n_samples
    idxs = range(init_idx, sim.num_conf)
    
    db_out = open(db_name, 'w')
    #write first character in the header
    name_str = ''
    
    for i in idxs:
        Yproton, xp, yp, zp, xn, yn, zn = sim.read_conf(i)
        n_max_p = len(xp)
        print('number of protons: {}'.format(n_max_p))
        n_max_n = len(xn)
        print('number of neutrons: {}'.format(n_max_n))
        if i == init_idx:
            for j in range(n_max_p):
                name_str= name_str + 'xp{0},yp{0},zp{0},ipp{0},inp{0},'.format(j)
            for j in range(n_max_n):
                if j == n_max_n -1:
                    name_str= name_str + 'xn{0},yn{0},zn{0},ipn{0},inn{0}'.format(j)
                else:
                    name_str= name_str + 'xn{0},yn{0},zn{0},ipn{0},inn{0},'.format(j)
            # write the header in the file
            db_out.write(name_str+'\n')
        #Write proton coordinates
        isproton = 1
        isneutron = 0
        #print(len(xp))
        for j in range(n_max_p):
           db_out.write('{0},{1},{2},{3},{4},'.format(xp[j]/sim.L_sim_box, yp[j]/sim.L_sim_box, zp[j]/sim.L_sim_box, isproton, isneutron))
        #write neutron coordinates
        isproton = 0
        isneutron = 1    
        for j in range(n_max_n):
            if j == n_max_n-1:
                db_out.write('{0},{1},{2},{3},{4}\n'.format(xn[j]/sim.L_sim_box, yn[j]/sim.L_sim_box, zn[j]/sim.L_sim_box, isproton, isneutron))
            else:
                db_out.write('{0},{1},{2},{3},{4},'.format(xn[j]/sim.L_sim_box, yn[j]/sim.L_sim_box, zn[j]/sim.L_sim_box, isproton, isneutron))
    db_out.close()
    

dump_names = ['dumpFile2sample/dump.pos-rho01',\
            'dumpFile2sample/dump.pos-rho02',\
                'dumpFile2sample/dump.pos-rho03',\
                    'dumpFile2sample/dump.pos-rho04',\
                        'dumpFile2sample/dump.pos-rho05',\
                            'dumpFile2sample/dump.pos-rho06']
    
db_names = ['rho01Sam2_db.csv','rho02Sam2_db.csv','rho03Sam2_db.csv','rho04Sam2_db.csv','rho05Sam2_db.csv','rho06Sam2_db.csv']
    
n_samp = 10

for m in range(6):
    db_file(dump_names[m] , db_names[m], n_samp)