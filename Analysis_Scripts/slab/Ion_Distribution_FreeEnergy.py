#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import glob

path1          = './'
path2          = '/data/HOME_BACKUP/pengchao/deepmd2/lmp_dp/slab/lmp/DatasetWaterionSCAN/03/1000w_water_ion_slab0.005/lmp/300K/5500H2O_24OH/15-20ns'
find_ion_file1 =  'h3o_o_xyz_dump1_cut1.2_slab_4k_5-20ns.txt'
find_ion_file2 = 'oh_o_xyz_dump1_cut1.275_slab_4k_5-20ns.txt'
find_ion_file1 = os.path.join(path1, find_ion_file1)
find_ion_file2 = os.path.join(path2, find_ion_file2)
z_low          = 94.5
z_high         = 205.5
T              = 300 # K
txt1           = np.loadtxt(find_ion_file1)
txt2           = np.loadtxt(find_ion_file2)
time_step      = 1e-6 # ns
mintime        = txt1[0][1]*time_step
maxtime        = txt2[-1][1]*time_step

def get_o_list(txt):
    total_list = []
    for i in range(len(txt)):
        for j in range(int(txt[i][0])):
            total_list.append(txt[i][5+4*j])
    return total_list

def draw_ion_distribution(list1, list2,mintime,maxtime,
                          z_low,
                          z_high,
                          z_step = 500):
    bins = np.linspace(z_low, z_high, z_step)
    fig1 = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig1.add_subplot(1, 1, 1)
    
    y1, x1 = np.histogram(list1, bins)
    y2, x2 = np.histogram(list2, bins)
    tot_num = len(list1)
    ax.plot(x1[:-1], y1/tot_num*100, alpha=0.8, lw=2, label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    ax.plot(x2[:-1], y2/tot_num*100, alpha=0.8, lw=2, label='$\mathregular{OH^-}$')
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Ion distribution (%)", fontsize = 12)
    fig1.savefig(os.path.join('./', "Distribution_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

def draw_free_energy(list1, list2, T,mintime,maxtime,
                          z_low,
                          z_high,
                          z_step = 500):
    bins = np.linspace(z_low, z_high, z_step)
    fig = plt.figure(figsize=(13,13), dpi=150, facecolor='white')
    tot_num = len(list1)
    R = 8.314462618E-3 # kJ/(mol K)
    kJ2kcal = 0.239005736   
    # Free energy
    ax = fig.add_subplot(2, 2, 1)    
    y1, x1 = np.histogram(list1, bins)
    y2, x2 = np.histogram(list2, bins) 
    ax.plot(x1[:-1], np.log(y1/tot_num)*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    ax.plot(x2[:-1], np.log(y2/tot_num)*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{OH^-}$')
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Free energy (kcal/mol)", fontsize = 12)

    # Relative free energy
    ax = fig.add_subplot(2, 2, 2)    
    y1_avrg = np.mean(y1[int(z_step*0.25):int(z_step*0.75)])
    y2_avrg = np.mean(y2[int(z_step*0.25):int(z_step*0.75)])
    ax.plot(x1[:-1], (np.log(y1/tot_num)-np.log(y1_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    ax.plot(x2[:-1], (np.log(y2/tot_num)-np.log(y2_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{OH^-}$')
    #ax.set_xlim(98, 127)
    #ax.set_ylim(-1, 2)
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Relative free energy (kcal/mol)", fontsize = 12)
    
    # Relative free energy
    ax = fig.add_subplot(2, 2, 3)    
    ax.plot(x1[:-1], (np.log(y1/tot_num)-np.log(y1_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    ax.plot(x2[:-1], (np.log(y2/tot_num)-np.log(y2_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{OH^-}$')
    ax.set_xlim(98, 127)
    #ax.set_ylim(-1, 2)
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Relative free energy (kcal/mol)", fontsize = 12)
    
    # Relative free energy
    ax = fig.add_subplot(2, 2, 4)    
    ax.plot(x1[:-1], (np.log(y1/tot_num)-np.log(y1_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    ax.plot(x2[:-1], (np.log(y2/tot_num)-np.log(y2_avrg/tot_num))*R*kJ2kcal*T*(-1), alpha=0.8, lw=2, label='$\mathregular{OH^-}$')
    ax.set_xlim(173, 202)
    #ax.set_ylim(-1, 2)
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Relative free energy (kcal/mol)", fontsize = 12)    
    fig.savefig(os.path.join('./', "FreeEnergy_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

list1 = get_o_list(txt1)
list2 = get_o_list(txt2)

draw_ion_distribution(list1,list2,mintime,maxtime,z_low,z_high,z_step = 500)
draw_free_energy(list1, list2, T, mintime,maxtime,z_low, z_high, z_step = 500)