#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.colors 
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator


path       = './'
log_file   = 'log.lammps'
ion_file   = 'h3o_coord_dump1_cut1.2_slab_4k.txt'
middle_z = 150.0
time_step = 1e-6 # ns

log_file   = os.path.join(path, log_file)
ion_file   = os.path.join(path, ion_file)
log_info   = np.loadtxt(log_file)
ion_info   = np.loadtxt(ion_file)

z_info_dict = {}
e_info_dict = {}

for i in range(1,ion_info.shape[0]):
    plus_z = 0.0
    for j in range(int(ion_info[i][0])):
        plus_z += np.square(ion_info[i][5+4*j]-middle_z)
    plus_z = np.sqrt(plus_z/ion_info[i][0])
    z_info_dict[int(ion_info[i][1])] = plus_z

min_t = min(z_info_dict.keys())*time_step
max_t = max(z_info_dict.keys())*time_step
color = [plt.get_cmap("viridis", 100)(int(float(i*time_step-min_t)/(max_t-min_t)*100)) for i in z_info_dict.keys()]

for i in range(log_info.shape[0]):
    if int(log_info[i][0]) in z_info_dict.keys():
        e_info_dict[z_info_dict[int(log_info[i][0])]] = log_info[i][2]

fig = plt.figure(figsize=(8,6.5), dpi=150, facecolor='white')
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("The distance of ions from the center"+r"$\ \rm (\AA)$", fontsize = 12)
ax.set_ylabel("Potential energy (eV)", fontsize = 12)
im = ax.scatter(e_info_dict.keys(), e_info_dict.values(), s=100, c=color, marker='.')    
# colorbar
fmt = matplotlib.ticker.FuncFormatter(lambda x,pos:'%.2f' %(x*(max_t-min_t)+min_t))
fig.colorbar(im, shrink=1.0,format = fmt,label='Simulation time (ns)')    
#fig.show()
fig.savefig(os.path.join(path, "Energy_vs_Distribution.png"), dpi=600)