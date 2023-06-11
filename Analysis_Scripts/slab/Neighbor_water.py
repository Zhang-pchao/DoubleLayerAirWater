#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import MDAnalysis as mda


####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/deepmd2/lmp_dp/slab/lmp/datafiles"
trj_path    = "./"
data_geo    = "5500H2O_2_24Hin_atomic_HOtype.data"
o_idx_file  = "h3o_o_xyz_dump1_cut1.2_slab_4k_5-25ns.txt"
trj_name    = "slab_4k_5-25ns.lammpstrj"
edge_start  = 95
edge_end    = 205
binss       = [150,15]
time_step   = 1e-6 # ns
trj_step    = 1
oo_cutoff   = 3.3
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
o_idx       = np.loadtxt(os.path.join(trj_path,o_idx_file))
mintime     = o_idx[0][1]*time_step
maxtime     = o_idx[-1][1]*time_step

u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')

oo_num_dict   = {}
oo_frm_dict   = {}
for j in range(0,len(u.trajectory),trj_step): # frames 
    u.trajectory[j]
    oo_idx_dict  = {}
    for k in range(int(o_idx[j][0])):         # o_index     
        target_o_idx = int(o_idx[j][2+k*4])-1
        target_o_z   = u.atoms[target_o_idx].position[2]
        o_neighbors  = 0
        oo_idx_list  = []
        for i in range(len(u.atoms)):         # o_index_neighbor
            if u.atoms[i].type == '2':
                if u.atoms[[target_o_idx,i]].bond.value() <= oo_cutoff:
                    o_neighbors += 1
                    oo_idx_list.append(i)
        oo_idx_dict[target_o_idx] = oo_idx_list
        oo_num_dict[target_o_z] = o_neighbors - 1 # minus target_o
    oo_frm_dict[j] = oo_idx_dict

#print(oo_num_dict)
#print(oo_frm_dict)

HBfile = open(os.path.join(trj_path, "Neighbor_%2d-%2dns.txt"%(mintime,maxtime)), 'w+')
HBfile.write('# %8s%16s%16s\n' %('frames','target_o_idx=','oo_idx_list'))
for j in range(0,len(u.trajectory),trj_step): # frames
    HBfile.write('%8d ' %(j))
    for i in oo_frm_dict[j].keys():
        HBfile.write('%8d=' %(i))
        for k in oo_frm_dict[j][i]:
            if k != i:
                HBfile.write('%8d' %(k))  
    HBfile.write('\n')             
HBfile.close()

# Fig
def get_scatter(oo_num_dict,edge_start,edge_end,trj_path,mintime,maxtime):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor="white")    
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(oo_num_dict.keys(),oo_num_dict.values())
    #ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.set_xlim(edge_start-2, edge_end+2)
    ax.set_ylim(min(oo_num_dict.values())-0.5, max(oo_num_dict.values())+0.5)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$",fontsize = 14)
    ax.set_ylabel("Number of adjacent water molecules",fontsize = 14)    
    fig.savefig(os.path.join(trj_path, "Neighbor_scatter_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

def get_hist2d(oo_num_dict,edge_start,edge_end,trj_path,mintime,maxtime,binss):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor="white")    
    ax = fig.add_subplot(1, 1, 1)
    h = ax.hist2d(oo_num_dict.keys(),oo_num_dict.values(),bins=binss,
                  range=[[edge_start-2, edge_end+2],[2.5, 7.5]],
                  cmap='Blues')
    fig.colorbar(h[3], shrink=1.0)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$",fontsize = 14)
    ax.set_ylabel("Number of adjacent water molecules",fontsize = 14)
    fig.savefig(os.path.join(trj_path, "Neighbor_hist2d_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

get_scatter(oo_num_dict,edge_start,edge_end,trj_path,mintime,maxtime)
get_hist2d(oo_num_dict,edge_start,edge_end,trj_path,mintime,maxtime,binss)