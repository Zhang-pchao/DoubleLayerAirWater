#!/usr/bin/env python
# coding: utf-8


import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/deepmd2/lmp_dp/slab/lmp/datafiles"
trj_path    = "./"
data_geo    = "5500H2O_2_24Hin_atomic_HOtype.data"    # 300K
#data_geo    = "5500H2O_2_24Hin_2_atomic_HOtype.data" # 330K
#data_geo    = "5500H2O_2_24OHin_atomic_HOtype.data"
#data_geo    = "5500H2O_2_atomic_HOtype.data"
trj_name    = "slab_4k.lammpstrj"
scheme      = "4430"
####################change above####################

if scheme == "4412":
    mda_step    = 10
    edge_start  = 35  
    edge_end    = 85  
    edge_step   = 0.2
    xfig        = MultipleLocator(5)
else:        #"4430"
    mda_step    = 4
    edge_start  = 95
    edge_end    = 205
    edge_step   = 0.2
    xfig        = MultipleLocator(15)

HBfile = open(os.path.join(trj_path, "HBnumber_step{0}_edge{1}_{2}_preH2O.txt".format(mda_step,edge_step,trj_name.split('.')[0])), 'w+')

# bins in z for the histogram
bin_edges   = np.linspace(edge_start, edge_end, int((edge_end-edge_start)/edge_step)+1)
bin_centers = bin_edges[:-1] + edge_step/2
data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)

u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')

# set hbonds
hbonds = HydrogenBondAnalysis(
    universe           = u,
    donors_sel         = "type 2", # O
    hydrogens_sel      = "type 1", # H
    acceptors_sel      = "type 2", # O
    d_a_cutoff         = 3.5,      # <3.5
    d_h_a_angle_cutoff = 140,      # >140
    update_selections  = False
)

trj_skip    = 0#int(len(u.trajectory)/10)   # skip 10% trajectorys
hbonds.run(
    start   = trj_skip,                   # skip 10% trajectorys, default:None
    stop    = None,
    step    = mda_step,
    verbose = True
)

counts = np.full(bin_centers.size, fill_value=0.0)
counts_list = counts
for i in range(hbonds.n_frames-1):
    counts_list = np.vstack((counts_list, counts))

o_counts = np.full(bin_centers.size, fill_value=0.0)
o_counts_list = o_counts
for i in range(hbonds.n_frames-1):
    o_counts_list = np.vstack((o_counts_list, o_counts))

per_o_counts = np.full(bin_centers.size, fill_value=0.0)
per_o_counts_list = per_o_counts
for i in range(hbonds.n_frames-1):
    per_o_counts_list = np.vstack((per_o_counts_list, per_o_counts))

for frame, donor_ix, *_ in hbonds.results.hbonds:
    u.trajectory[frame.astype(int)]
    donor = u.atoms[donor_ix.astype(int)]
    zpos = donor.position[2]
    hist, *_ = np.histogram(zpos, bins=bin_edges)
    # multiply by two as each hydrogen bond involves two water molecules  
    idx = int((frame-trj_skip)/mda_step)
    counts_list[idx] += hist * 2

for i in hbonds.frames:
    u.trajectory[i]
    for k in range(len(u.atoms)):
        if u.atoms[k].type =='2':
            z_atom = u.atoms[k].position[2]
            #print(u.atoms[k].type)
            hist2, xxx = np.histogram(z_atom, bins=bin_edges)
            #print('hist2',hist2)
            idx = int((i-trj_skip)/mda_step)
            o_counts_list[idx] += hist2
    
    for l in range(len(per_o_counts_list[idx])):
        if o_counts_list[idx][l] == 0:
            per_o_counts_list[idx][l] = 0
        else:
            per_o_counts_list[idx][l]=counts_list[idx][l]/o_counts_list[idx][l]

counts_avg = np.average(per_o_counts_list, axis=0)
counts_std = np.std(per_o_counts_list, axis=0, ddof=1)
counts_up   = counts_avg + counts_std
counts_down = counts_avg - counts_std

HBfile.write('# %16s%16s%16s%16s%16s\n' %('Z(Angstrom)','HBnum_avrg','HBnum_up','HBnum_down','HBnum_std'))
for i in range(len(bin_centers)):
    HBfile.write('%16.4f%16.4f%16.4f%16.4f%16.4f\n' %(bin_centers[i],counts_avg[i],counts_up[i],counts_down[i],counts_std[i]))
HBfile.close()

# HBnumber_Fig
fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
ax = fig.add_subplot(1, 1, 1)

ax.plot(bin_centers, counts_avg, lw=2, label='Average value', c='firebrick')
ax.fill_between(bin_centers,counts_up,counts_down,facecolor = 'gold', alpha = 0.6, label='± Standard deviation')
ax.legend(fontsize = 12)

#ax.xaxis.set_major_locator(xfig)
#ax.yaxis.set_major_locator(MultipleLocator(5))

#ax.set_xlim(edge_start-2, edge_end+2)
#ax.set_ylim(-2, 32)

ax.set_xlabel("z-axis coordinates of donor atoms "+r"$\ \rm (\AA)$",fontsize = 14)
ax.set_ylabel("Number of hydrogen bonds per H"+'$\mathregular{_2O}$',fontsize = 14)
fig.savefig(os.path.join(trj_path, "HBnumber_step{0}_edge{1}_{2}_perH2O.png".format(mda_step,edge_step,trj_name.split('.')[0])), dpi=600)