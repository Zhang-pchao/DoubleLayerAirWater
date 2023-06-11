#!/usr/bin/env python
# coding: utf-8


import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.offsetbox import AnchoredText
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/deepmd2/lmp_dp/slab/lmp/datafiles"
trj_path    = "./"
data_geo    = "Sphere25_10H_atomic_HOtype.data"
#data_geo    = "Sphere25_10OH_atomic_HOtype.data"
#data_geo    = "Sphere25_eq_atomic.data"
trj_name    = "slab_1k.lammpstrj"
####################change above####################

if data_geo.split('_')[0] == "Sphere25":
    mda_step    = 4
    edge_start  = 0  
    edge_end    = 40  
    edge_step   = 0.2
    box_center  = 40
    xfig        = MultipleLocator(5)
else:        #"Sphere9"
    mda_step    = 500
    edge_start  = 0
    edge_end    = 20
    edge_step   = 1
    box_center  = 20
    xfig        = MultipleLocator(5)

HBfile = open(os.path.join(trj_path, "HBnumber_step{0}_edge{1}_{2}_perH2O.txt".format(mda_step,edge_step,trj_name.split('.')[0])), 'w+')

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

trj_skip    = 0 #int(len(u.trajectory)/10)   # skip 10% trajectorys
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

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)

for frame, donor_ix, *_ in hbonds.results.hbonds:
    u.trajectory[frame.astype(int)]
    donor = u.atoms[donor_ix.astype(int)]
    zpos = donor.position
    hist, *_ = np.histogram(get_distance(zpos, [box_center,box_center,box_center]), bins=bin_edges)
    # multiply by two as each hydrogen bond involves two water molecules  
    idx = int((frame-trj_skip)/mda_step)
    counts_list[int((frame-trj_skip)/mda_step)] += hist * 2

for i in hbonds.frames:
    u.trajectory[i]
    for k in range(len(u.atoms)):
        if u.atoms[k].type =='2':
            zpos2 = u.atoms[k].position
            #print(u.atoms[k].type)
            hist2, xxx = np.histogram(get_distance(zpos2, [box_center,box_center,box_center]), bins=bin_edges)
            #print(hist)
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

#ax.set_xlim(edge_start, edge_end)
#ax.set_ylim(-2, 32)

#at = AnchoredText("(a)", prop=dict(size=15), frameon=True, loc='upper left')
#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
#ax.add_artist(at)

ax.set_xlabel("Radial coordinates of donor atoms "+r"$\ \rm (\AA)$",fontsize = 15)
ax.set_ylabel("Number of hydrogen bonds per H"+'$\mathregular{_2O}$',fontsize = 15)
fig.savefig(os.path.join(trj_path, "HBnumber_step{0}_edge{1}_{2}_perH2O.png".format(mda_step,edge_step,trj_name.split('.')[0])), dpi=600)