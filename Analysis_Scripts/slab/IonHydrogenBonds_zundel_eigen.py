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
data_geo    = "5500H2O_2_24Hin_atomic_HOtype.data"
#data_geo    = "5500H2O_2_24Hin_2_atomic_HOtype.data" # 330K
trj_name    = "slab_4k_10-30ns.lammpstrj"
o_xyz       = "h3o_o_xyz_dump1_cut1.2_slab_4k_10-30ns.txt"
o_xyz_dump  = 1
density_file= 'Density_step1_dr0.1_slab_4k_10-30ns.txt'
edge_start  = -50
edge_end    = 10
trj_step    = 1
#300K h3o
interface   =  [102.27374921887291, 197.75251728473117] 
#330K h3o
#interface   =  [103.0267134849186, 197.00162580433377] 
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
oid_file    = os.path.join(trj_path,o_xyz)
trj_file    = os.path.join(trj_path,trj_name)
density_file= os.path.join(trj_path,density_file)


time_step   = 1e-6 # ns
oid_txt     = np.loadtxt(oid_file)
mintime     = oid_txt[0][1]*time_step
maxtime     = oid_txt[-1][1]*time_step

u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')

def get_mid_density(density_file):
    coord_list   = np.loadtxt(density_file,usecols=0)
    density_list = np.loadtxt(density_file,usecols=1)
    start = int(len(density_list)*0.25)
    end   = int(len(density_list)*0.75)
    mid_density = np.mean(density_list[start:end])/2
    #print(mid_density)
    for i in range(len(density_list)-1):
        if density_list[i] < mid_density and density_list[i+1] > mid_density:
            #print(i,coord_list[i],density_list[i],coord_list[i+1],density_list[i+1])
            interface1 = (mid_density-density_list[i])/(density_list[i+1]-density_list[i])*(coord_list[i+1]-coord_list[i])+coord_list[i]
        if density_list[i] > mid_density and density_list[i+1] < mid_density:
            #print(i,coord_list[i],density_list[i],coord_list[i+1],density_list[i+1])
            interface2 = (mid_density-density_list[i])/(density_list[i+1]-density_list[i])*(coord_list[i+1]-coord_list[i])+coord_list[i]
    
    return [interface1,interface2],mid_density

# set hbonds
def get_hbonds(u,index,scheme):
    if scheme == "donor":
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "index {0}".format(index), # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "type 2", # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )
    else:  # scheme == "acceptor"
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "type 2", # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "index {0}".format(index), # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )    
    return hbonds

def get_ion_hb_dict(u,o_id,scheme,frame,hb_dict,name):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_zcoord = u.atoms[o_id].position[2]
    o_d_or_a = hbonds.results.hbonds.shape[0]
    hb_dict[o_zcoord] = o_d_or_a
    #print(name,scheme,o_id,o_zcoord,o_d_or_a)    
    return hb_dict,o_d_or_a,hbonds.results.hbonds

def get_zundel(u,a_id,h_id,frame,zundel_dict,cutoff=1.3):
    """
    a_id: acceptor O
    """
    u.trajectory[frame]
    o_coord = u.atoms[a_id].position
    h_coord = u.atoms[h_id].position
    oh_distance = get_distance(o_coord, h_coord)
    if oh_distance <= cutoff:
        zundel_dict[h_coord[2]] = 1
    return zundel_dict

def get_eigen(u,eigen__dict,hbonds_results,cutoff=1.4):
    """
    d_id: donor O
    a_id: acceptor O
    """
    ah_distance = []
    for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
        u.trajectory[frame.astype(int)]
        a_coord = u.atoms[acceptor_index.astype(int)].position
        h_coord = u.atoms[hydrogen_index.astype(int)].position
        ah_distance.append(get_distance(a_coord, h_coord))
        
    d_coord = u.atoms[donor_index.astype(int)].position
    
    if ah_distance[0] > cutoff and ah_distance[1] > cutoff and ah_distance[2] > cutoff:
        eigen__dict[d_coord[2]] = 2
    return eigen__dict

def get_neg_hb_dict(u,o_id,scheme,frame,hb_dict,name):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_zcoord = u.atoms[o_id].position[2]
    o_d_or_a = hbonds.results.hbonds.shape[0]
    hb_dict[o_zcoord] = o_d_or_a
    #print(name,scheme,o_id,o_zcoord,o_d_or_a)    
    return hb_dict

def get_hist2d(hb_dict,edge_start,edge_end,mintime,maxtime,binss,name):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')    
    ax = fig.add_subplot(1, 1, 1)
    h = ax.hist2d(hb_dict.keys(),hb_dict.values(),bins=binss,
                  range=[[edge_start-2, edge_end+2],[-0.75, 7.5]],
                  cmap='Blues')
    fig.colorbar(h[3], shrink=1.0)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$",fontsize = 14)
    ax.set_ylabel("Number of adjacent water molecules",fontsize = 14)
    #fig.savefig(os.path.join('./', "HBNeg_%3s_%2d-%2dns.png"%(name,mintime,maxtime)), dpi=600, bbox_inches='tight')

def get_interface_hist2d(hb_dict,edge_start,edge_end,mintime,maxtime,binss,name,path,interface):
    fig = plt.figure(figsize=(15.8,6.5), dpi=150, facecolor='white')    
    ax = fig.add_subplot(121)
    h1 = ax.hist2d((list(hb_dict.keys())-interface[0])*(-1),hb_dict.values(),bins=binss,
                  range=[[-50, 10],[-0.75, 7.5]],
                  cmap='Blues')
    fig.colorbar(h1[3], shrink=1.0)
    ax.set_xlabel("Distance to lower interface "+r"$\ \rm (\AA)$",fontsize = 15)
    ax.set_ylabel("Number of adjacent water molecules",fontsize = 15)
    at = AnchoredText("(a)", prop=dict(size=15), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    
    ax = fig.add_subplot(122)
    h2 = ax.hist2d((list(hb_dict.keys())-interface[1]),hb_dict.values(),bins=binss,
                  range=[[-50, 10],[-0.75, 7.5]],
                  cmap='Blues')
    fig.colorbar(h2[3], shrink=1.0)
    ax.set_xlabel("Distance to upper interface "+r"$\ \rm (\AA)$",fontsize = 15)
    ax.set_ylabel("Number of adjacent water molecules",fontsize = 15)
    at = AnchoredText("(b)", prop=dict(size=15), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)      
    
    fig.savefig(os.path.join(path, "HBNeg_interface_%3s_%2d-%2dns.png"%(name,mintime,maxtime)), dpi=600, bbox_inches='tight')

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path) 

    if not isExists:
        os.makedirs(path) 

def write_dict(path,ion_dict,name,trj_step):
    file = open(os.path.join(path, "Configuration_{0}_step{1}.txt".format(name,trj_step)), 'w+')
    file.write('# %16s%16s\n' %('D(Angstrom)','Number'))
    for i in ion_dict.keys():
        file.write('%16.6f%16d\n' %(i,ion_dict[i]))
    file.close()   

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)

ion_donor_dict = {} # find ion as donor
zundel_dict    = {}
eigen__dict    = {}
for i in range(1,len(oid_txt),trj_step):
#for i in [1]:
    for j in range(int(oid_txt[i][0])):
        ion_o_id = int(oid_txt[i][2+4*j]-1)
        
        ion_donor_dict,ion_o_donors,hbonds_results = get_ion_hb_dict(u,ion_o_id,"donor",i*o_xyz_dump,ion_donor_dict,"ion_donor_dict")        
        if ion_o_donors > 0:
            for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
                ion_d_a = acceptor_index.astype(int)
                ion_d_h = hydrogen_index.astype(int)
                zundel_dict = get_zundel(u,ion_d_a,ion_d_h,i,zundel_dict)

        if ion_o_donors == 3:
            eigen__dict = get_eigen(u,eigen__dict,hbonds_results)

def zundel_eigen_hist2d(zundel_dict,eigen__dict,edge_start,edge_end,mintime,maxtime,path,interface,trj_step,binss=[150,5]):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')    
    
    ax = fig.add_subplot(211)
    x1 = []
    x2 = []
    for i in zundel_dict.keys():
        x1.append(i-interface)
    for i in eigen__dict.keys():
        x2.append(i-interface)
    
    h1 = ax.hist2d(x1,zundel_dict.values(),bins=binss,
                  range=[[edge_start, edge_end],[0,2]],
                  cmap='Blues')
    fig.colorbar(h1[3], shrink=1.0)
    #ax.set_xlabel("Distance to interface "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Zundel configuration",fontsize = 12)
    #ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_yticks([])
    at = AnchoredText("(a)", prop=dict(size=15), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)  
    
    ax = fig.add_subplot(212)
    h2 = ax.hist2d(x2,eigen__dict.values(),bins=binss,
                  range=[[edge_start, edge_end],[1,3]],
                  cmap='Blues')
    fig.colorbar(h2[3], shrink=1.0)
    ax.set_xlabel("Distance to interface "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Eigen configuration",fontsize = 12)
    #ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_yticks([])
    at = AnchoredText("(b)", prop=dict(size=15), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)  

    fig.savefig(os.path.join(path, "ZundelEigen_%4d_%2d-%2dns.png"%(trj_step,mintime,maxtime)), dpi=600, bbox_inches='tight')

def zundel_divide_eigen(zundel_dict,eigen__dict,edge_start,edge_end,mintime,maxtime,path,interface,trj_step):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')     
    ax = fig.add_subplot(111)
    z_low          = 94.5
    z_high         = 205.5
    delta_r        = 0.2
    bins   = np.linspace(z_low, z_high, int((z_high-z_low)/delta_r)+1)
    bin_centers = bins[:-1] + delta_r/2
    y1, x1 = np.histogram(list(zundel_dict.keys()), bins)
    y2, x2 = np.histogram(list(eigen__dict.keys()), bins)
    
    y1dy2 = np.full(bin_centers.size, fill_value=-1.0)
    for i in range(len(bin_centers)):
        if y1[i] != 0 and y2[i] != 0:
            y1dy2[i] = (y1[i]/y2[i]*100)
        
    file = open(os.path.join(path, "Zundel_divide_Eigen_step{0}.txt".format(trj_step)), 'w+')
    file.write('# %16s%16s\n' %('D(A)','Z/E(%)'))
    for i in range(len(bin_centers)):
        file.write('%16.6f%16.6f\n' %(bin_centers[i],y1dy2[i]))
    file.close()  
    
    ax.plot(bin_centers-interface, y1dy2, alpha=0.8, lw=2, c='dodgerblue')
    ax.set_xlabel("Distance to interface "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Zundel / Eigen (%)",fontsize = 12)
    ax.set_xlim(edge_start, edge_end)
    #at = AnchoredText("(b)", prop=dict(size=15), frameon=True, loc='upper left')
    #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    #ax.add_artist(at)  

    fig.savefig(os.path.join(path, "Zundel_divide_Eigen_%4d_%2d-%2dns.png"%(trj_step,mintime,maxtime)), dpi=600, bbox_inches='tight')

savepath = os.path.join(trj_path,"ZundelEigen")
mkdir(savepath)
zundel_eigen_hist2d(zundel_dict,eigen__dict,edge_start,edge_end,mintime,maxtime,savepath,interface[1],trj_step)

write_dict(savepath,zundel_dict,"Zundel",trj_step)
write_dict(savepath,eigen__dict,"Eigen",trj_step)
zundel_divide_eigen(zundel_dict,eigen__dict,edge_start,edge_end,mintime,maxtime,savepath,interface[1],trj_step)