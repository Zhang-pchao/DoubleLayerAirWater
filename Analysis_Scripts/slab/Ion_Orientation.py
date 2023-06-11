#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import glob

path           = './'
find_ion_file  = 'h3o_h_xyz_dump1_cut1.2_slab_4k_5-25ns.txt'
scheme         = find_ion_file.split('_')[0]
find_ion_file  = os.path.join(path, find_ion_file)
if scheme == 'h3o':
    binss = 150
else:
    binss = 150
z_low          = 94.5
z_high         = 205.5
time_step      = 1e-6 # ns
txt            = np.loadtxt(find_ion_file)
mintime        = txt[0][1]*time_step
maxtime        = txt[-1][1]*time_step

def get_coord(txt,scheme):
    h_to_o_dict = {}    
    for i in range(len(txt)):
        for j in range(int(txt[i][0])):
            if scheme == 'h3o':
                o_coord=[ txt[i][3+16*j], txt[i][4+16*j], txt[i][5+16*j]]
                h_coord=[[txt[i][7+16*j], txt[i][8+16*j], txt[i][9+16*j]],
                        [txt[i][11+16*j],txt[i][12+16*j],txt[i][13+16*j]],
                        [txt[i][15+16*j],txt[i][16+16*j],txt[i][17+16*j]]]
                center = getCenterOfCircle(h_coord)
            else: # scheme == 'oh'
                o_coord=[txt[i][3+8*j], txt[i][4+8*j], txt[i][5+8*j]]
                center =[txt[i][7+8*j], txt[i][8+8*j], txt[i][9+8*j]] # center = h_coord
            h_to_o = [o_coord[0]-center[0],o_coord[1]-center[1],o_coord[2]-center[2]]
            cos    = getCOS(h_to_o,[0,0,1])
            h_to_o_dict[o_coord[2]] = cos                
    return h_to_o_dict

def draw_ion_orientation(mintime,maxtime,h_to_o_dict,scheme):
    fig  = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    if scheme == 'h3o':
        ax.scatter(h_to_o_dict.keys(), h_to_o_dict.values(),s=5,label='$\mathregular{H_3}$'+'$\mathregular{O^+}$')
    else: #'oh'
        ax.scatter(h_to_o_dict.keys(), h_to_o_dict.values(),s=5,label='$\mathregular{OH^-}$')
    ax.legend(fontsize = 12)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("cos "+r'$\theta$', fontsize = 12)
    fig.savefig(os.path.join('./', "Orientation_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

def draw_ion_orientation_hist2d(mintime,maxtime,h_to_o_dict,binss,scheme):
    fig  = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    if scheme == 'h3o':
        h = ax.hist2d(h_to_o_dict.keys(), h_to_o_dict.values(), bins=binss,range=[[97,203], [-1.05,1.05]],cmap='Blues')
    else: # 'oh'
        h = ax.hist2d(h_to_o_dict.keys(), h_to_o_dict.values(), bins=binss,range=[[97,203], [-1.05,1.05]],cmap='Blues')    
    fig.colorbar(h[3], shrink=1.0)
    ax.set_xlabel("Z-axis coordinates of ions "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("cos "+r'$\theta$', fontsize = 12)
    fig.savefig(os.path.join('./', "Orientation_hist2d_%2d-%2dns.png"%(mintime,maxtime)), dpi=600, bbox_inches='tight')

def getCenterOfCircle(points):
    x1 = points[0][0]
    x2 = points[1][0]
    x3 = points[2][0]
    y1 = points[0][1]
    y2 = points[1][1]
    y3 = points[2][1]
    z1 = points[0][2]
    z2 = points[1][2]
    z3 = points[2][2]
    
    a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2)
    b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2)
    c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
    d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1)
    
    a2 = 2 * (x2 - x1)
    b2 = 2 * (y2 - y1)
    c2 = 2 * (z2 - z1)
    d2 = x1*x1 + y1*y1 + z1*z1 - x2*x2 - y2*y2 - z2*z2

    a3 = 2 * (x3 - x1)
    b3 = 2 * (y3 - y1)
    c3 = 2 * (z3 - z1)
    d3 = x1*x1 + y1*y1 + z1*z1 - x3*x3 - y3*y3 - z3*z3

    cx = -(b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1)/(a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)
    cy =  (a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1)/(a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)
    cz = -(a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1)/(a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)
    
    return [cx,cy,cz]

def getCOS(x,y):
    x = np.asarray(x)
    y = np.asarray(y)
    l_x=np.sqrt(x.dot(x))
    l_y=np.sqrt(y.dot(y))
    dian=x.dot(y)
    cos_=dian/(l_x*l_y)
    return cos_

h_to_o_dict = get_coord(txt,scheme)
draw_ion_orientation(mintime,maxtime,h_to_o_dict,scheme)
draw_ion_orientation_hist2d(mintime,maxtime,h_to_o_dict,binss,scheme)