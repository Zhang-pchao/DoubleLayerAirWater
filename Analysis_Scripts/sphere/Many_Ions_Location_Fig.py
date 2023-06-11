#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.colors 
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

def ManyIons_3D_fig(z_major, zlim, lable_s, x, y, z, color,max_T,min_T,name):
    fig = plt.figure(figsize=(8,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(111,projection='3d')
    plt.set_cmap(plt.get_cmap("viridis", 100))
    fmt = matplotlib.ticker.FuncFormatter(lambda x,pos:'%.2f' %((x*(max_T-min_T)+min_T)))
    for i in range(len(x)):
        im = ax.scatter(x[i], y[i], z[i], s=lable_s,c=color[i],marker='.')
    fig.colorbar(im, shrink=0.7,format = fmt,label='Simulation time (ns)')            
    
    ax.set_xlabel('X'+r"$\ \rm (\AA)$", fontsize = 10)
    ax.set_ylabel('Y'+r"$\ \rm (\AA)$", fontsize = 10)
    ax.set_zlabel('Z'+r"$\ \rm (\AA)$", fontsize = 10)
    
    # coordinate spacing
    x_major_locator = MultipleLocator(5)
    y_major_locator = MultipleLocator(5)
    z_major_locator = MultipleLocator(z_major)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    ax.zaxis.set_major_locator(z_major_locator)
    
    ax.set_xlim(-0.5, 41)
    ax.set_ylim(-0.5, 41)
    ax.set_zlim(zlim[0], zlim[1])

    fig.savefig(os.path.join('./', "O_%s_%.2fns.png" %(name,max_T)), dpi=600, bbox_inches='tight')

def ManyIons_2D_fig(z_major, zlim, lable_s, X,Y,color,max_T,min_T,name):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    plt.set_cmap(plt.get_cmap("viridis", 100)) 
    for i in range(len(X)):
        for yi in Y[i]:
            ax.scatter(X[i], yi, s=lable_s,c=color[i],marker='.')    
    
    if name.split('_')[0] == 'Z':
        ax.set_ylim(zlim[0], zlim[1])
        y_major_locator = MultipleLocator(z_major)
        ax.yaxis.set_major_locator(y_major_locator)
    else: #name == 'X' or 'Y':
        ax.set_ylim(-0.5, 41)
        y_major_locator = MultipleLocator(5)
        ax.yaxis.set_major_locator(y_major_locator)
    
    ax.set_xlabel("Simulation time (ns)", fontsize = 12)
    ax.set_ylabel("%s-axis coordinates of ions " %name.split('_')[0]+r"$\ \rm (\AA)$", fontsize = 12)    
    fig.savefig(os.path.join('./', "O_%s_%.2fns.png" %(name,max_T)), dpi=600, bbox_inches='tight')

def OneIon_3D_fig(z_major, zlim, lable_s, x, y, z, color, max_T, min_T, name):
    fig = plt.figure(figsize=(8,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(111,projection='3d')
    plt.set_cmap(plt.get_cmap("viridis", 100))
    im = ax.scatter(x, y, z, s=lable_s,c=color,marker='.')
    
    # colorbar
    fmt = matplotlib.ticker.FuncFormatter(lambda x,pos:'%.2f' %(x*(max_T-min_T)+min_T))
    fig.colorbar(im, shrink=0.7,format = fmt,label='Simulation time (ns)')
    
    ax.set_xlabel('X'+r"$\ \rm (\AA)$", fontsize = 10)
    ax.set_ylabel('Y'+r"$\ \rm (\AA)$", fontsize = 10)
    ax.set_zlabel('Z'+r"$\ \rm (\AA)$", fontsize = 10)
    
    # coordinate spacing
    x_major_locator = MultipleLocator(5)
    y_major_locator = MultipleLocator(5)
    z_major_locator = MultipleLocator(z_major)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    ax.zaxis.set_major_locator(z_major_locator)
    
    ax.set_xlim(-0.5, 41)
    ax.set_ylim(-0.5, 41)
    ax.set_zlim(zlim[0], zlim[1])

    fig.savefig(os.path.join('./', "O_%s_%.2fns.png" %(name,max_T)), dpi=600, bbox_inches='tight')

def OneIon_2D_fig(z_major, zlim, lable_s, x, y, color, max_T, min_T, name):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    plt.set_cmap(plt.get_cmap("viridis", 100))
    im = ax.scatter(x, y, s=lable_s,c=color,marker='.')
        
    ax.set_xlabel("Simulation time (ns)", fontsize = 12)
    ax.set_ylabel("%s-axis coordinates of ions " %name.split('_')[0]+r"$\ \rm (\AA)$", fontsize = 12)    
    
    if name.split('_')[0] == 'Z':
        ax.set_ylim(zlim[0], zlim[1])
        y_major_locator = MultipleLocator(z_major)
        ax.yaxis.set_major_locator(y_major_locator)
    else: #name == 'X' or 'Y':
        ax.set_ylim(-0.5, 41)
        y_major_locator = MultipleLocator(5)
        ax.yaxis.set_major_locator(y_major_locator)

    fig.savefig(os.path.join('./', "O_%s_%.2fns.png" %(name,max_T)), dpi=600, bbox_inches='tight')

def ManyIons_2D_fig_hist2d(X,Y,zlim,binss,name):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    xxx = []
    yyy = []
    for i in range(len(X)):
        for yi in Y[i]:
            xxx.append(X[i])
            yyy.append(yi)
    h = ax.hist2d(xxx, yyy, bins=binss,range=[[min(xxx)-0.3,max(xxx)+0.3], [zlim[0], zlim[1]]],cmap='Blues')
    fig.colorbar(h[3], shrink=1.0)
    #ax.yaxis.set_major_locator(MultipleLocator(z_major))    
    ax.set_xlabel("Simulation time (ns)", fontsize = 12)
    ax.set_ylabel("Z-axis coordinates of ions"+r"$\ \rm (\AA)$", fontsize = 12)    
    fig.savefig(os.path.join('./', "Distribution_%3s_hist2d_%2d-%2dns.png"%(name,min(xxx),max(xxx))), dpi=600, bbox_inches='tight')

def Sphere_2D_fig_hist2d(T,X,Y,Z,center,binss,name):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    xxx = []
    yyy = []
    for i in range(len(T)):
        for j in range(len(X[i])):
            d = np.sqrt(get_distance([X[i][j],Y[i][j],Z[i][j]], center))
            xxx.append(T[i])
            yyy.append(d)
    h = ax.hist2d(xxx, yyy, bins=binss,range=[[min(xxx)-0.3,max(xxx)+0.3], [0.0, 30.0]],cmap='Blues')
    fig.colorbar(h[3], shrink=1.0)
    #ax.yaxis.set_major_locator(MultipleLocator(z_major))    
    ax.set_xlabel("Simulation time (ns)", fontsize = 12)
    ax.set_ylabel("Distances between ions and the center of droplet "+r"$\ \rm (\AA)$", fontsize = 12)    
    fig.savefig(os.path.join('./', "Distribution_%3s_sphere_%2d-%2dns.png"%(name,min(xxx),max(xxx))), dpi=600, bbox_inches='tight')

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return distance


if __name__ == '__main__':

    path          = './'
    find_ion_file = 'h3o_o_xyz_dump1_cut1.2_slab_1k.txt'
    ion_file_name = find_ion_file.split("_")[0]
    find_ion_file = os.path.join(path, find_ion_file)
    time_step     = 1e-6 # ns
    ion_num       = -2 # -1 for slab hist; -2 for sphere; 1 for one ion; >1 for scatter figs
    sphere_center = [40.0, 40.0, 40.0]
    zlim          = [94.5, 205.5] # [34.5, 85.5]
    lable_s       = 20
    z_major       = 20
    binss         = 100
    txt           = np.loadtxt(find_ion_file)
    
    T = []
    X = []
    Y = []
    Z = []
    if ion_num == 1:        
        for i in range(len(txt)):
            if int(txt[i][0]) == 1:
                T.append(txt[i][1]*time_step)
                X.append(txt[i][3])
                Y.append(txt[i][4])
                Z.append(txt[i][5])
            elif int(txt[i][0]) == 2:
                oo_distance = get_distance([txt[i][3],txt[i][4],txt[i][5]], [txt[i][7],txt[i][8],txt[i][9]])
                if oo_distance < 12:
                    T.append(txt[i][1]*time_step)
                    X.append((txt[i][3]+txt[i][7])/2)
                    Y.append((txt[i][4]+txt[i][8])/2)
                    Z.append((txt[i][5]+txt[i][9])/2)
            elif int(txt[i][0]) > 2:                
                print("Find more than 2 ions after step {0}".format(txt[i-1][1]))
            else:
                print("Find 0 ions after step {0}".format(txt[i-1][1]))
    
    else: # ion_num > 1
        for i in range(len(txt)):
            T.append(txt[i][1]*time_step)
            x = []
            y = []
            z = []
            for j in range(int(txt[i][0])):
                x.append(txt[i][3+4*j])
                y.append(txt[i][4+4*j])
                z.append(txt[i][5+4*j])
            X.append(x)
            Y.append(y)
            Z.append(z)
    
    min_T = min(T)
    max_T = max(T)
    color = [plt.get_cmap("viridis", 100)(int(float(i-min_T)/(max_T-min_T)*100)) for i in T]
    
    if ion_num == 1:
        OneIon_2D_fig(z_major, zlim, lable_s, T,X,color,max_T,min_T,'X_'+ion_file_name)
        OneIon_2D_fig(z_major, zlim, lable_s, T,Y,color,max_T,min_T,'Y_'+ion_file_name)
        OneIon_2D_fig(z_major, zlim, lable_s, T,Z,color,max_T,min_T,'Z_'+ion_file_name)
        OneIon_3D_fig(z_major, zlim, lable_s, X,Y,Z,color,max_T,min_T,'XYZ_'+ion_file_name)
    elif ion_num == -1:
        ManyIons_2D_fig_hist2d(T,Z,zlim,binss,ion_file_name)
    elif ion_num == -2:    
        Sphere_2D_fig_hist2d(T,X,Y,Z,sphere_center,binss,ion_file_name)
    else: # ion_num > 1
        ManyIons_2D_fig(z_major, zlim, lable_s, T,X,color,max_T,min_T,'X_'+ion_file_name)
        ManyIons_2D_fig(z_major, zlim, lable_s, T,Y,color,max_T,min_T,'Y_'+ion_file_name)
        ManyIons_2D_fig(z_major, zlim, lable_s, T,Z,color,max_T,min_T,'Z_'+ion_file_name)
        ManyIons_3D_fig(z_major, zlim, lable_s, X,Y,Z,color,max_T,min_T,'XYZ_'+ion_file_name) 