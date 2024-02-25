#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 08:56:28 2023

Creates a Modflow6 groundwater model with 2 pumping wells using
a triangular mesh domain. The domain of this run is a single domain.
Plots are generated of the results.

@author: jrevier
"""

import os
import time

import numpy as np
import matplotlib.pyplot as plt

import flopy
from flopy.utils.triangle import Triangle as Triangle

starttime = time.time()
cpustart= time.process_time()

name = "tri_proj_nonest"
workspace = 'PATH_TO_SIMULATION_OUTPUT'
tri_exe = 'PATH_TO_triangle_EXECUTABLE'
mf6_exe = 'PATH_TO_mf6_EXECUTABLE'


#Creating the Aquifer Domain

#Step 1: Define the corner of each region
#ORDER MATTERS!!! BottomLeft, BottomRight, TopRight, TopLeft 
#(basically  Counter-Clockwise) if making a rectangular region

#full domain
active_domain = [(0,0), (1000,0),(1000,1000), (0,1000)]


#initialize Triangle Module to create the domain workspace
tri = Triangle(angle=25, model_ws= workspace, exe_name= tri_exe)

#add the region polygon to the workspace
tri.add_polygon(active_domain)


#add properties to the regions defined
#THE POINT THAT IS DEFINED MUST BE WITHIN THE REGION BEING ADDED!!!
#order added respective to order added to the workspace

tri.add_region((1,29), 0, maximum_area= 1109) #active domain

#build the domain
tri.build()


#create centeroid vertice plot of the domain
#this will be useful later when defining
#where the well is located within the domain

fig = plt.figure(figsize=(25, 15))
ax = plt.subplot(1,1,1, aspect='equal')
tri.plot(ax=ax, edgecolor='gray')
tri.plot_centroids(ax=ax, marker='o', color='red')
tri.label_cells(ax=ax, fontsize=10, color='red')
ax.set_xlabel('East to West (200m Intervals)', fontsize='20')
ax.set_ylabel('South to North (200m Intervals)', fontsize='20')
ax.set_title('Triangular Mesh Domain : 2 Wells', fontsize='20')
plt.show()


#building the MODFLOW Model
sim = flopy.mf6.MFSimulation(sim_name=name, version='mf6', 
                             exe_name=mf6_exe, 
                             sim_ws=workspace)

tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',nper=1, 
                             perioddata=[[1.0,1,1.]])

#call the Modflow Groundwater Flow Model
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)


#set head change calculations
ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', 
                           complexity='simple', outer_hclose=1.e-8, 
                           inner_hclose=1.e-8)


#model domain
cell2d = tri.get_cell2d()
vertices = tri.get_vertices()
xcyc = tri.get_xcyc()
ncpl = tri.ncpl
nvert = tri.nvert
ztop = 0.0
zbot = -150.0
nlay = 1
nrow = 30
ncol = 30
botm = 0


#setup the simulation using the .ModflowGwfdisv
#(Discretization with Vertices) for triangular mesh grids
dis = flopy.mf6.ModflowGwfdisv(gwf, length_units='METERS', nlay=nlay, 
                               ncpl=ncpl, 
                               nvert=nvert, top=ztop, botm=zbot,
                               vertices=vertices, cell2d=cell2d)

npf = flopy.mf6.ModflowGwfnpf(gwf, k=0.000009)

ic = flopy.mf6.ModflowGwfic(gwf, strt=0)

#create the constant head boundaries for the domain
chdList = []

leftcells= tri.get_edge_cells(4)
rightcells= tri.get_edge_cells(2)

for icpl in leftcells:
    chdList.append([(0,icpl), 100])

for icpl in rightcells:
    chdList.append([(0, icpl), 50])
    
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdList)

#create the wells

#For a grid that uses the DISV input file, 
#Well Location is (centeroid number)-1.

#see centeroid vertice plot created earlier
#CELLID will change if mesh shape is changed

#well1 location and pump rate
well1_loc = (1281) #centeroid vertice
well1_rate = -0.37

#create the well using ModflowGwfwel module
well1 = flopy.mf6.ModflowGwfwel(gwf, 
                              stress_period_data=[(0, 
                                                   well1_loc, 
                                                   well1_rate)])

#well2 location and pump rate 
well2_loc = (474)
well2_rate = -0.37
well2 = flopy.mf6.ModflowGwfwel(gwf, 
                              stress_period_data=[(0, 
                                                   well2_loc, 
                                                   well2_rate)])



#create the output file
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord='{}.cbc'.format(name),
                            head_filerecord='{}.hds'.format(name),
                            
                            saverecord=[('HEAD', 'LAST'),
                                        ('BUDGET', 'LAST')],
                            printrecord=[('HEAD', 'LAST'),
                                         ('BUDGET', 'LAST')])

sim.write_simulation()

success, buff = sim.run_simulation()






#create the plot (code from Hatari Labs)
fname = os.path.join(workspace, name + '.hds')
hdobj = flopy.utils.HeadFile(fname, precision='double')
head = hdobj.get_data()



from scipy.interpolate import griddata
x, y = tri.get_xcyc()[:,0],tri.get_xcyc()[:,1]
xG = np.linspace(x.min(),x.max(),100)
yG = np.linspace(y.min(),y.max(),100)
X, Y = np.meshgrid(xG, yG,)
z = head[0,0,]

fig = plt.figure(figsize=(15, 15))
ax = plt.subplot(1, 1, 1, aspect='equal')
img=tri.plot(ax=ax, edgecolor='lightsteelblue', 
             a=head[0, 0, :], 
             cmap='RdBu', vmin=-100, alpha=0.8)

colorbar = fig.colorbar(img, fraction=0.02)
Ti = griddata((x, y), z, (X, Y), method='cubic')
CS = ax.contour(X,Y,Ti,levels=np.arange(-150,100,15), 
                colors='Black', linewidths=2.0)

ax.clabel(CS, inline=1, fontsize=15, fmt='%1.0f')



ax.set_xlabel('East to West (200m Intervals)', fontsize='20')
ax.set_ylabel('South to North (200m Intervals)', fontsize='20')
ax.set_title('Triangular Mesh Domain: 2 Wells', fontsize='20')
colorbar.set_label('Head Elevations (m)', fontsize='15')

endtime = time.time()
cpuend= time.process_time()

plt.savefig('PATH_TO_....', 
            format='png')


computetime = cpuend - cpustart
time = endtime - starttime
print(ncpl)
print(computetime, time)
