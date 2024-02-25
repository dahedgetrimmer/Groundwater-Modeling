#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 21:32:44 2023

Creates a Modflow6 groundwater model with 2 pumping wells using
a gridded domain. Plots are generated of the results.

@author: jrevier
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import flopy

starttime = time.time()
cpustart= time.process_time()

workspace = 'PATH_TO_SIMULATION_OUTPUTS'
mf6_exe = 'PATH_TO_mf6_EXECUTABLE'
name = "gridded"
sim = flopy.mf6.MFSimulation(name, 
                           exe_name=mf6_exe, 
                           )


# Domain parameters
Lx = 1000.0
Ly = 1000.0
ztop = 0.0
zbot = -150.0
nlay = 1
nrow = 30
ncol = 30
delr = Lx / ncol
delc = Ly / nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)

# Create a MODFLOW 6 model and simulation

sim = flopy.mf6.MFSimulation(sim_name=name,version='mf6', exe_name= mf6_exe, sim_ws=workspace)
tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',nper=1, perioddata=[[1.0,1,1.]])
ims = flopy.mf6.ModflowIms(sim, pname="ims", complexity="SIMPLE")

# Create the MODFLOW 6 groundwater flow model (GWF)
model_nam_file = "{}.nam".format(name)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=model_nam_file, save_flows=True)

# Create the discretization package (DIS) with StructuredGrid
dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc,
                              top=ztop, botm=botm[1:])

npf = flopy.mf6.ModflowGwfnpf(gwf, k=0.000009)

# Set up the constant head 
chdspd = []

# Loop through rows to set the constant head boundary on the left side
for i in range(nrow):
    chdspd.append([(0, i, 0), 100])  # [(layer, row, column), head value]
    chdspd.append([(0,i,29),50])

# Create the constant head package (CHD)
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

# Set up initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)  # Initial conditions with zero head


#create the wells
well_loc = (0,15, 15)  # (layer, row, column)
pumping_rate = -0.37  # Negative value indicates pumping


well2_loc = (0, 4, 24)  # (layer, row, column)
pumping_rate = -0.37  # Negative value indicates pumping


wel = flopy.mf6.ModflowGwfwel(gwf, 
                              stress_period_data=[(well_loc, pumping_rate)])

wel = flopy.mf6.ModflowGwfwel(gwf, 
                              stress_period_data=[(well2_loc, pumping_rate)])

# Create the output control (`OC`) Package
headfile = "{}.hds".format(name)
head_filerecord = [headfile]
budgetfile = "{}.cbc".format(name)
budget_filerecord = [budgetfile]
saverecord = [("HEAD", "LAST"), ("BUDGET", "LAST")]
printrecord = [("HEAD", "LAST")]
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    saverecord=saverecord,
    head_filerecord=head_filerecord,
    budget_filerecord=budget_filerecord,
    printrecord=printrecord,
)
# Write the MODFLOW 6 input files
sim.write_simulation()

# Run the MODFLOW 6 simulation
success, buff = sim.run_simulation()



   #create plots
import flopy.utils.binaryfile as bf

hds = bf.HeadFile(workspace + "gridded.hds")
head = hds.get_data(totim=1.0)


# Contour the heads

extent = (delr / 2.0, Lx - delr / 2.0, Ly - delc / 2.0, delc / 2.0)
mapview = flopy.plot.PlotMapView(model=gwf, layer=0)
quadmesh = mapview.plot_ibound()
linecollection = mapview.plot_grid()
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
contour = ax.contour(head[0, :, :], levels=np.arange(-150, 100, 15), extent=extent, cmap='RdBu')
ax.set_xlabel('East to West (200m Intervals)', fontsize='15')
ax.set_ylabel('South to North (200m Intervals)', fontsize='15')
ax.set_title('Grid Domain: 2 Wells', fontsize='15')
ax.clabel(contour, inline=1, fontsize=14, fmt='%1.0f')
plt.grid(True)
plt.savefig('PATH_TO_.....', format='png')


# Plot color grid of heads
extent = (delr / 2.0, Lx - delr / 2.0, Ly - delc / 2.0, delc / 2.0)
mapview = flopy.plot.PlotMapView(model=gwf, layer=0)
quadmesh = mapview.plot_ibound()
linecollection = mapview.plot_grid()
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
plt.imshow(head[0, :, :], extent=extent, cmap='RdBu', vmin=-100)
colorbar = plt.colorbar(fraction=0.035, orientation='horizontal')
ax.set_xlabel('East to West (200m Intervals)', fontsize='15')
ax.set_ylabel('South to North (200m Intervals)', fontsize='15')
ax.set_title('Grid Domain: 2 Wells', fontsize='15')
colorbar.set_label('Head Elevation (m)')
plt.savefig('PATH_TO_.....', format='png')

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
mapview = flopy.plot.PlotMapView(model=gwf, layer=0)
quadmesh = mapview.plot_ibound()
linecollection = mapview.plot_grid()
contour_set = mapview.contour_array(head, 
                                    levels=np.arange(-150, 100, 15), 
                                    colors="r")

plt.clabel(contour_set, inline=1, fontsize=14, fmt='%1.0f');
ax.set_xlabel('East to West (200m Intervals)', fontsize='15')
ax.set_ylabel('South to North (200m Intervals)', fontsize='15')
ax.set_title('Grid Domain: 2 Wells', fontsize='15')


endtime = time.time()
cpuend= time.process_time()

plt.savefig('PATH_TO_.....')

computetime = cpuend - cpustart
time = endtime - starttime

print(computetime, time)
