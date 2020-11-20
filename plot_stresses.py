import numpy as np
from generate_exp2 import *
import glob
from matplotlib.backends.backend_pdf import PdfPages
from loadmodel import *
from plotmodel import *
from GetAreas import *
from scipy.optimize import curve_fit
import gc
from mpl_toolkits.mplot3d import Axes3D
from bamg import *
from InterpFromMeshToMesh2d import *
import os
from copy import *
from plotting import *
from celluloid import Camera
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from basaldrag import *
from drivingstress import *
from slope import *
from matplotlib import colors

p='*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*.nc'
p='*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*asy*.nc'
p='*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*'
p='*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*off.nc'
p='*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*'
p='*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*'

modpath='./Models/embayments/'+p
mod=glue_runs_md(modpath)


fig = plt.figure()
gs = GridSpec(nrows=21, ncols=4)
triangles=mpl.tri.Triangulation(mod.mesh.x, mod.mesh.y, mod.mesh.elements-1)
x=(mod.mesh.x[mod.mesh.elements[:,0]-1]+mod.mesh.x[mod.mesh.elements[:,1]-1]+mod.mesh.x[mod.mesh.elements[:,2]-1])/3
y=(mod.mesh.y[mod.mesh.elements[:,0]-1]+mod.mesh.y[mod.mesh.elements[:,1]-1]+mod.mesh.y[mod.mesh.elements[:,2]-1])/3
countr=0
countc=0
vmin=-50000
vmax=50000
for i in range(205,225,1):
    ax=fig.add_subplot(gs[countr:countr+4, countc])
    if countr != 0:
        ax.set_xticklabels([])
    else:
        ax.xaxis.tick_top()
        xticks([35000,45000,55000,65000],[35,45,55,65])
    yticks([12000,15000,18000], [2,5,8])
    if countc != 0:
        ax.set_yticklabels([])
    countc+=1
    if countc==4:
        countr+=4
        countc=0
    P=mod.results.TransientSolution[i].Pressure
    ds=drivingstress(mod, mod.results.TransientSolution[i].Surface, mod.results.TransientSolution[i].Thickness)[-1]
    dsf=drivingstress(mod, mod.results.TransientSolution[i-1].Surface, mod.results.TransientSolution[i-1].Thickness)[-1]
    bd=basal_drag(mod, mod.results.TransientSolution[i].Thickness, mod.results.TransientSolution[i].Base, mod.results.TransientSolution[i].Vel)[0]
    bdf=basal_drag(mod, mod.results.TransientSolution[i-1].Thickness, mod.results.TransientSolution[i-1].Base, mod.results.TransientSolution[i-1].Vel)[0]
    #r=bd/(ds+1)
    #f=bdf/(dsf+1)
    ls=ds-bd
    lsf=dsf-bdf
    tripcolor(triangles, (bd).flatten(), vmin=vmin, vmax=vmax, cmap='RdYlBu_r')
    xlim(35000,70000)
    ylim(11000,19000)
    plotcontour(mod, mod.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors='green')
    plotcontour(mod, mod.results.TransientSolution[i].Thickness, levels=[10], colors='black')
ax21 = fig.add_subplot(gs[20,:])
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=mpl.cm.RdYlBu_r, norm=norm, label='Longitudinal Stress', orientation='horizontal')




fig = plt.figure()
gs = GridSpec(nrows=1, ncols=11)
triangles=mpl.tri.Triangulation(mod.mesh.x, mod.mesh.y, mod.mesh.elements-1)
vmin=-200000
vmax=200000
i=120    ### 120 for embayment, bump and depression, 12 for bottleneck
ds=drivingstress(mod, mod.results.TransientSolution[i].Surface, mod.results.TransientSolution[i].Thickness)[0]
bd=basal_drag(mod, mod.results.TransientSolution[i].Thickness, mod.results.TransientSolution[i].Base, mod.results.TransientSolution[i].Vx)[0]
ls=ds-bd
ax=fig.add_subplot(gs[:,1:10])
tripcolor(triangles, (ls).flatten(), vmin=vmin, vmax=vmax, cmap='RdYlBu_r')
xlim(35000,75000)
ylim(11000,19000)
plotcontour(mod, mod.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors='green', linewidths=2)
plotcontour(mod, mod.results.TransientSolution[i].Thickness, levels=[10], colors='black', linewidths=2)
yticks([12000,15000,18000], [2,5,8])
xticks([35000,45000,55000,65000,75000],[35,45,55,65,75])
xlabel('along-flow Distance [km]')
ylabel('across-flow Distance [km]')
ax21 = fig.add_subplot(gs[:,10])
norm = colors.Normalize(vmin=vmin/1000, vmax=vmax/1000)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=mpl.cm.RdYlBu_r, norm=norm, label='Longitudinal Stress [kPa]', orientation='vertical')
