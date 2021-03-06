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

fs=16
ls=12

md=loadmodel('Models/GeomProj_SpinUp_load_SG_spcvx50_NL-450_FrM200_FlMreal30_FC40_FT100_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT875_inH0_smb55_smbPos30000_funnel300_FullMelt.nc')
md2=loadmodel('Models/embayments/GeomProj_extenddomain_SG_spcvx50_NL-450_FrM1200_FlMreal180_FC40_FT300_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH1350_ByP55000_ByS20000_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT1575_inH0_smb55_smbPos30000_funnel300_FullMelt_noCutoff.nc')

triangles=mpl.tri.Triangulation(md.mesh.x, md.mesh.y, md.mesh.elements-1)
x=(md.mesh.x[md.mesh.elements[:,0]-1]+md.mesh.x[md.mesh.elements[:,1]-1]+md.mesh.x[md.mesh.elements[:,2]-1])/3
y=(md.mesh.y[md.mesh.elements[:,0]-1]+md.mesh.y[md.mesh.elements[:,1]-1]+md.mesh.y[md.mesh.elements[:,2]-1])/3

triangles2=mpl.tri.Triangulation(md2.mesh.x, md2.mesh.y, md2.mesh.elements-1)
x2=(md2.mesh.x[md2.mesh.elements[:,0]-1]+md2.mesh.x[md2.mesh.elements[:,1]-1]+md2.mesh.x[md2.mesh.elements[:,2]-1])/3
y2=(md2.mesh.y[md2.mesh.elements[:,0]-1]+md2.mesh.y[md2.mesh.elements[:,1]-1]+md2.mesh.y[md2.mesh.elements[:,2]-1])/3


vmin=-600
vmax=2100
gs = GridSpec(nrows=2, ncols=16)
fig = plt.figure()
ax=fig.add_subplot(gs[0,1:15])
#triplot(triangles, lw=0.1, color='white', alpha=0.5)
tripcolor(triangles, (md.geometry.bed).flatten(),cmap='terrain', vmin=vmin, vmax=vmax)
tricontourf(md.mesh.x,md.mesh.y,np.squeeze(md.results.TransientSolution[-1].Thickness), levels=[10,10000], colors='Lightblue', alpha=.8)
plotcontour(md, md.results.TransientSolution[-1].Thickness, levels=[10], colors='Black')
plotcontour(md, md.results.TransientSolution[-1].MaskGroundediceLevelset, levels=[0], colors='Green') 

#axvspan(0,10000, color='lightgrey', alpha=.5)
yticks([10000,12500,15000,17500,20000], [0,2.5,5,7.5,10])
#xticks([0,10000,20000,30000,40000,50000,60000,70000,80000], [0,10,20,30,40,50,60,70,80])
xticks([], [])
ylabel('y-coordinates [km]', fontsize=fs)
ax21 = fig.add_subplot(gs[0,15])
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=mpl.cm.terrain, norm=norm, orientation='vertical')
ylabel('Bed elevation [m]', fontsize=fs)

ax2=fig.add_subplot(gs[1,1:15])
along(md, md.results.TransientSolution[-1].Surface, color='Blue')
along(md, md.results.TransientSolution[-1].Base, color='Blue')
ylabel('z [m]', color='Blue', fontsize=fs)
ax2.tick_params(axis='y', labelcolor='Blue')
xlabel('x-coordinates [km]', fontsize=fs)
ax3 = ax2.twinx()
along_vel(md, md.results.TransientSolution[-1].Vel, color='red')
xlim(0,85000)
xticks([0,10000,20000,30000,40000,50000,60000,70000,80000], [0,10,20,30,40,50,60,70,80])
ylabel('Velocity [m/yr]', color='Red', fontsize=fs)
ax3.tick_params(axis='y', labelcolor='Red')

ax.tick_params(axis='both', which='major', labelsize=ls)
ax21.tick_params(axis='both', which='major', labelsize=ls)
ax2.tick_params(axis='both', which='major', labelsize=ls)

plt.savefig("./Figures/Thesis/reference_glacier.svg")




acc=np.where(md2.mesh.x<10000)
triangles_acc=mpl.tri.Triangulation(md2.mesh.x[acc], md2.mesh.y[acc], (md2.mesh.elements-1)[acc])

plot(wets(55000, 0,0,0,0)[2])


ax2=fig.add_subplot(gs[1,1:15])
triplot(triangles2, lw=0.1, color='white')
tripcolor(triangles2, (md2.geometry.bed).flatten(), cmap='terrain', vmin=vmin, vmax=vmax)
