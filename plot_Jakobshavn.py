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

md=loadmodel('Models/embayments/GeomProj_extenddomain_SG_spcvx50_NL-450_FrM1200_FlMreal180_FC40_FT300_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH1800_ByP55000_ByS20000_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT1275_inH0_smb55_smbPos30000_funnel300_FullMelt_asymmetric.nc')
md=loadmodel('Models/embayments/GeomProj_extenddomain_SG_spcvx50_NL-450_FrM1200_FlMreal180_FC40_FT300_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH901_ByP55000_ByS20000_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT1225_inH0_smb55_smbPos30000_funnel300_FullMelt_asymmetric.nc')

triangles=mpl.tri.Triangulation(md.mesh.x, md.mesh.y, md.mesh.elements-1)
x=(md.mesh.x[md.mesh.elements[:,0]-1]+md.mesh.x[md.mesh.elements[:,1]-1]+md.mesh.x[md.mesh.elements[:,2]-1])/3
y=(md.mesh.y[md.mesh.elements[:,0]-1]+md.mesh.y[md.mesh.elements[:,1]-1]+md.mesh.y[md.mesh.elements[:,2]-1])/3


colors=getcolors(51, 'hsv')  

vmin=-600
vmax=2100
gs = GridSpec(nrows=1, ncols=16)
fig = plt.figure()
ax=fig.add_subplot(gs[0,1:15])
tripcolor(triangles, (md.geometry.bed).flatten(),cmap='terrain', vmin=vmin, vmax=vmax)
tricontourf(md.mesh.x,md.mesh.y,np.squeeze(md.results.TransientSolution[230].Thickness), levels=[10,10000], colors='Lightblue', alpha=.8)
for i in range(230,280):
    plotcontour(md, md.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors=colors[i-230].reshape(-1,4)) 

#axvspan(0,10000, color='lightgrey', alpha=.5)
yticks([10000,12500,15000,17500,20000], [0,2.5,5,7.5,10])
#xticks([0,10000,20000,30000,40000,50000,60000,70000,80000], [0,10,20,30,40,50,60,70,80])
xticks([], [])
ylabel('y-coordinates [km]')
ax21 = fig.add_subplot(gs[0,15])
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=mpl.cm.hsv, norm=norm, label='Elevation [m]', orientation='vertical')








ax2=fig.add_subplot(gs[1,1:15])
along(md, md.results.TransientSolution[-1].Surface, color='Blue')
along(md, md.results.TransientSolution[-1].Base, color='Blue')
ylabel('Elevation [m]', color='Blue')
ax2.tick_params(axis='y', labelcolor='Blue')
xlabel('x-coordinates [km]')
ax3 = ax2.twinx()
along_vel(md, md.results.TransientSolution[-1].Vel, color='red')
xlim(0,85000)
xticks([0,10000,20000,30000,40000,50000,60000,70000,80000], [0,10,20,30,40,50,60,70,80])
ylabel('Velocity [m/yr]', color='Red')
ax3.tick_params(axis='y', labelcolor='Red')



acc=np.where(md2.mesh.x<10000)
triangles_acc=mpl.tri.Triangulation(md2.mesh.x[acc], md2.mesh.y[acc], (md2.mesh.elements-1)[acc])

plot(wets(55000, 0,0,0,0)[2])


ax2=fig.add_subplot(gs[1,1:15])
triplot(triangles2, lw=0.1, color='white')
tripcolor(triangles2, (md2.geometry.bed).flatten(), cmap='terrain', vmin=vmin, vmax=vmax)
