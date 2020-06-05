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

AOI=False

mods=['GeomProj_SpinUp_load_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT50_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT725_inH0_smb55_smbPos30000_funnel300_FullMelt.nc', 'GeomProj_SpinUp_load_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT200_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT825_inH0_smb55_smbPos30000_funnel300_FullMelt.nc', 'GeomProj_extenddomain_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT200_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim65000_accT1025_inH0_smb55_smbPos30000_funnel300_FullMelt.nc','GeomProj_extenddomain_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT100_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim52000_accT1125_inH0_smb55_smbPos30000_funnel300_FullMelt.nc','GeomProj_extenddomain_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT100_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim47000_accT1225_inH0_smb55_smbPos30000_funnel300_FullMelt.nc','GeomProj_extenddomain_SG_spcvx50_NL-450_FrM800_FlMreal120_FC40_FT100_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim43000_accT1325_inH0_smb55_smbPos30000_funnel300_FullMelt.nc']

modpath=[]
for p in mods:
    modpath.append('./Models/linear/'+p)

mod=glue_runs_md('dummy', manual_mods=modpath)
all_values=glue_runs('dummy', manual_mods=modpath)

cut=0
markdot=[]
colors_w=getcolors(len(mod.results.TransientSolution)-cut, 'viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=len(mod.results.TransientSolution)-cut)
colors_s=getcolors(len(mod.results.TransientSolution), 'autumn')
intervall=1
fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI=False)

### animation
fs=16
ls=12

fig = plt.figure()#dpi=400)
gs = GridSpec(nrows=3, ncols=20)
camera=Camera(fig)
ax0 = fig.add_subplot(gs[0,:19])
ax1 = fig.add_subplot(gs[1,:19])
ax2 = fig.add_subplot(gs[2,:19])
ax0.set_xticklabels([])
ax1.set_xticklabels([])
ax2.set_xticklabels(range(0,85,10))
ax3 = fig.add_subplot(gs[:,-1])
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.viridis, norm=norm, orientation='vertical')
ylabel('Years', fontsize=fs)
ax0.tick_params(axis='both', which='major', labelsize=ls)
ax1.tick_params(axis='both', which='major', labelsize=ls)
ax2.tick_params(axis='both', which='major', labelsize=ls)
ax3.tick_params(axis='both', which='major', labelsize=ls)

for q in range(0,len(mod.results.TransientSolution)-cut,intervall):
    plt.sca(ax1)
    along(mod, mod.results.TransientSolution[q].Surface, color=colors_w[q], linewidth=0.8)
    along(mod, mod.results.TransientSolution[q].Base, color=colors_w[q], linewidth=0.8)
    hlines(0,0,85000, color='lightgrey')
    xlim(0,85000)
    ylabel('z [m]', fontsize=fs)
    plt.sca(ax0)
    plotcontour(mod, mod.results.TransientSolution[q].MaskGroundediceLevelset, levels=[0], colors=colors_w[q].reshape(-1,4))
    ylim(10000,20000)
    ylabel('y-coordinates [km]', fontsize=fs)
    ax0.set_yticklabels(range(0,11,2))
    plt.sca(ax2)
    along_vel(mod, mod.results.TransientSolution[q].Vel, color=colors_w[q], linewidth=0.8)
    xlabel('x-coordinates [km]', fontsize=fs)
    xlim(0,85000)
    ylabel('V [m/yr]', fontsize=fs)
plt.sca(ax0)
plotcontour(mod, mod.geometry.bed, levels=[0], colors='black', linewidths=2)

fig = plt.figure(2)
plt.tight_layout()
gs = GridSpec(nrows=6, ncols=2)
ax0 = fig.add_subplot(gs[0:2,0])
ax1 = fig.add_subplot(gs[2:4,0])
ax2 = fig.add_subplot(gs[4:6,0])
ax3 = fig.add_subplot(gs[0:2,1])
ax4 = fig.add_subplot(gs[2:4,1])
ax5 = fig.add_subplot(gs[4:6,1])
ax0.tick_params(axis='both', which='major', labelsize=ls)
ax1.tick_params(axis='both', which='major', labelsize=ls)
ax2.tick_params(axis='both', which='major', labelsize=ls)
ax3.tick_params(axis='both', which='major', labelsize=ls) 
ax4.tick_params(axis='both', which='major', labelsize=ls) 
ax5.tick_params(axis='both', which='major', labelsize=ls)

palette='viridis'

plt.sca(ax0)
val_evol(palette, '','no', all_values['GroundinglineMassFlux'])
ylabel('$\it{\mathregular{Q_{GL}}}$ \n [km\u00b3/yr]', fontsize=fs)
plt.plot([], label='a)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

plt.sca(ax1)
val_evol(palette,'', 'no',all_values['dGL'])
ylabel('$\it{dGL}$ [m/yr]', fontsize=fs)
plt.plot([], label='b)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

plt.sca(ax2)
val_evol(palette, '','no',all_values['GLvel'])
ylabel('$\it{\mathregular{V_{GL}}}$ [m/yr]', fontsize=fs)
xlabel('Years', fontsize=fs)
plt.plot([], label='c)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

plt.sca(ax3)
val_evol(palette,'','no', ((np.array(all_values['TotalCalvingFluxLevelset'])/917)*31536000)/1e9)
ylabel('$\it{C}$ \n [km\u00b3/yr]', fontsize=fs)
plt.plot([], label='d)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

plt.sca(ax4)
val_evol(palette, '','no',np.array(all_values['IceVolume'])/1e9)
ylabel('$\it{I}$ [km\u00b3]', fontsize=fs)
plt.plot([], label='e)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

plt.sca(ax5)
val_evol(palette, '','no',np.array(all_values['GLval'])/1000)
ylabel('$\it{\mathregular{x_{GL}}}$  [km]', fontsize=fs)
xlabel('Years', fontsize=fs)
plt.plot([], label='f)')
plt.legend(handlelength=0, frameon=False, fontsize=fs)

ax0.set_xticklabels([])
ax1.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])

ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.yaxis.tick_right()
ax5.yaxis.set_label_position("right")
ax5.yaxis.tick_right()    

for a in gcf().axes:
    a.yaxis.set_major_locator(plt.MaxNLocator(4))
subplots_adjust(hspace=0.5)


fs=15
l=12
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
ylim(10000,19000)
plotcontour(mod, mod.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors='green', linewidths=2)
plotcontour(mod, mod.results.TransientSolution[i].Thickness, levels=[10], colors='black', linewidths=2)
yticks([12000,15000,18000], [2,5,8], fontsize=l)
xticks([35000,45000,55000,65000,75000],[35,45,55,65,75], fontsize=l)
xlabel('along-flow Distance [km]', fontsize=fs)
ylabel('across-flow Distance [km]', fontsize=fs)
ax21 = fig.add_subplot(gs[:,10])
norm = colors.Normalize(vmin=vmin/1000, vmax=vmax/1000)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=mpl.cm.RdYlBu_r, norm=norm, orientation='vertical')
xticks(fontsize=l)
ylabel('Longitudinal Stress [kPa]', fontsize=fs)
yticks(fontsize=l)
