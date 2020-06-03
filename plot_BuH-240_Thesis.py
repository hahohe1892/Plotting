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
from mpl_toolkits.axes_grid.inset_locator import inset_axes, InsetPosition

pattern=['*FrM1000*FlMreal150*BuH-240*.nc']

AOI=False

for p in pattern:
    modpath='./Models/depressions/'+p

    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)


    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)
    
markdot=[]
colors_w=getcolors(len(mod.results.TransientSolution)-20, 'viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=len(mod.results.TransientSolution)-20)
colors_s=getcolors(len(mod.results.TransientSolution), 'autumn')
intervall=1
#all_values=getallpars(mod, cut='all')[0]
#fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI=False)
plt.close('all')
begin=45000
end=65000
alongr=np.where(np.logical_and(mod.mesh.y<=15050, mod.mesh.y>=14950)) 
array=np.array((mod.mesh.x[alongr], np.squeeze(mod.geometry.bed[alongr]))) 
ind=np.argsort(array[0])  
array=array[:,ind]
fs=16
ls=12

triangles=mpl.tri.Triangulation(mod.mesh.x, mod.mesh.y, mod.mesh.elements-1)
fig= plt.figure()
camera = Camera(fig)
r=s=1
for i in range(0,len(mod.results.TransientSolution)-40,intervall):
    #plotcontour(mod, mod.geometry.bed, levels=[0], colors='black')
    #plotcontour(mod, mod.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors=colors_w[i].reshape(-1,4))
    #plotcontour(mod, mod.results.TransientSolution[i].MaskIceLevelset, levels=[0], colors=colors_s[i].reshape(-1,4))
    text(75000, 20100, 'Year {}'.format(i), fontsize=14)
    ds=drivingstress(mod, mod.results.TransientSolution[i].Surface, mod.results.TransientSolution[i].Thickness)
    bd=basal_drag(mod, mod.results.TransientSolution[i].Thickness, mod.results.TransientSolution[i].Base, mod.results.TransientSolution[i].Vel)
    tripcolor(triangles, (bd/(ds[-1]+1)).flatten(), vmin=0, vmax=2, cmap='RdYlBu_r')
    camera.snap()

title('Glacier Velocity evolution')
cbar=plt.colorbar(orientation='horizontal')
cbar.ax.set_xlabel('Velocity')
animation = camera.animate()
animation.save('celluloid_minimal.mp4')



### for depressions, calculate equi as equi=array[1]*((1023-917)/1023)*-1 where array comes from along; floating starts at 111 yrs for 360 depression

fig, ax = plt.subplots()
camera = Camera(fig)
for i in range(0,len(mod.results.TransientSolution)-0,intervall):
    for q in range(0,i, intervall):
        along(mod, mod.results.TransientSolution[q].Surface, color=colors_w[q])
        along(mod, mod.results.TransientSolution[q].Base, color=colors_w[q])
    text(75000, 1200, 'Year {}'.format(i), fontsize=14)
    greybox()
    hlines(0,0,85000, color='lightgrey')
    plot(array[0], equi, color='red')
    camera.snap()

title('Glacier profile time evolution', fontsize=20)
text(0, 70, 'Flotation Height', color='red')
text(55000,1200, 'Depression', color='darkgrey',horizontalalignment='center', fontsize=14)
cbaxes = inset_axes(ax,width="20%", height="5%", loc=1, borderpad=3)
cb1 = mpl.colorbar.ColorbarBase(cbaxes, cmap=mpl.cm.winter, norm=norm, label='Years', orientation='horizontal')
animation=camera.animate()
animation.save('test.mp4')

fig = plt.figure()
gs = GridSpec(nrows=3, ncols=3)
camera=Camera(fig)
ax0 = fig.add_subplot(gs[0,:])
ax1 = fig.add_subplot(gs[1,:])
ax2 = fig.add_subplot(gs[2,0])
ax3 = fig.add_subplot(gs[2,1])
ax4 = fig.add_subplot(gs[2, 2])
for i in range(0,len(mod.results.TransientSolution)-0,intervall):
    for q in range(0,i,intervall):
        plt.sca(ax0)
        along(mod, mod.results.TransientSolution[q].Surface, color=colors_w[q])
        along(mod, mod.results.TransientSolution[q].Base, color=colors_w[q])
        greybox()
        plt.sca(ax1)
        along(mod, mod.results.TransientSolution[q].Vel, color=colors_w[q])
        greybox()
        plt.sca(ax2)
        plt.scatter(rfj_chars_GL['dP'][q], all_values['GLvel'][q], color=colors_w[q])
        plt.sca(ax3)
        plt.scatter(rfj_chars_GL['dP'][q], all_values['GroundinglineMassFlux'][q],color=colors_w[q])
        plt.sca(ax4)
        try:
            plt.scatter(rfj_chars_GL['P'][q], all_values['dGL'][q], color=colors_w[q])
        except:
            continue
    camera.snap()

animation=camera.animate()


fig = plt.figure()
plt.tight_layout()
gs = GridSpec(nrows=3, ncols=20)
camera=Camera(fig)
ax0 = fig.add_subplot(gs[0,:19])
ax1 = fig.add_subplot(gs[1,:19])
ax2 = fig.add_subplot(gs[2,:19])


### stagnant
plt.sca(ax1)
along_evol(mod, 'viridis','','no',-1,'Surface','Base', linewidth=0.8)
hlines(0,0,85000, color='lightgrey')
axvspan(begin,end, color='gainsboro', alpha=0.5)
xlim(0,85000)
ylabel('z [m]', fontsize=fs)

axins = inset_axes(ax1, width=1.3, height=0.9)
#ax10 = plt.axes([0,0,1,1])
#ip = InsetPosition(ax1, [0.75,0.62,0.22,0.35])
#ax10.set_axes_locator(ip)

along(mod, mod.results.TransientSolution[217].Surface, color='blue', linewidth=1)
along(mod, mod.results.TransientSolution[217].Base, color='blue', linewidth=1)
along(mod, mod.geometry.bed, color='black', linewidth=1)
plot(np.linspace(0,85000,532),array[1]*((1023-917)/1023)*-1, color='red', linewidth=1)
hlines(0,0,85000, color='lightgrey')
xlim(40000,75000)
axvspan(begin,end, color='gainsboro', alpha=0.5)
ylim(-850,500)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
yticks([-500,0,500], [-500,0,500])
xticks([45000,55000,65000,75000],[45,55,65,75])

plt.sca(ax2)
for q in range(0,len(mod.results.TransientSolution)-20,intervall):
    along_vel(mod, mod.results.TransientSolution[q].Vel, color=colors_w[q], linewidth=0.8)
xlabel('x-coordinates [km]', fontsize=fs)
xlim(0,85000)
ylabel('V [m/yr]', fontsize=fs)
axvspan(begin,end, color='gainsboro', alpha=0.5)
plt.sca(ax0)
cont_all(mod, 'GL',1, 'nobar', linewidths=0.8)
plotcontour(mod, mod.geometry.bed, levels=[0], linewidths=0.8)
axvspan(begin,end, color='gainsboro', alpha=0.5)
ylim(10000,20000)
ylabel('y-coordinates [km]', fontsize=fs)
ax0.set_yticklabels(range(0,11,2))

axins1 = ax0.inset_axes([0.025, 0.3, 0.25, 0.4])
plt.axes(axins1)
plotcontour(mod, mod.results.TransientSolution[217].MaskGroundediceLevelset, levels=[0], colors='blue')
plotcontour(mod, mod.geometry.bed, levels=[0], colors='black')
axvspan(begin,end, color='gainsboro', alpha=0.5)
axins1.yaxis.tick_right()
yticks([])
xticks([0,25000,50000,75000],[0,25,50,75],fontsize=6)
#ax0.text(0.9, 0.85, 'Grounding Line', color='black', transform=ax0.transAxes,horizontalalignment='center', weight='bold')
#ax1.text(0.9, 0.85, 'Shape Profile', color='black', transform=ax1.transAxes,horizontalalignment='center', weight='bold')
#ax2.text(0.9, 0.85, 'Velocity Profile', color='black', transform=ax2.transAxes, horizontalalignment='center', weight='bold')
ax3 = fig.add_subplot(gs[:,-1])
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.viridis, norm=norm, orientation='vertical')
ylabel('Years', fontsize=fs)
ax0.set_xticklabels([])
ax1.set_xticklabels([])
ax2.set_xticklabels(range(0,85,10))
ax0.tick_params(axis='both', which='major', labelsize=ls)
ax1.tick_params(axis='both', which='major', labelsize=ls)
ax2.tick_params(axis='both', which='major', labelsize=ls)
ax3.tick_params(axis='both', which='major', labelsize=ls) 
### arrange window first
plt.savefig("./Figures/Thesis/depression_overview_new.svg")

### animation
for i in range(0,len(mod.results.TransientSolution)-0,intervall):
    for q in range(0,i,intervall):
        plt.sca(ax1)
        along(mod, mod.results.TransientSolution[q].Surface, color=colors_w[q], linewidth=0.8)
        along(mod, mod.results.TransientSolution[q].Base, color=colors_w[q], linewidth=0.8)
        hlines(0,0,85000, color='lightgrey')
        greybox()
        xlim(0,85000)
        ylabel('z [m]')
        plt.sca(ax2)
        along(mod, mod.results.TransientSolution[q].Vel, color=colors_w[q], linewidth=0.8)
        greybox()
        xlabel('x-coordinates [m]')
        xlim(0,85000)
        ylabel('V [m/a]')
    plt.sca(ax0)
    plotcontour(mod, mod.geometry.bed, levels=[0], colors='black')
    plotcontour(mod, mod.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors='Blue')
    plotcontour(mod, mod.results.TransientSolution[i].MaskIceLevelset, levels=[0], colors='red')
    ylim(10000,20000)
    ylabel('y-coordinates [m]')
    ax0.text(75000, 20100, 'Year {}'.format(i), fontsize=14)
    camera.snap()
    
ax0.text(0.02, 0.5, 'Grounding Line and Front Positions', color='black', transform=ax0.transAxes)
ax1.text(0.02, 0.5, 'Shape Profile', color='black', transform=ax1.transAxes)
ax2.text(0.02, 0.5, 'Velocity Profile', color='black', transform=ax2.transAxes)
ax3 = fig.add_subplot(gs[1,-1])
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.winter, norm=norm, label='Years', orientation='vertical')
animation=camera.animate()
