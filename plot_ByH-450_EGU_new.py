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
pattern=['*FrM800*FlMreal120*ByH-450*.nc']

AOI=False

for p in pattern:
    modpath='./Models/bottlenecks/'+p

    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)

markdot=[]
cut=124
colors_w=getcolors(len(mod.results.TransientSolution)-cut, 'viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=len(mod.results.TransientSolution)-cut)
colors_s=getcolors(len(mod.results.TransientSolution), 'autumn')
intervall=1
fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI=False)

### animation

fig = plt.figure()#dpi=400)
gs = GridSpec(nrows=2, ncols=20)
camera=Camera(fig)
ax0 = fig.add_subplot(gs[0,:19])
ax1 = fig.add_subplot(gs[1,:19])
ax0.set_xticklabels([])
ax1.set_xticklabels(range(0,85,10))
ax3 = fig.add_subplot(gs[:,-1])
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.viridis, norm=norm, label='Years', orientation='vertical')

for q in range(0,len(mod.results.TransientSolution)-cut,intervall):
    plt.sca(ax1)
    along(mod, mod.results.TransientSolution[q].Surface, color=colors_w[q], linewidth=0.8)
    along(mod, mod.results.TransientSolution[q].Base, color=colors_w[q], linewidth=0.8)
    along(mod, mod.geometry.bed, color='black')
    hlines(0,0,85000, color='lightgrey')
    greybox()
    xlim(0,85000)
    ylabel('z [m]')
    xlabel('x [km]')
    plt.sca(ax0)
    plotcontour(mod, mod.results.TransientSolution[q].MaskGroundediceLevelset, levels=[0], colors=colors_w[q].reshape(-1,4))
    plotcontour(mod, mod.geometry.bed, levels=[0], colors='black')
    ylim(10000,20000)
    ylabel('y [km]')
    ax0.set_yticklabels(range(0,11,2))
    ax0.text(75000, 20100, 'Year {}'.format(q), fontsize=14)
    greybox()
    plt.sca(ax3)
    cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.viridis, norm=norm, label='Years', orientation='vertical')
    camera.snap()
    
plt.sca(ax1)
along_evol(mod, 'viridis','','no',-cut,'Surface','Base', linewidth=0.8)
along(mod, mod.geometry.bed, color='black')
hlines(0,0,85000, color='lightgrey')
greybox()
xlim(0,85000)
ylabel('z [m]')
xlabel('x [km]')
plt.sca(ax0)
cont_all(mod, 'GL',1, 'nobar', linewidths=0.8, cut=cut)
plotcontour(mod, mod.geometry.bed, levels=[0], colors='black')
ylim(10000,20000)
ylabel('y [km]')
ax0.set_yticklabels(range(0,11,2))
greybox()
along(mod, mod.geometry.bed, color='black')
plt.sca(ax3)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=mpl.cm.viridis, norm=norm, label='Years', orientation='vertical')
camera.snap()
animation=camera.animate()
