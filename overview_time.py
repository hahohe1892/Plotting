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


pattern=['*FrM1200*FlMreal180*ByH900*ByS20*.nc']

pattern2=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']

AOI=False
kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']
fs=15
ls=1

for p in pattern:
    modpath='./Models/embayments/'+p
    mod=glue_runs_md(modpath)
    all_values_long=glue_runs(modpath)
    all_values=cutallpars(all_values_long, kw, 30000)


    newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
    lon = np.array([int(i) for i in newstr.split()])
    wet_pat=[0,0,900,20000]
    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat)

cut=90
markdot=[]
colors_w=getcolors(len(mod.results.TransientSolution)-cut, 'viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=len(mod.results.TransientSolution)-cut)
colors_s=getcolors(len(mod.results.TransientSolution), 'autumn')
intervall=1

fig = plt.figure(figsize=(14,3))
gs = GridSpec(nrows=2, ncols=20)
ax0 = fig.add_subplot(gs[:,:12])
ax1 = fig.add_subplot(gs[0,13:20])
ax2 = fig.add_subplot(gs[1,13:20])
plt.sca(ax0)
for q in range(0,len(mod.results.TransientSolution)-cut,intervall):
    plotcontour(mod, mod.results.TransientSolution[q].MaskGroundediceLevelset, levels=[0], colors=colors_w[q].reshape(-1,4), linewidths=ls)
plotcontour(mod, mod.geometry.bed, levels=[0], colors='black', linewidths=ls)
ylim(10000,20000)
xlim(10000,85000)
ax0.set_yticklabels(range(0,11,2))
ax0.set_xticklabels(range(10,85,10))
ylabel('y-coordinates [km]', fontsize=fs)
xlabel('x-coordinates [km]', fontsize=fs)

plt.sca(ax1)
ran=range(0,len(all_values['GroundinglineMassFlux']))
for t in ran:
    plt.scatter(ran[t], all_values['GroundinglineMassFlux'][t], edgecolor='cornflowerblue', facecolor='none')
plt.plot(all_values['GroundinglineMassFlux'], color='cornflowerblue')
ylabel('$\it{\mathregular{Q_{GL}}}$ [km\u00b3/yr]', fontsize=fs)

plt.sca(ax2)
ran=range(0,len(all_values['dGL']))
for t in ran:
    plt.scatter(ran[t], all_values['dGL'][t], edgecolor='cornflowerblue', facecolor='none')
plt.plot(all_values['dGL'], color='cornflowerblue')
ylabel('$\it{dGL}$ [m/yr]', fontsize=fs)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()

palette2=['lightblue', 'DarkBlue']
for r,q in enumerate(pattern2):
    modpath2='./Models/embayments/'+q
    mod2=glue_runs_md(modpath2)
    all_values_long2=glue_runs(modpath2)
    all_values2=cutallpars(all_values_long2, kw, 30000)


    newstr2 = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
    lon2 = np.array([int(i) for i in newstr2.split()])
    wet_pat2=[0,0,900,20000]
    fj_chars2, rfj_chars_GL2, rfj_chars_mval2, inds_dic_GL2, inds_dic_mval2 = get_fjord(mod2, all_values2, AOI, pattern=wet_pat2)
    ran=range(0, len(all_values2['GroundinglineMassFlux']))
    for t in ran:
        plt.sca(ax1)
        plt.scatter(ran[t], all_values2['GroundinglineMassFlux'][t], edgecolor=palette2[r], facecolor='none')
        plt.sca(ax2)
        plt.scatter(ran[t], all_values2['dGL'][t], edgecolor=palette2[r], facecolor='none')
    plt.sca(ax1)
    plt.grid(True)
    plt.plot(all_values2['GroundinglineMassFlux'], color=palette2[r])
    plt.sca(ax2)
    plt.grid(True)
    plt.plot(all_values2['dGL'], color=palette2[r])

pyplot.locator_params(axis='y', nbins=6)
ax1.set_xticklabels([])
