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


pattern=['*FrM1200*FlMreal180*ByH900*ByS20*.nc','*FrM800*FlMreal120*BuH-360*BuS20*.nc','*FrM800*FlMreal120*ByH-450*ByS20*.nc','*FrM1200*FlMreal180*BuH120*BuS20*.nc']
paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']
cuts=[90,90,90,90]
wet_pats=[[0,0,900,20000], [-360,20000,0,0],[0,0,-450,20000],[120,20000,0,0]]
colors=['cornflowerblue', 'DarkGreen', 'Gold', 'lightgray']

pattern2=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*', '*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
palette2=['lightblue','DarkBlue','palegreen','mediumseagreen','Orange','Saddlebrown','Gray','Black']
AOI=False
kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']
fs=15
ls=1

fig=plt.figure(figsize=(14,15))
gs = GridSpec(nrows=8, ncols=20)
ax=[]
for t in range(0,8,2):
    ax.append(fig.add_subplot(gs[t:t+2,:12]))
for z in range(8):
    ax.append(fig.add_subplot(gs[z,13:20]))

for v,p in enumerate(pattern):
    modpath='./Models/'+paths[v]+p
    mod=glue_runs_md(modpath)

    colors_w=getcolors(len(mod.results.TransientSolution)-cuts[v], 'viridis')
    norm = mpl.colors.Normalize(vmin=0, vmax=len(mod.results.TransientSolution)-cuts[v])
    colors_s=getcolors(len(mod.results.TransientSolution), 'autumn')
    intervall=1

    plt.sca(ax[v])
    for q in range(0,len(mod.results.TransientSolution)-cuts[v],intervall):
        plotcontour(mod, mod.results.TransientSolution[q].MaskGroundediceLevelset, levels=[0], colors=colors_w[q].reshape(-1,4), linewidths=ls)
    plotcontour(mod, mod.geometry.bed, levels=[0], colors='black', linewidths=ls)

    ylim(10000,20000)
    xlim(10000,85000)
    ylabel('y-coordinates [km]', fontsize=fs)
    ax[v].set_yticklabels(range(0,11,2))
    if v==3:
        ax[v].set_xticklabels(range(10,85,10))
        xlabel('x-coordinates [km]', fontsize=fs)
    else:
        ax[v].set_xticklabels([])
        
    all_values_long=glue_runs(modpath)
    all_values=cutallpars(all_values_long, kw, 30000)

    wet_pat=wet_pats[v]
    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat)

    ran=range(0,len(all_values['GroundinglineMassFlux']))
    for t in ran:
        plt.sca(ax[4+v*2])
        plt.scatter(ran[t], all_values['GroundinglineMassFlux'][t], edgecolor=colors[v], facecolor='none')
        plt.sca(ax[4+v*2+1])
        plt.scatter(ran[t], all_values['dGL'][t], edgecolor=colors[v], facecolor='none')
    plt.sca(ax[4+v*2])
    plt.plot(all_values['GroundinglineMassFlux'], color=colors[v])
    plt.sca(ax[4+v*2+1])
    plt.plot(all_values['dGL'], color=colors[v])



    
for r,q in enumerate(pattern2):
    modpath2='./Models/'+paths[int(floor(r/2))]+q
    mod2=glue_runs_md(modpath2)
    all_values_long2=glue_runs(modpath2)
    all_values2=cutallpars(all_values_long2, kw, 30000)


    newstr2 = ''.join((ch if ch in '0123456789-' else ' ') for ch in q)
    lon2 = np.array([int(i) for i in newstr2.split()])
    wet_pat2=lon2[[2,4,5,7]]
    fj_chars2, rfj_chars_GL2, rfj_chars_mval2, inds_dic_GL2, inds_dic_mval2 = get_fjord(mod2, all_values2, AOI, pattern=wet_pat2)
    ran=range(0, len(all_values2['GroundinglineMassFlux']))
    for t in ran:
        plt.sca(ax[4+int(floor(r/2))*2])
        plt.scatter(ran[t], all_values2['GroundinglineMassFlux'][t], edgecolor=palette2[r], facecolor='none')
        plt.sca(ax[4+int(floor(r/2))*2+1])
        plt.scatter(ran[t], all_values2['dGL'][t], edgecolor=palette2[r], facecolor='none')
    plt.sca(ax[4+int(floor(r/2))*2])
    plt.grid(True)
    plt.plot(all_values2['GroundinglineMassFlux'], color=palette2[r])
    plt.sca(ax[4+int(floor(r/2))*2+1])
    plt.grid(True)
    plt.plot(all_values2['dGL'], color=palette2[r])

pyplot.locator_params(axis='y', nbins=6)
ax1.set_xticklabels([])
