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
import pickle
from scipy.stats import zscore
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

a=open('OneFileCont.pkl', 'rb')
all_vals=pickle.load(a)
a.close()

interest=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off*','*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*', '*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*']

begin=45000
end=65000

cuts=[90, 124, 0, 0]

fig = plt.figure()
gs = GridSpec(nrows=4, ncols=20)
ax=[]
bar=[]
ax.append(fig.add_subplot(gs[0,:19]))
ax.append(fig.add_subplot(gs[1,:19]))
ax.append(fig.add_subplot(gs[2,:19]))
ax.append(fig.add_subplot(gs[3,:19]))
bar.append(fig.add_subplot(gs[0,19:20]))
bar.append(fig.add_subplot(gs[1,19:20]))
bar.append(fig.add_subplot(gs[2,19:20]))
bar.append(fig.add_subplot(gs[3,19:20]))

for i,element in enumerate(interest):
    if i < 2:
        cont=all_vals[element]

        colors_w=getcolors(len(cont)-cuts[i], 'viridis')
        norm = mpl.colors.Normalize(vmin=0, vmax=len(cont)-cuts[i])

        plt.sca(ax[i])
        for q in range(len(cont)-cuts[i]):
            plot(cont[q][0], cont[q][1], color=colors_w[q])
        plot(cont[len(cont)-1][0], cont[len(cont)-1][1], color='black')
        plot(cont[len(cont)-2][0], cont[len(cont)-2][1], color='black')
        ylim(10000,20000)
        yticks([10000,15000,20000],[0,5,10])
        ylabel('y [km]')

    if i == 2:
        m='./Models/depressions/'+element+'.nc'
    if i == 3:
        m='./Models/bumps/'+element+'.nc'
    if i ==2 or i==3:
        mod=glue_runs_md(m)
        plt.sca(ax[i])
        along_evol(mod, 'viridis','','no',-1,'Surface','Base', linewidth=0.8)
        hlines(0,0,85000, color='lightgrey')
        along(mod, mod.geometry.bed, color='black')
        ylabel('z [m]')
        yticks([-500,0,1000,2000],[-500,0,1000,2000])

    xlim(20000,85000)
    axvspan(begin,end, color='gainsboro', alpha=0.5)


    if i == 2:
        alongr=np.where(np.logical_and(mod.mesh.y<=15050, mod.mesh.y>=14950))
        array=np.array((mod.mesh.x[alongr], np.squeeze(mod.geometry.bed[alongr])))
        ind=np.argsort(array[0])
        array=array[:,ind]
        axins = inset_axes(ax[i], width=1.3, height=0.7)
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
        

    if i<3:
        ax[i].set_xticklabels([])
    else:
        xlabel('x [km]')
        xticks([20000,30000,40000,50000,60000,70000,80000],[20,30,40,50,60,70,80])

    b1 = mpl.colorbar.ColorbarBase(bar[i], cmap=mpl.cm.viridis, norm=norm, orientation='vertical')
    bar[i].set_yticks([0,100,200],[0,100,200])


