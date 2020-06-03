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

bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*off.nc']#, '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*']
embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off.nc','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']#, '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
megapat=bottlenecks+embayments+depressions+bumps
paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']

loadpath='./Models/'

fs=16
ls=12

fig, ax = plt.subplots(3,1)
ax[0].tick_params(axis='both', which='major', labelsize=ls)
ax[1].tick_params(axis='both', which='major', labelsize=ls)
ax[2].tick_params(axis='both', which='major', labelsize=ls)
for p in megapat:

    z=np.nonzero(megapat==np.intersect1d(megapat, p))[0]

    for r in paths:
        try:
            mods=glob.glob(loadpath+r+p)
            mods.sort(key=os.path.getmtime)
            md=loadmodel(mods[0])

        except:
            continue
        

    if z==0:
        palette='Gold'
        marker='-'
    if z==1:
        palette='Orange'
    if z==2:
        palette='Saddlebrown'
    #if z==3:
    #    palette='LightPink'
    #    marker='--'
    #if z==4:
    #    palette='DeepPink'
    #    marker='--'
    if z==3:
        palette='SkyBlue'
        marker='-'
    if z==4:
        palette='MediumBlue'
    if z==5:
        palette='DarkBlue'
    #if z==8:
    #    palette='SkyBlue'
    #    marker='--'
    #if z==9:
    #    palette='DarkBlue'
    #    marker='--'
    if z==6:
        palette='LightGreen'
        marker='-'
    if z==7:
        palette='SeaGreen'
    if z==8:
        palette='DarkGreen'
    if z==9:
        palette='Silver'
    if z==10:
        palette='Gray'
    if z==11:
        palette='Black'
    #if z==12:
    #    palette='Gold'
    #if z==17:
    #    palette='Orange'
    #if z==18:
    #    palette='SaddleBrown'

    plt.sca(ax[0])
    plotcontour(md, md.mask.groundedice_levelset, levels=[0], colors=palette, linestyles=marker)#, linewidths=0.8)
    plotcontour(md, md.geometry.bed, levels=[0], colors=palette, linestyles=marker)#, linewidths=0.8)
    plt.sca(ax[1])
    along(md, md.geometry.surface, color=palette, linestyle=marker)#, linewidth=0.8)
    along(md, md.geometry.base, color=palette, linestyle=marker)#, linewidth=0.8)
    plt.sca(ax[2])
    along_vel(md, np.sqrt((md.initialization.vx)**2+(md.initialization.vy)**2), color=palette, linestyle=marker, x_area=70000)#, linewidth=0.8)


plt.sca(ax[0])
xticks(range(0,85000,10000), [])
ylabel('y-coordinates [km]', fontsize=fs)
yticks(range(10000,22000,2000), range(0,24,2))
xlim(20000,85000)
ylim(10000,20000)
greybox()

plt.sca(ax[1])
hlines(0,0,85000, color='lightgrey')
xticks(range(0,85000,10000), [])
ylabel('z [m]', fontsize=fs)
xlim(20000,85000)
greybox()

plt.sca(ax[2])
xlabel('x-coordinates [km]', fontsize=fs)
ylabel('V [m/yr]', fontsize=fs)
xlim(20000,85000)
xticks(range(20000,85000,10000), range(20,85,10))
greybox()

plt.savefig("./Figures/Thesis/ss_glaciers_overview.svg")
