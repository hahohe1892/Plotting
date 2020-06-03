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

embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']

bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*']

bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']


megapat=[embayments, depressions,bottlenecks, bumps]

paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']

kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']

AOI=True

markdot=[]

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
        
fig = plt.figure(1)
xlabel('dS [m\u00b2/100m]')
ylabel('$\mathregular{Q_{GL}}$ [km\u00b3/yr]')
fig2 = plt.figure(2)
xlabel('dS [m\u00b2/100m]')
ylabel('$\mathregular{V_{GL}}$ [m/yr]')

fig3, ax = plt.subplots()
fig3.subplots_adjust(right=0.75)
fig3.subplots_adjust(left=0.25)
xlabel('S [km\u00b2]')
#ax.yaxis.label.set_color('blue')
ax.tick_params(axis='y', colors='blue')
par1 = ax.twinx()
#par1.yaxis.label.set_color('green')
par1.tick_params(axis='y', colors='green')
par2 = ax.twinx()
#par2.yaxis.label.set_color('orange')
par2.tick_params(axis='y', colors='orange')
ylabel('dGL [m/yr]')
par3 = ax.twinx()
#par3.yaxis.label.set_color('gray')
par3.tick_params(axis='y', colors='gray')
ylabel('dGL [m/yr]')

par2.spines["right"].set_position(("axes", 1.2))
par3.spines["left"].set_position(("axes", -0.2))

make_patch_spines_invisible(par2)
make_patch_spines_invisible(par3)
par2.spines["right"].set_visible(True)
par3.spines["left"].set_visible(True)
par3.yaxis.set_label_position('left')
par3.yaxis.set_ticks_position('left')

for pattern in megapat:
    n=int(np.nonzero([x==pattern for x in megapat])[0])

    geompath='./Models/'+paths[n]
    
    for p in pattern:

        modpath=geompath+p
        mod=glue_runs_md(modpath)
        all_values_long=glue_runs(modpath)
        all_values=cutallpars(all_values_long, kw, 25000)

        
        newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
        lon = np.array([int(i) for i in newstr.split()])
        wet_pat=lon[[2,4,5,7]]
        fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat)
        
        r=np.nonzero(pattern==np.intersect1d(pattern, p))[0]
        z=n*3+r

        ### plotting ###
        
        if z==0:
            palette='SkyBlue'
        if z==1:
            palette='MediumBlue'
        if z==2:
            palette='DarkBlue'
        if z==3:
            palette='LightGreen'
        if z==4:
            palette='SeaGreen'
        if z==5:
            palette='DarkGreen'
        if z==6:
            palette='Gold'
        if z==7:
            palette='Orange'
        if z==8:
            palette='Saddlebrown'
        if z==9:
            palette='Silver'
        if z==10:
            palette='Gray'
        if z==11:
            palette='Black'

        plt.figure(1)
        inds_dic_GL['dWA']=np.array(inds_dic_GL['dWA'])[np.array(inds_dic_GL['dWA'])<len(all_values['GroundinglineMassFlux'])]
        par=np.array(all_values['GroundinglineMassFlux'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'],par,color=palette)
        grid(False)

        plt.figure(2)
        par=np.array(all_values['GLvel'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'], par,color=palette)
        grid(False)

        plt.figure(3)
        if n==0:
            plt.sca(ax)
        if n==1:
            plt.sca(par1)
        if n==2:
            plt.sca(par2)
        if n==3:
            plt.sca(par3)

        inds_dic_GL['WA']=np.array(inds_dic_GL['WA'])[np.array(inds_dic_GL['WA'])<len(all_values['dGL'])]
        par=np.array(all_values['dGL'])[inds_dic_GL['WA']]
        plt.scatter(np.array(rfj_chars_GL['WA'])/1e6,par,color=palette)
        grid(False)
        


        
