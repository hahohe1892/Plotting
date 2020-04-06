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

#pattern=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH-900*','*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*','*FrM1200*FlMreal180*ByH1350*','*FrM1200*FlMreal180*BuH120*','*FrM1200*FlMreal180*BuH180*','*FrM1200*FlMreal180*BuH240*','*FrM1200*FlMreal180*BuH-120*','*FrM1000*FlMreal150*BuH-240*','*FrM1200*FlMreal180*BuH-360*']
#pattern=['*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*off.nc','*FrM1200*FlMreal180*ByH1350*']

bottlenecks=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH-900*off.nc', '*FrM1200*ByH-900*asy*', '*FrM1200*ByH-1800*asy*']
embayments=['*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*off.nc','*FrM1200*FlMreal180*ByH1350*', '*FrM1200*ByH900*asy*', '*FrM1200*ByH1800*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*','*FrM1000*FlMreal150*BuH-240*','*FrM1200*FlMreal180*BuH-360*']
bumps=['*FrM1200*FlMreal180*BuH120*','*FrM1200*FlMreal180*BuH180*','*FrM1200*FlMreal180*BuH240*']

name_dict=['bottlenecks','embayments','depressions','bumps']

megapat=[bottlenecks, embayments, depressions, bumps]


AOI=True

for pattern in megapat:

    n=int(np.nonzero([x==pattern for x in megapat])[0])
    fig = plt.figure(n, figsize=(30,10))
    gs = GridSpec(nrows=2, ncols=3)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2= fig.add_subplot(gs[0,2])

    s=45
    fs=16
    ls=14

    for p in pattern:
        modpath='./Models/'+p
        mod=glue_runs_md(modpath)
        all_values=glue_runs(modpath)

        crap= plt.figure(10)
        #all_values=getallpars(mod, cut=(0,0))[0]

        fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)
        plt.close(10)
        markdot=[]


        r=z=np.nonzero(pattern==np.intersect1d(pattern, p))[0]
        z=n*3+z
        if z==0:
            palette='Gold'
        if z==1:
            palette='Orange'
        if z==2:
            palette='Saddlebrown'
        if z==3:
            palette='SkyBlue'
        if z==4:
            palette='DodgerBlue'
        if z==5:
            palette='DarkBlue'
        if z==6:
            palette='LightGreen'
        if z==7:
            palette='MediumSeaGreen'
        if z==8:
            palette='DarkGreen'
        if z==9:
            palette='Silver'
        if z==10:
            palette='Gray'
        if z==11:
            palette='Black'       

        if 'asy' in p and r==3:
            palette='LightPink'
            s=20
        if 'asy' in p and r==4:
            palette='Violet'
            s=20

            
        plt.sca(ax0)
        inds_dic_GL['dP']=np.array(inds_dic_GL['dP'])[np.array(inds_dic_GL['dP'])<len(all_values['GroundinglineMassFlux'])]
        par=np.array(all_values['GroundinglineMassFlux'])[inds_dic_GL['dP']]
        plt.scatter(rfj_chars_GL['dP'],par,color=palette, s=s)
        ylabel('GL Mass Flux \n [km\u00b3/yr]', fontsize=fs)
        xlabel('dP [m\u00b2]', fontsize=fs)
        plt.sca(ax1)
        par=np.array( all_values['GLvel'])[inds_dic_GL['dP']]
        plt.scatter(rfj_chars_GL['dP'], par,color=palette, s=s)
        ylabel('$\mathregular{V_{GL}}$ [m/yr]', fontsize=fs)
        xlabel('dP [m\u00b2]', fontsize=fs)
        plt.sca(ax2)
        inds_dic_GL['P']=np.array(inds_dic_GL['P'])[np.array(inds_dic_GL['P'])<len(all_values['dGL'])]
        par=np.array(all_values['dGL'])[inds_dic_GL['P']]
        try:
            plt.scatter(np.array(rfj_chars_GL['P'])/1e6,par,color=palette, s=s)
        except:
            plt.scatter(np.array(rfj_chars_GL['P'][:-1])/1e6,par,color=palette, s=s)        
        ylabel('dGL [m/yr]', fontsize=fs)
        xlabel('P [km\u00b2]', fontsize=fs)
        
    subplots_adjust(wspace=0.3)
    for a in gcf().axes:
        a.tick_params(axis='x', labelsize=ls)
        a.tick_params(axis='y', labelsize=ls)
    #plt.savefig(str(name_dict[n])+'.png')

        
