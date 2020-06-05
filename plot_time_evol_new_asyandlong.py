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

embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH901*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP50000*ByS30000*']
depressions=['*FrM1200*FlMreal180*BuH-240*BuP50000*BuS30000*ByH0*ByP0*ByS0*']
bottlenecks=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP50000*ByS30000*']

bumps=['*FrM1200*FlMreal180*BuH180*BuP50000*BuS30000*ByH0*ByP0*ByS0*']

megapat=[embayments, depressions,bottlenecks, bumps]

paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']

kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']

AOI=False

markdot=[]

fs=15
ls=12
#rc('text', usetex=True)
for pattern in megapat:
    
    n=int(np.nonzero([x==pattern for x in megapat])[0])
    
    fig = plt.figure(n)
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

    
    geompath='./Models/'+paths[n]
    
    for p in pattern:

        modpath=geompath+p
        mod=glue_runs_md(modpath)
        all_values_long=glue_runs(modpath)
        all_values=cutallpars(all_values_long, kw, 30000)

        r=np.nonzero(pattern==np.intersect1d(pattern, p))[0]
        
        newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
        lon = np.array([int(i) for i in newstr.split()])
        wet_pat_all=lon[[2,3,4,5,6,7]]
        wet_pat=lon[[2,4,5,7]]
        if n == 0 or n == 2:
            if r==0 or r == 2:
                asy='yes'
            else:
                asy='no'
        else:
            asy='no'
        if n == 0 or n == 2:
            AOI=[wet_pat_all[4]+(wet_pat_all[5]/2), wet_pat_all[4]-(wet_pat_all[5]/2)]
        if n == 1 or n ==3:
            AOI=[wet_pat_all[1]+(wet_pat_all[2]/2), wet_pat_all[1]-(wet_pat_all[2]/2)]
            
        fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat_all, asy=asy)
        
        z=n*10+r

        ### plotting ###
        
        if z==0:
            palette='SkyBlue'
        if z==1:
            palette='MediumBlue'
        if z==2:
            palette='DarkBlue'
        if z ==10:
            palette='SeaGreen'
        if z==20:
            palette='Gold'
        if z==21:
            palette='Orange'
        if z==22:
            palette='Saddlebrown'
        if z==30:
            palette='Gray'


        plt.sca(ax0)
        val_evol(palette, '','no', all_values['GroundinglineMassFlux'], marker='+')
        ylabel('$\it{\mathregular{Q_{GL}}}$ [km\u00b3/yr]', fontsize=fs)
        if r==2:
            plt.plot([], label='a)')
            plt.legend(handlelength=0, frameon=False, fontsize=fs)
            
        plt.sca(ax1)
        val_evol(palette,'', 'no',all_values['dGL'], marker='+')
        ylabel('$\it{dGL}$ [m/yr]', fontsize=fs)
        if r==2:
            plt.plot([], label='b)')
            plt.legend(handlelength=0, fontsize=fs, frameon=False)
            
        plt.sca(ax2)
        val_evol(palette, '','no',all_values['GLvel'], marker='+')
        ylabel('$\it{\mathregular{V_{GL}}}$ [m/yr]', fontsize=fs)
        xlabel('Years', fontsize=fs)
        if r==2:
            plt.plot([], label='c)')
            plt.legend(handlelength=0, fontsize=fs, frameon=False)
            
        plt.sca(ax3)
        val_evol(palette,'','no', ((np.array(all_values['TotalCalvingFluxLevelset'])/917)*31536000)/1e9, marker='+')
        ylabel('$\it{C}$ [km\u00b3/yr]', fontsize=fs)
        if r==2:
            plt.plot([], label='d)')
            plt.legend(handlelength=0, fontsize=fs, frameon=False)
            
        plt.sca(ax4)
        val_evol(palette, '','no',np.array(all_values['IceVolume'])/1e9, marker='+')
        ylabel('$\it{I}$ [km\u00b3]', fontsize=fs)
        if r==2:
            plt.plot([], label='e)')
            plt.legend(handlelength=0, fontsize=fs, frameon=False)
            
        plt.sca(ax5)
        val_evol(palette, '','no',np.array(all_values['GLval'])/1000, marker='+')
        ylabel('$\it{\mathregular{x_{GL}}}$  [km]', fontsize=fs)
        xlabel('Years', fontsize=fs)
        if r==2:
            plt.plot([], label='f)')
            plt.legend(handlelength=0, fontsize=fs, frameon=False)
            
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
    
    ax5.axhspan(45,65, color='lightgrey', alpha=0.5, zorder=-1)

### arrange window first
plt.savefig("./Figures/Thesis/embayments_time.svg")
    
