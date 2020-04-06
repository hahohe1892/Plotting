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
pattern=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH-900*off.nc']



bottlenecks=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH-900*off.nc']#, '*FrM1200*ByH-900*asy*', '*FrM1200*ByH-1800*asy*']
embayments=['*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*off.nc','*FrM1200*FlMreal180*ByH1350*']#, '*FrM1200*ByH900*asy*', '*FrM1200*ByH1800*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*','*FrM1000*FlMreal150*BuH-240*','*FrM1200*FlMreal180*BuH-360*']
bumps=['*FrM1200*FlMreal180*BuH120*','*FrM1200*FlMreal180*BuH180*','*FrM1200*FlMreal180*BuH240*']

name_dict=['bottlenecks','embayments','depressions','bumps']

megapat=[bottlenecks, embayments, depressions, bumps]

AOI=False

for pattern in megapat:

    n=int(np.nonzero([x==pattern for x in megapat])[0])
    fig = plt.figure(n)
    plt.tight_layout()
    gs = GridSpec(nrows=11, ncols=2)
    ax0 = fig.add_subplot(gs[5:7,0])
    ax1 = fig.add_subplot(gs[7:9,0])
    ax2 = fig.add_subplot(gs[9:11,0])
    ax3 = fig.add_subplot(gs[5:7,1])
    ax4 = fig.add_subplot(gs[7:9,1])
    #ax5 = fig.add_subplot(gs[7:9,1])
    ax5 = fig.add_subplot(gs[9:11,1])
    #ax7 = fig.add_subplot(gs[1:3,1])
    #ax8 = fig.add_subplot(gs[4:6,1])
    #ax9 = fig.add_subplot(gs[11:13,1])
    ax6 = fig.add_subplot(gs[0:2,:])
    ax7 = fig.add_subplot(gs[2:4,:])


    s=20

    ax5.axhspan(45,65, color='lightgrey', alpha=0.5, zorder=-1)

    for p in pattern:
        modpath='./Models/'+p
        mod=glue_runs_md(modpath)
        #all_values=glue_runs(modpath)

        crap= plt.figure(10)
        all_values=getallpars(mod, cut=(0,-30))[0]

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
        val_evol(palette, '','no', all_values['GroundinglineMassFlux'], s=s)
        ylabel('GL Mass Flux \n [km\u00b3/yr]')
        plt.sca(ax1)
        val_evol(palette,'', 'no',all_values['dGL'],s=s)
        ylabel('dGL [m/yr]')
        plt.sca(ax2)
        val_evol(palette, '','no',all_values['GLvel'], s=s)
        ylabel('$\mathregular{V_{GL}}$ [m/yr]')
        xlabel('Years')
        plt.sca(ax3)
        val_evol(palette,'','no', ((np.array(all_values['TotalCalvingFluxLevelset'])/917)*31536000)/1e9, s=s)
        ylabel('Calving Flux \n [km\u00b3/yr]')
        #xlabel('Years')
        plt.sca(ax4)
        val_evol(palette, '','no',np.array(all_values['IceVolume'])/1e9, s=s)
        ylabel('Volume [km\u00b3]')
        #plt.sca(ax5)
        #val_evol(palette, '','no',np.array(all_values['shelf_length'])/1e3, s=s)
        #ylabel('Shelf Length \n [km]')
        plt.sca(ax5)
        val_evol(palette, '','no',np.array(all_values['GLval'])/1000, s=s)
        ylabel('GL Position [km]')
        xlabel('Years')
        plt.sca(ax6)
        plot(fj_chars['P'][1], np.array(fj_chars['P'][0])/1e6, color=palette)
        xlim(20000,85000)
        ylim(1,2.5)
        ylabel('P [km\u00b2]')
        plt.sca(ax7)
        plot(np.array(fj_chars['dP'][1])/1e3, fj_chars['dP'][0], color=palette)
        ylabel('dP [m\u00b2]')
        xlabel('Distance [km]')
        xlim(20,85)
        ylim(-200,200)

    #ax9.remove()


    ax0.set_xticklabels([])
    ax1.set_xticklabels([])
    #ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax4.set_xticklabels([])
    #ax5.set_xticklabels([])
    ax6.set_xticklabels([])

    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    ax5.yaxis.tick_right()
    #ax6.yaxis.set_label_position("right")
    #ax6.yaxis.tick_right()
    #ax7.yaxis.set_label_position("right")
    #ax7.yaxis.tick_right()

    for a in gcf().axes:
        a.yaxis.set_major_locator(plt.MaxNLocator(4))
    subplots_adjust(hspace=0.5)



