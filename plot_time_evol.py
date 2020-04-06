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
pattern=['*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*off.nc','*FrM1200*FlMreal180*ByH1350*']

AOI=False

fig = plt.figure()
plt.tight_layout()
gs = GridSpec(nrows=7, ncols=2)
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0])
ax2 = fig.add_subplot(gs[2,0])
ax3 = fig.add_subplot(gs[3,0])
ax4 = fig.add_subplot(gs[4,0])
ax5 = fig.add_subplot(gs[5,0])
ax6 = fig.add_subplot(gs[6,0])
ax7 = fig.add_subplot(gs[1:3,1])
ax8 = fig.add_subplot(gs[4:6,1])
ax9 = fig.add_subplot(gs[6:7,1])

s=20

for p in pattern:
    modpath='./Models/'+p
    mod=glue_runs_md(modpath)
    #all_values=glue_runs(modpath)

    all_values=getallpars(mod, cut=(0,-30))[0]

    plt.sca(ax9)
    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)

    markdot=[]


    z=np.nonzero(pattern==np.intersect1d(pattern, p))[0]

    if z==0:
        palette='SkyBlue'
    if z==1:
        palette='MediumBlue'
    if z==2:
        palette='DarkBlue'

    plt.sca(ax0)
    val_evol(palette, '','no', all_values['GroundinglineMassFlux'], s=s)
    ylabel('GL Mass Flux \n [km\u00b3/yr]')
    plt.sca(ax1)
    val_evol(palette,'', 'no',all_values['dGL'],s=s)
    ylabel('dGL [m/yr]')
    plt.sca(ax2)
    val_evol(palette, '','no',all_values['GLvel'], s=s)
    ylabel('$\mathregular{V_{GL}}$ [m/yr]')
    plt.sca(ax3)
    val_evol(palette, '','no',np.array(all_values['IceVolume'])/1e9, s=s)
    ylabel('Volume [km\u00b3]')
    plt.sca(ax4)
    val_evol(palette,'','no', ((np.array(all_values['TotalCalvingFluxLevelset'])/917)*31536000)/1e9, s=s)
    ylabel('Calving Flux \n [km\u00b3/yr]')
    plt.sca(ax5)
    val_evol(palette, '','no',np.array(all_values['shelf_length'])/1e3, s=s)
    ylabel('Shelf Length \n [km]')
    plt.sca(ax6)
    val_evol(palette, '','no',np.array(all_values['GLval'])/1000, s=s)
    ylabel('GL Position [km]')
    xlabel('Years')
    plt.sca(ax7)
    val_evol(palette, '','',np.array(rfj_chars_GL['P'])/1e6, s=s)
    ylabel('P [km\u00b2]')
    plt.sca(ax8)
    val_evol(palette, '','',rfj_chars_GL['dP'], s=s)
    ylabel('dP [m\u00b2]')
    xlabel('Years')

ax9.remove()


ax0.set_xticklabels([])
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])
ax5.set_xticklabels([])

ax7.yaxis.set_label_position("right")
ax7.yaxis.tick_right()
ax8.yaxis.set_label_position("right")
ax8.yaxis.tick_right()

for a in gcf().axes:
    a.yaxis.set_major_locator(plt.MaxNLocator(4))





