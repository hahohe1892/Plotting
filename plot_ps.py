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

fig, ax = plt.subplots(2,2)
#ax1 = ax[0].twinx()
#ax2 = ax[1].twinx()
x=np.linspace(0,85,425)
for p in megapat:

    z=np.nonzero(megapat==np.intersect1d(megapat, p))[0]

    newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
    lon = np.array([int(i) for i in newstr.split()])
    wet_pat=lon[[2,4,5,7]]
    wp=[]
    wa=[]
    #if z==3 or z==4 or z==8 or z==9:
    #    continue
    #else:
    #    asy='no'
    for i in np.linspace(0,85000,425):
        per, a = wets(i, wet_pat[0], wet_pat[1], wet_pat[2], wet_pat[3])[0:2]
        wp.append(per)
        wa.append(a)

    dwp=deltaval(wp)
    dwa=deltaval(wa)
        
    if z==0:
        palette='Gold'
        marker='-'
    if z==1:
        palette='Orange'
    if z==2:
        palette='Saddlebrown'
    if z==3:
        palette='lightblue'
        marker='--'
    if z==4:
        palette='cornflowerblue'
        marker='--'
    if z==5:
        palette='DarkBlue'
        marker='-'
    if z==6:
        palette='palegreen'
        marker='-'
    if z==7:
        palette='mediumseagreen'
    if z==8:
        palette='DarkGreen'
    if z==9:
        palette='lightgray'
    if z==10:
        palette='Gray'
    if z==11:
        palette='Black'

    if z<6:
        plt.sca(ax[0][0])
        ax[0][0].yaxis.set_label_position("left")
        ax[0][0].yaxis.tick_left()        
    else:
        plt.sca(ax[0][1])
        ax[0][1].yaxis.set_label_position("right")
        ax[0][1].yaxis.tick_right()
        
    plot(x, np.array(wa)/1e6, color=palette, linestyle='-')
    xlim(20,85)
    ylim(1,3.8)
    ylabel('S [km\u00b2]')
#    plt.sca(ax1)
#    plot(x, wp, linestyle='--', color=palette)
#    ylim(0,10000)
    if z<6:
        plt.sca(ax[1][0])
        ax[1][0].yaxis.set_label_position("left")
        ax[1][0].yaxis.tick_left()
    else:
        plt.sca(ax[1][1])
        ax[1][1].yaxis.set_label_position("right")
        ax[1][1].yaxis.tick_right()
    plot(x[:-1], np.array(dwa)*-1, color=palette, linestyle='-')
    ylabel('dS [m\u00b2]')
    xlabel('x-coordinates [km]')
    xlim(20,85)
    ylim(-55000,55000)
#    plt.sca(ax2)
#    plot(x[:-1], dwp, linestyle='--', color=palette)
#    ylim(-120,120)
#    xlim(20,85)


#ax1.yaxis.set_label_text('WP [m]')
#ax2.yaxis.set_label_text('dWP [m]')
