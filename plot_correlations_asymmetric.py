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


embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH901*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*']
bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*']




megapat=[embayments,bottlenecks]

paths=[ 'embayments/', 'bottlenecks/']

kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']

AOI=True

markdot=[]

ax=[]
fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2)
ax.append(fig.add_subplot(gs[0,0]))
ax.append(fig.add_subplot(gs[0,1]))
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position("right")

ax2=[]
fig2 = plt.figure(2)
gs2 = GridSpec(nrows=1, ncols=2)
ax2.append(fig2.add_subplot(gs2[0,0]))
ax2.append(fig2.add_subplot(gs2[0,1]))
ax2[1].yaxis.tick_right()
ax2[1].yaxis.set_label_position("right")

ax3=[]
fig3 = plt.figure(3)
gs3 = GridSpec(nrows=1, ncols=2)
ax3.append(fig3.add_subplot(gs3[0,0]))
ax3.append(fig3.add_subplot(gs3[0,1]))
ax3[1].yaxis.tick_right()
ax3[1].yaxis.set_label_position("right")


for pattern in megapat:
    n=int(np.nonzero([x==pattern for x in megapat])[0])

    geompath='./Models/'+paths[n]
    
    for p in pattern:

        modpath=geompath+p
        mod=glue_runs_md(modpath)
        all_values_long=glue_runs(modpath)
        all_values=cutallpars(all_values_long, kw, 25000)

        r=np.nonzero(pattern==np.intersect1d(pattern, p))[0]

        newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
        lon = np.array([int(i) for i in newstr.split()])
        wet_pat=lon[[2,4,5,7]]
        if r>2:
            asy='yes'
        else:
            asy='no'
        fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat, asy=asy)
        
        z=n*5+r

        if z==0:
            palette='SkyBlue'
        if z==1:
            palette='MediumBlue'
        if z==2:
            palette='DarkBlue'
        if z==3:
            palette='LightPink'
        if z==4:
            palette='Violet'
        if z==5:
            palette='Gold'
        if z==6:
            palette='Orange'
        if z==7:
            palette='Saddlebrown'
        if z==8:
            palette='LightPink'
        if z==9:
            palette='Violet'



        ### GL Flux over dWA ### 
        plt.sca(ax[n])
        
        xlabel('dS [m\u00b2/100m]')
        ylabel('$\mathregular{Q_{GL}}$ [km\u00b3/yr]')
        
        inds_dic_GL['dWA']=np.array(inds_dic_GL['dWA'])[np.array(inds_dic_GL['dWA'])<len(all_values['GroundinglineMassFlux'])]
        par=np.array(all_values['GroundinglineMassFlux'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'],par,color=palette)
        grid(True)

        ### GLvel over dWA ### 
        plt.sca(ax2[n])
        
        xlabel('dS [m\u00b2/100m]')
        ylabel('$\mathregular{V_{GL}}$ [m/yr]')
        
        par=np.array(all_values['GLvel'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'], par,color=palette)
        grid(True)        

        ### dGL over WA ###
        plt.sca(ax3[n])
        
        xlabel('S [km\u00b2]')
        ylabel('dGL [m/yr]')
        
        inds_dic_GL['WA']=np.array(inds_dic_GL['WA'])[np.array(inds_dic_GL['WA'])<len(all_values['dGL'])]
        par=np.array(all_values['dGL'])[inds_dic_GL['WA']]
        plt.scatter(np.array(rfj_chars_GL['WA'])/1e6,par,color=palette)
        grid(True)




        
