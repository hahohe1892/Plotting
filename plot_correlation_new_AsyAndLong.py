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

embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH901*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP50000*ByS30000*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*', '*FrM1200*FlMreal180*BuH-240*BuP50000*BuS30000*ByH0*ByP0*ByS0*']

bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP50000*ByS30000*']

bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP50000*BuS30000*ByH0*ByP0*ByS0*']


megapat=[embayments, depressions,bottlenecks, bumps]

paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']

kw=['GroundinglineMassFlux', 'dGL', 'GLvel', 'TotalCalvingFluxLevelset', 'IceVolume', 'GLval']

markdot=[]

ax=[]
fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)
ax.append(fig.add_subplot(gs[0,0]))
ax.append(fig.add_subplot(gs[0,1]))
ax.append(fig.add_subplot(gs[1,0]))
ax.append(fig.add_subplot(gs[1,1]))
ax[1].yaxis.tick_right()
ax[3].yaxis.tick_right()
ax[1].yaxis.set_label_position("right")
ax[3].yaxis.set_label_position("right")

ax2=[]
fig2 = plt.figure(2)
gs2 = GridSpec(nrows=2, ncols=2)
ax2.append(fig2.add_subplot(gs2[0,0]))
ax2.append(fig2.add_subplot(gs2[0,1]))
ax2.append(fig2.add_subplot(gs2[1,0]))
ax2.append(fig2.add_subplot(gs2[1,1]))
ax2[1].yaxis.tick_right()
ax2[3].yaxis.tick_right()
ax2[1].yaxis.set_label_position("right")
ax2[3].yaxis.set_label_position("right")

ax3=[]
fig3 = plt.figure(3)
gs3 = GridSpec(nrows=2, ncols=2)
ax3.append(fig3.add_subplot(gs3[0,0]))
ax3.append(fig3.add_subplot(gs3[0,1]))
ax3.append(fig3.add_subplot(gs3[1,0]))
ax3.append(fig3.add_subplot(gs3[1,1]))
ax3[1].yaxis.tick_right()
ax3[3].yaxis.tick_right()
ax3[1].yaxis.set_label_position("right")
ax3[3].yaxis.set_label_position("right")

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
        wet_pat_all=lon[[2,3,4,5,6,7]]
        wet_pat=lon[[2,4,5,7]]
        if n == 0 or n == 2:
            if r==3 or r == 4:
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
        if z==3:
            palette='LightPink'
        if z==4:
            palette='Violet'
        if z==5:
            palette='Red'        
        if z==10:
            palette='LightGreen'
        if z==11:
            palette='SeaGreen'
        if z==12:
            palette='DarkGreen'
        if z==13:
            palette='Red'
        if z==20:
            palette='Gold'
        if z==21:
            palette='Orange'
        if z==22:
            palette='Saddlebrown'
        if z==23:
            palette='LightPink'
        if z==24:
            palette='Violet'
        if z==25:
            palette='Red'        
        if z==30:
            palette='Silver'
        if z==31:
            palette='Gray'
        if z==32:
            palette='Black'
        if z==33:
            palette='Red'

        ### GL Flux over dWA ### 
        plt.sca(ax[n])
        #xlim(-26000,26000)
        #ylim(5,14.5)
        
        if n == 2 or n == 3:
            xlabel('dS [m\u00b2/100m]')
        ylabel('$\mathregular{Q_{GL}}$ [km\u00b3/yr]')
        
        inds_dic_GL['dWA']=np.array(inds_dic_GL['dWA'])[np.array(inds_dic_GL['dWA'])<len(all_values['GroundinglineMassFlux'])]
        par=np.array(all_values['GroundinglineMassFlux'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'],par,color=palette)
        grid(True)

        ### GLvel over dWA ### 
        plt.sca(ax2[n])
        #xlim(-26000,26000)
        #ylim(1500,7500)
        
        if n == 2 or n == 3:
            xlabel('dS [m\u00b2/100m]')
        ylabel('$\mathregular{V_{GL}}$ [m/yr]')
        
        par=np.array(all_values['GLvel'])[inds_dic_GL['dWA']]
        plt.scatter(rfj_chars_GL['dWA'], par,color=palette)
        grid(True)

        ### dGL over WA ###
        plt.sca(ax3[n])
        #xlim(1,3.5)
        #ylim(-9000,500)
        
        if n == 2 or n == 3:
            xlabel('S [km\u00b2]')
        ylabel('dGL [m/yr]')
        
        inds_dic_GL['WA']=np.array(inds_dic_GL['WA'])[np.array(inds_dic_GL['WA'])<len(all_values['dGL'])]
        par=np.array(all_values['dGL'])[inds_dic_GL['WA']]
        plt.scatter(np.array(rfj_chars_GL['WA'])/1e6,par,color=palette)
        grid(True)
    
