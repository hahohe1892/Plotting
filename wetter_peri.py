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

#pattern=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH-900*','*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*','*FrM1200*FlMreal180*ByH1350*','*FrM1200*FlMreal180*BuH120*','*FrM1200*FlMreal180*BuH180*','*FrM1200*FlMreal180*BuH240*','*FrM1200*FlMreal180*BuH-120*','*FrM1000*FlMreal150*BuH-240*','*FrM1200*FlMreal180*BuH-360*']

bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*off.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*']
embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off.nc','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']

megapat=bottlenecks+embayments+depressions+bumps

AOI=True

namestring='all_comparison_WP'
compare_path='./Figures/'+str(AOI)+'/'+namestring+str(AOI)+'.pdf'
if path.exists(compare_path):
    raise KeyError('File allready exists, skipping')

for p in megapat:
    modpath='./Models/'+p
    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)
    
    fig = plt.figure(20)
    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)
    newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
    lon = np.array([int(i) for i in newstr.split()])
    wp=[[],np.linspace(0,85000,85000)]
    #wa=[[],np.linspace(0,85000,85000)]
    for q in range(85000):
        per, area=wets(q, lon[2], lon[4], lon[5], lon[7])[0:2]
        wp[0].append(per)
     #   wa[0].append(area)
    
    fj_chars['WPer']=wp
    #fj_chars['WA']=wa
    #fj_chars['dWA']=[deltaval(wa[0]), wa[1]]
    rfj_chars_GL, inds_dic_GL=common(all_values['GLval'], fj_chars, AOI,  **fj_chars)
    markdot=[]
    plt.close(20)


    z=np.nonzero(megapat==np.intersect1d(megapat, p))[0]
    marker='o'
    
    if z==0:
        palette='LightPink'
    if z==1:
        palette='HotPink'
    if z==2:
        palette='DeepPink'
    if z==3:
        palette='LightPink'
        marker='+'
    if z==4:
        palette='DeepPink'
        marker='+'
    if z==5:
        palette='SkyBlue'
    if z==6:
        palette='MediumBlue'
    if z==7:
        palette='DarkBlue'
    if z==8:
        palette='SkyBlue'
        marker='+'
    if z==9:
        palette='DarkBlue'
        marker='+'
    if z==10:
        palette='LightGreen'
    if z==11:
        palette='SeaGreen'
    if z==12:
        palette='DarkGreen'
    if z==13:
        palette='Silver'
    if z==14:
        palette='Gray'
    if z==15:
        palette='Black'


    ### plot topo correlation ### 
    countt=1
    for t in ['dGL', 'GroundinglineMassFlux', 'GLvel', 'TotalCalvingFluxLevelset', 'max Vx']:
        for i in ['WPer']:
            fig=plt.figure(countt)
            countt+=1
            inds_dic_GL[i]=np.array(inds_dic_GL[i])[np.array(inds_dic_GL[i])<len(all_values[t])]
            par=np.array(all_values[t])[inds_dic_GL[i]]
            scatter_pos(rfj_chars_GL[i], palette, markdot, par, marker=marker)#all_values[t][:len(rfj_chars_GL[i])])
            plt.title(t+' over '+i)
