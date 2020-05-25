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

bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*off.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*']
embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off.nc','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
linear=['*FrM1000*FlMreal150*BuH0*BuP0*BuS0*ByH0*ByP0*ByS0*','*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH0*ByP0*ByS0*' ]
megapat=bottlenecks+embayments+depressions+bumps

AOI=False
countg=1

for p in megapat:
    compare_path='./Figures/scaled/'+str(AOI)+'/'+p+'.pdf'
    if path.exists(compare_path):
        print('File allready exists, skipping')
        continue
    
    modpath='./Models/'+p
    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)

    z=np.nonzero(megapat==np.intersect1d(megapat, p))[0]
    if z==3 or z==4 or z==8 or z==9:
        asy='yes'
    else:
        asy='no'
    newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
    lon = np.array([int(i) for i in newstr.split()])
    wet_pat=lon[[2,4,5,7]]
    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat, GL=all_values['GLval'],asy=asy)

    if z==0:
        linear='*FrM800*FlMreal120*BuP0*BuS0*ByH0*ByP0*ByS0*'
        modpathl='./Models/'+linear
        nm=glob.glob(modpathl)
        mods=[nm[5],nm[1], nm[0], nm[2], nm[3], nm[4]]
        modl=glue_runs_md(modpathl, manual_mods=mods)
        all_valuesl=glue_runs(modpathl, manual_mods=mods)
    else:
        linear='*FrM1000*FlMreal150*BuP0*BuS0*ByH0*ByP0*ByS0*'
        modpathl='./Models/'+linear
        nm=glob.glob(modpathl)
        mods=[nm[1], nm[3], nm[0], nm[2]]
        modl=glue_runs_md(modpathl, manual_mods=mods)
        all_valuesl=glue_runs(modpathl, manual_mods=mods)
    fj_charsl, rfj_chars_GLl, rfj_chars_mvall, inds_dic_GLl, inds_dic_mvall = get_fjord(modl, all_valuesl, AOI, pattern=[0,0,0,0], GL=all_valuesl['GLval'], asy=asy)
    
    markdot=[]
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
    if z==16:
        palette='Gold'
    if z==17:
        palette='Orange'
    if z==18:
        palette='SaddleBrown'

    Fluxratio=np.array(all_values['GroundinglineMassFlux'])/np.array(all_values['IceVolume'])
    Fluxratiol=np.array(all_valuesl['GroundinglineMassFlux'])/np.array(all_valuesl['IceVolume'])
    all_values['Fluxratio']=Fluxratio
    all_valuesl['Fluxratio']=Fluxratiol
    
    for n in range(1,15):
        rs=[]
        matches=[]
        try:
            for r in range(n,len(rfj_chars_GL['x'])):
                x_match=np.squeeze(np.where(np.logical_and(rfj_chars_GL['x'][r]-201<rfj_chars_GLl['x'],rfj_chars_GL['x'][r]+201>rfj_chars_GLl['x'] )))
                try:
                    x_match=x_match[int(len(x_match)/2)]
                    x_match=int(x_match)
                except:
                    x_match=int(x_match)
                matches.append(x_match)
                rs.append(r)
            break
        except:
            print('unsuccessfull matching for first {} years'.format(r))
            continue

    #inds_dic_GL['WA']=np.array(inds_dic_GL['WA'])[np.array(inds_dic_GL['WA'])<len(all_values['dGL'])]
    #inds_dic_GL['dWA']=np.array(inds_dic_GL['dWA'])[np.array(inds_dic_GL['dWA'])<len(all_values['GroundinglineMassFlux'])]
    #par=np.array(all_values['dGL'])[inds_dic_GL['WA']]
    countt=countg+1
    for t in ['Fluxratio', 'GroundinglineMassFlux','GLvel','dGL']:
        for i in ['WA', 'dWA', 'WPer', 'dWPer']:
            fig = plt.figure(countt)
            countt+=1
            if t == 'dGL':
                rsx=rs[:-1]
                rmatches=matches[:-1]
                clean=np.array(all_values[t])[np.array(inds_dic_GL['x'])[rsx]]-np.array(all_valuesl[t])[np.array(inds_dic_GLl['x'])[rmatches]]
            else:
                clean=np.array(all_values[t])[np.array(inds_dic_GL['x'])[rs]]-np.array(all_valuesl[t])[np.array(inds_dic_GLl['x'])[matches]]
            scatter_pos(rfj_chars_GL[i][n:], 'winter', markdot, clean, marker=marker)
            title(t+' over '+i)
    countg+=16
    pp=PdfPages(compare_path)
    figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
    plt.close('all')

   
