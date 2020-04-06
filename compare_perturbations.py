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
pattern=['*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH900*','*FrM1200*FlMreal180*ByH1350*']

AOI=True

namestring='all_comparison'
compare_path='./Figures/'+str(AOI)+'/'+namestring+str(AOI)+'.pdf'
if path.exists(compare_path):
    raise KeyError('File allready exists, skipping')

for p in pattern:
    modpath='./Models/'+p
    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)

    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)    
    markdot=[]


    z=np.nonzero(pattern==np.intersect1d(pattern, p))[0]+3

    if z==0:
        palette='LightPink'
    if z==1:
        palette='HotPink'
    if z==2:
        palette='DeepPink'
    if z==3:
        palette='SkyBlue'
    if z==4:
        palette='MediumBlue'
    if z==5:
        palette='DarkBlue'
    if z==6:
        palette='LightGreen'
    if z==7:
        palette='SeaGreen'
    if z==8:
        palette='DarkGreen'
    if z==9:
        palette='Silver'
    if z==10:
        palette='Gray'
    if z==11:
        palette='Black'
    
        ### retreat parameters over time ###
    countr=1
    for i in range(0, ceil(len(list(all_values.keys()))/3)):
        fig=plt.figure(countr)
        countr+=1
        counts=1
        for q in range(i*3, i*3+3):
            if q<len(list(all_values.keys())):
                plt.subplot(3,1,counts)
                par=all_values[list(all_values.keys())[q]]
                val_evol(palette, list(all_values.keys())[q], 'no',par)
                counts+=1
                
        ### fjord parameters ###
    countf=countr+1
    for t in range(0, ceil(len(list(fj_chars.keys()))/3)):
        fig=plt.figure(countf)
        countf+=1
        countl=1
        for r in range(t*3, t*3+3):
            if r<len(list(fj_chars.keys())):
                plt.subplot(3,1,countl)
                par=np.array(fj_chars[list(fj_chars.keys())[r]])
                plot(par[1],par[0], palette)
                greybox()
                xlim=(40000,85000)
                aoi=np.where(np.logical_and(par[1]<65000, par[1]>45000))
                plt.ylim(np.min(par[0, aoi])-np.nanstd(par[0, aoi]), np.max(par[0, aoi])+np.nanstd(par[0,aoi]))
                plt.title(list(fj_chars.keys())[r])
                grid(axis='x')
                countl+=1

    ### plot topo correlation ### 
    countt=countf+1
    for t in all_values:
        for i in rfj_chars_GL:
            fig=plt.figure(countt)
            countt+=1
            inds_dic_GL[i]=np.array(inds_dic_GL[i])[np.array(inds_dic_GL[i])<len(all_values[t])]
            par=np.array(all_values[t])[inds_dic_GL[i]]
            scatter_pos(rfj_chars_GL[i], palette, markdot, par)#all_values[t][:len(rfj_chars_GL[i])])
            plt.title(t+' over '+i)


pp=PdfPages(compare_path)
figs = [plt.figure(n) for n in plt.get_fignums()]
for fig in figs:
    fig.savefig(pp, format='pdf')
pp.close()
plt.close('all')
