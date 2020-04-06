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

pattern=['*FrM800*FlMreal120*ByH-450*','*FrM1200*FlMreal180*ByH450*','*FrM1200*FlMreal180*ByH-900*Cutoff.nc','*FrM1200*FlMreal180*ByH900*Cutoff.nc','*FrM1200*FlMreal180*ByH-675*','*FrM1200*FlMreal180*ByH1350*','*FrM1200*FlMreal180*BuH120*','*FrM1200*FlMreal180*BuH-120*','*FrM1000*FlMreal150*BuH-240*','*FrM1200*FlMreal180*BuH240*','*FrM1200*FlMreal180*BuH180*','*FrM1200*FlMreal180*BuH-360*', '*FrM1200*FlMreal180*ByH900*asy*','*FrM1200*FlMreal180*ByH-900*asy*','*FrM1200*FlMreal180*ByH-1800*asy*']

AOI=False

for p in pattern:
    modpath='./Models/'+p
    namestring='analysis'+p
    analysis_path='./Figures/'+str(AOI)+'/'+namestring+str(AOI)+'.pdf'
    if path.exists(analysis_path):
        print('Analysis allready exists, skipping')
        continue

    mod=glue_runs_md(modpath)
    all_values=glue_runs(modpath)


    fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI)
    
    plt.close('all')
    markdot=[]
    
    pp=PdfPages(analysis_path)

    cont_all(mod, ['GL'],1)
    pp.savefig()
    plt.close('all')
    cont_all(mod, ['Front'],1)
    pp.savefig()
    plt.close('all')
    cont_all(mod, ['GL'],10)
    pp.savefig()
    plt.close('all')

    fig=plt.figure(2)
    plt.subplot(2,1,1)
    along_evol(mod, 'winter', 'yes', '', 'Surface', 'Base')
    greybox()
    plt.subplot(2,1,2)
    along_evol(mod, 'winter', 'yes', '', 'Vel')
    plt.xlim(0,85000)
    greybox()
    pp.savefig(2)


    ### retreat parameters over time ###
    for i in range(0, ceil(len(list(all_values.keys()))/3)):
        fig=plt.figure()
        counts=1
        for q in range(i*3, i*3+3):
            if q<len(list(all_values.keys())):
                plt.subplot(3,1,counts)
                par=all_values[list(all_values.keys())[q]]
                val_evol('winter', list(all_values.keys())[q],'yes', par)
                counts+=1
        pp.savefig()

    ### fjord parameters ###
    for t in range(0, ceil(len(list(fj_chars.keys()))/3)):
        fig=plt.figure()
        countl=1
        for r in range(t*3, t*3+3):
            if r<len(list(fj_chars.keys())):
                plt.subplot(3,1,countl)
                par=np.array(fj_chars[list(fj_chars.keys())[r]])
                plot(par[1],par[0], 'blue')
                greybox()
                xlim=(40000,85000)
                aoi=np.where(np.logical_and(par[1]<65000, par[1]>45000))
                plt.ylim(np.min(par[0, aoi])-np.nanstd(par[0, aoi]), np.max(par[0, aoi])+np.nanstd(par[0,aoi]))
                plt.title(list(fj_chars.keys())[r])
                grid(axis='x')
                countl+=1
        pp.savefig()

    ### plot topo correlation ### 
    count=1
    for t in all_values:
        for i in rfj_chars_GL:
            fig=plt.figure()
            inds_dic_GL[i]=np.array(inds_dic_GL[i])[np.array(inds_dic_GL[i])<len(all_values[t])]
            par=np.array(all_values[t])[inds_dic_GL[i]]
            scatter_pos(rfj_chars_GL[i], 'winter', markdot, par)#all_values[t][:len(rfj_chars_GL[i])])
            plt.title(t+' over '+i)
            count+=1
            pp.savefig()
    plt.close('all')    
    pp.close()

    ### topo over time ###
    fjord_string='fjord'+p
    fjord_path='./Figures/'+str(AOI)+'/'+fjord_string+str(AOI)+'.pdf'
    if path.exists(fjord_path):
        print('Fjord allready exists, skipping')
        continue

    pp=PdfPages(fjord_path)
    for i in range(ceil(len(list(rfj_chars_GL.keys()))/3)):
        fig=plt.figure()
        countf=1
        for r in range(i*3, i*3+3):
            if r<len(list(rfj_chars_GL.keys())):
                plt.subplot(3,1,countf)
                val_evol('winter',list(rfj_chars_GL.keys())[r],'no',rfj_chars_GL[list(rfj_chars_GL.keys())[r]])
                countf+=1
        pp.savefig()
    plt.close('all')
    pp.close()
