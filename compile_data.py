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
import pickle

paths=[ 'embayments/', 'depressions/', 'bottlenecks/','bumps/']

embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH901*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP50000*ByS30000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']
depressions=['*FrM1200*FlMreal180*BuH-240*BuP50000*BuS30000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bottlenecks=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP50000*ByS30000*','*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*']
bumps=['*FrM1200*FlMreal180*BuH180*BuP50000*BuS30000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']

megapat=[embayments, depressions,bottlenecks, bumps]

AOI=True

one_file_val={}
one_file_cont={}

for i,path in enumerate(paths):
    geompath='./Models/'+path
    for element in megapat[i]:
        modpath=geompath+element+'.nc'
        mod_name=glob.glob(modpath)[0].split('.nc')[0]
        #if os.path.exists(mod_name+'_cont.pkl'):
        #    print(mod_name+' exists already, continue..')
        #    continue
        
        mod=glue_runs_md(modpath)
        all_values=glue_runs(modpath)

        newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in element)
        lon = np.array([int(i) for i in newstr.split()])
        wet_pat=lon[[2,4,5,7]]
        fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat)

        all_dict={'fj_chars':fj_chars,'all_values': all_values, 'rfj_chars_GL': rfj_chars_GL, 'rfj_chars_mval': rfj_chars_mval, 'inds_dic_GL': inds_dic_GL, 'inds_dic_mval': inds_dic_mval}

        #f=open(mod_name+'.pkl', 'wb')
        #pickle.dump(all_dict, f)
        #f.close()

        one_file_val[element]=all_dict
        
        cont=[]

        if hasattr(mod.results.TransientSolution[-1], 'catch_mesh'):
            lens=[]
            for q in range(len(mod.results.TransientSolution[-1].catch_mesh)):
                lens.append(len(mod.results.TransientSolution[-1].catch_mesh[q][0]))
        for q in range(len(mod.results.TransientSolution)):
            par=mod.results.TransientSolution[q].MaskGroundediceLevelset
            try:
                cont_x=tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(par), levels=[0]).allsegs[0][0][:,0]
                cont_y= tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(par), levels=[0]).allsegs[0][0][:,1]
            except:
                right_index=int(np.nonzero(lens==np.intersect1d(lens,len(par)))[0])
                cont_x=tricontour(mod.results.TransientSolution[-1].catch_mesh[right_index][0], mod.results.TransientSolution[-1].catch_mesh[right_index][1], np.squeeze(par), levels=[0]).allsegs[0][0][:,0]
                cont_y=tricontour(mod.results.TransientSolution[-1].catch_mesh[right_index][0], mod.results.TransientSolution[-1].catch_mesh[right_index][1], np.squeeze(par), levels=[0]).allsegs[0][0][:,1]

            cont.append([cont_x, cont_y])
            
        b_x1=tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(mod.geometry.bed), levels=[0]).allsegs[0][0][:,0]
        b_y1=tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(mod.geometry.bed), levels=[0]).allsegs[0][0][:,1]
        b_x2=tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(mod.geometry.bed), levels=[0]).allsegs[0][1][:,0]
        b_y2=tricontour(mod.mesh.x, mod.mesh.y, np.squeeze(mod.geometry.bed), levels=[0]).allsegs[0][1][:,1]     
        cont.append([b_x1, b_y1])
        cont.append([b_x2, b_y2])## so last element of cont is bed

        plt.gcf()
        close()
        
        #c=open(mod_name+'_cont.pkl', 'wb')
        #pickle.dump(cont,c)
        #c.close()

        one_file_cont[element]=cont
        
for q in range(len(cont)):
    plot(cont[q][0], cont[q][1])

plot(cont[-1].allsegs[0][0][:,0], cont_b.allsegs[0][0][:,1])
plot(cont[-1].allsegs[0][1][:,0], cont_b.allsegs[0][1][:,1])


v=open('OneFileVal_true.pkl', 'wb')
pickle.dump(one_file_val, v)
v.close()

c=open('OneFileCont.pkl', 'wb')
pickle.dump(one_file_cont, c)
c.close()


v=open('.pkl', 'wb')
data=pickle.dump(data_new, v)
v.close()

data_new={}
interest=[y for x in megapat for y in x]
for l,element in enumerate(data):
    data_new[interest[l]]=data[element]
