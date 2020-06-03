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
from basaldrag import *
from drivingstress import *
from slope import *
from matplotlib import colors


bottlenecks=['*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*off.nc']# '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*']
embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off.nc','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']#, '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*']
depressions=['*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bumps=['*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*MSF200000*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']



megapat=[bottlenecks, embayments, depressions, bumps]

AOI=False

modpath='./Models/'+p
mod=glue_runs_md(modpath)
all_values=glue_runs(modpath)
#all_values=getallpars(mod, cut='all')[0]

newstr = ''.join((ch if ch in '0123456789-' else ' ') for ch in p)
lon = np.array([int(i) for i in newstr.split()])
wet_pat=lon[[2,4,5,7]]
fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval = get_fjord(mod, all_values, AOI, pattern=wet_pat)

drag=[]
for i in range(len(all_values['GLval'])):
    GL=all_values['GLval'][i]
    Fr=all_values['mval'][i]
    bd=basal_drag(mod, mod.results.TransientSolution[i].Thickness, mod.results.TransientSolution[i].Base, mod.results.TransientSolution[i].Vx)[1]
    area=np.where(np.logical_and(mod.mesh.x<(GL+(Fr-GL)), mod.mesh.x>GL-10000))
    thk=np.where(mod.results.TransientSolution[i].Thickness>5)[0]
    aoi=np.intersect1d(area, thk)
    drag.append(np.sum(abs(bd[area])))
    
