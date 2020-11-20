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
from scipy.stats import zscore
from socket import gethostname
from solve import *
from transientrestart2 import transientrestart
import shutil


m=glob.glob('./Jakobshavn/*.nc')

GLMF=[]
for q in m:
    md=loadmodel(q)
    
    md.miscellaneous.name='GLMF_'+q.split('/')[2]
    
    clustername=gethostname()
    cluster=generic('name', clustername, 'np', 10)
    md.cluster=cluster

    md2=deepcopy(md)
    
    for i in range(len(md.results.TransientSolution)):
        md2=transientrestart(md2, md, i)

        md2.stressbalance.requested_outputs=['GroundinglineMassFlux']
    
        md2=solve(md2, 'Stressbalance')

        ex_path=glob.glob('/home/thomas/issm/trunk-jpl/execution/GLMF*')
        shutil.rmtree(ex_path[0])
        
        GLMF.append(md2.results.StressbalanceSolution.GroundinglineMassFlux)
        print('step {} from {}'.format(i, len(md.results.TransientSolution)))


np.save('GLMF_m0.npy', GLMF)
