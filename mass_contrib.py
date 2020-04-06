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
#from asymetric_geometries2 import *
import os
from copy import *
from matplotlib.pyplot import tricontour
from matplotlib import cm
from matplotlib.pyplot import axvspan
from matplotlib.pyplot import grid
from plotting import *



models=glob.glob('./Models/*')
q=1
for mod in models:
    md=loadmodel(mod)
    q+=1 
    plt.figure(q)
    for i in range(len(md.results.TransientSolution)):
        calving=(md.results.TransientSolution[i].TotalCalvingFluxLevelset/917)*md.constants.yts
        floating=md.results.TransientSolution[i].TotalFloatingBmb*1e09
        smb=md.results.TransientSolution[i].TotalSmb*1e9
        influx=(np.mean(md.results.TransientSolution[i].Thickness[np.where(md.mesh.x==0)])*50*10000)
        if i>0:
            dV=md.results.TransientSolution[i].IceVolume-md.results.TransientSolution[i-1].IceVolume
        if i==0:
            dV=0
        grounded=-(calving+floating)+smb+influx-dV
        plot(i,grounded/(calving+floating),'rd')
