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

a=open('OneFileVal_true.pkl', 'rb')
all_vals=pickle.load(a)
a.close()

interest=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']

embayments=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH901*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP50000*ByS30000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP55000*ByS20000*off*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*']
depressions=['*FrM1200*FlMreal180*BuH-240*BuP50000*BuS30000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH-120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1000*FlMreal150*BuH-240*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*']
bottlenecks=['*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-1800*ByP55000*ByS20000*asy*', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP50000*ByS30000*','*FrM800*FlMreal120*BuH0*BuP0*BuS0*ByH-450*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP55000*ByS20000*','*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*']
bumps=['*FrM1200*FlMreal180*BuH180*BuP50000*BuS30000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH120*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH180*BuP55000*BuS20000*ByH0*ByP0*ByS0*','*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*']

colors=['red', 'Gold','lightblue','cornflowerblue', 'DarkBlue']
colors=['lightblue','cornflowerblue','cornflowerblue','lightblue','cornflowerblue','DarkBlue','mediumseagreen','palegreen','mediumseagreen','DarkGreen','Gold','Saddlebrown','Orange','Gold','Orange','Saddlebrown','Gray','lightgray','Gray','Black']
markers=['*','*','+','o','o','o','+','o','o','o','*','*','+','o','o','o','+','o','o','o']
labels=['Embayment', 'Depression', 'Bottleneck', 'Bump']
intercepts=[]

megapat=[embayments,depressions,bottlenecks,bumps]
interest=[y for x in megapat for y in x]

allpoint=[[],[]]

fig, ax =plt.subplots(figsize=(8,7))
for i,element in enumerate(interest):
    all_values=all_vals[element]['all_values']
    inds_dic_GL=all_vals[element]['inds_dic_GL']
    rfj_chars_GL=all_vals[element]['rfj_chars_GL']

    dP=np.array(rfj_chars_GL['dP'])
    P=np.array(rfj_chars_GL['P'])
    inds_ind=np.array(inds_dic_GL['P'])[np.array(inds_dic_GL['P'])<len(all_values['GLvel'])]
    par=np.array(all_values['GLvel'])[inds_ind]
    #if i in [6,7,8,9,16,17,18,19]:
    #    par=par*1.3
    #inds_ind2=np.array(inds_dic_GL['dP'])[np.array(inds_dic_GL['dP'])<len(all_values['dGL'])-1]
    #par2=np.array(all_values['dGL'])[inds_ind2]
    flux=np.array(all_values['Flux'])[inds_ind]
    
    fig=plt.figure(1)
    if markers[i]!='+':
        facecolor='none'
    else:
        facecolor=colors[i]
    try:
        ax.scatter(dP, par, color=colors[i], marker=markers[i], facecolors=facecolor)
    except:
        ax.scatter(dP[:-1], par, color=colors[i], marker=markers[i], facecolors=facecolor)
    xlabel('dS [m\u00b2/100m]')
    ylabel('$\it{\mathregular{v_{GL}}}$ [m/yr]')



    
    ax.scatter(dP,-P, par, color=colors[i])

    for e,t in enumerate(par):
        allpoint[1].append(t/P[e])
    for b in dP:
        allpoint[0].append(b)
    
    fig2=plt.figure(2)
    plt.scatter(dP, par2, color=colors[i])

z = np.polyfit(allpoint[0], allpoint[1],1)

p=np.poly1d(z)

o=np.linspace(-200,200,len(allpoint[0]))

plt.plot(allpoint[0], allpoint[1], '.')
plt.plot(o, p(o))

    intercepts.append(z[0])



fig, ax =plt.subplots(figsize=(8,7))
for i,element in enumerate(interest):
    if i in [6,7,8,9]:
        continue
    all_values=all_vals[element]['all_values']
    inds_dic_GL=all_vals[element]['inds_dic_GL']
    rfj_chars_GL=all_vals[element]['rfj_chars_GL']

    dP=np.array(rfj_chars_GL['dP'])
    P=np.array(rfj_chars_GL['P'])
    inds_ind=np.array(inds_dic_GL['P'])[np.array(inds_dic_GL['P'])<len(all_values['dGL'])]
    par=np.array(all_values['dGL'])[inds_ind]
    #if i in [6,7,8,9,16,17,18,19]:
    #    par=par*1.3
    #inds_ind2=np.array(inds_dic_GL['dP'])[np.array(inds_dic_GL['dP'])<len(all_values['dGL'])-1]
    #par2=np.array(all_values['dGL'])[inds_ind2]
    flux=np.array(all_values['Flux'])[inds_ind]
    
    fig=plt.figure(1)
    if markers[i]!='+':
        facecolor='none'
    else:
        facecolor=colors[i]
    try:
        ax.scatter(P/1e6, par, color=colors[i], marker=markers[i], facecolors=facecolor)
    except:
        ax.scatter(P[:-1]/1e6, par, color=colors[i], marker=markers[i], facecolors=facecolor)
    xlabel('dS [m\00b2/100m]')
    ylabel('$\it{\mathregular{_{GL}]}$ [km\u00b2]')
