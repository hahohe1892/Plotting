import numpy as np
from loadmodel import *
import math
from mechanicalproperties import *
from plotting import *
from matplotlib import colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from averaging import averaging
from drivingstress import drivingstress
from basaldrag import basal_drag
from slope import slope
import scipy
top = cm.get_cmap('Blues_r', 128)
bottom = cm.get_cmap('Oranges', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

interest=[ '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1800*ByP55000*ByS20000*asy*.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH900*ByP50000*ByS30000*.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-901*ByP55000*ByS20000*asy*.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-675*ByP50000*ByS30000*.nc','*FrM1200*FlMreal180*BuH-240*BuP50000*BuS30000*ByH0*ByP0*ByS0*.nc','*FrM1200*FlMreal180*BuH180*BuP50000*BuS30000*ByH0*ByP0*ByS0*.nc']

time_instance=[10,10,4,4,300,5]

geompath=['embayments/', 'embayments/','bottlenecks/','bottlenecks/','depressions/', 'bumps/']

ax=[]
fig = plt.figure()
gs = GridSpec(nrows=33, ncols=21)
ax.append(fig.add_subplot(gs[0:8,0:10]))
ax.append(fig.add_subplot(gs[0:8,11:21]))
ax.append(fig.add_subplot(gs[10:18,0:10]))
ax.append(fig.add_subplot(gs[10:18,11:21]))
ax.append(fig.add_subplot(gs[20:28,0:10]))
ax.append(fig.add_subplot(gs[20:28,11:21]))

for i,q in enumerate(interest):
    md=glue_runs_md('./Models/'+geompath[i]+q)
    Vx=md.results.TransientSolution[time_instance[i]].Vx
    Vy=md.results.TransientSolution[time_instance[i]].Vy
    stress=mechanicalproperties(md, Vx, Vy)

    sigma_xy=md.results.deviatoricstress.xy
    sigma_xy_a=averaging(md, sigma_xy, 0)
    H=np.squeeze(md.results.TransientSolution[time_instance[i]].Thickness)
    shear=slope(md, H*sigma_xy_a)[1]
    shear_a=averaging(md, shear, 0)*-1

    sigma_xx=md.results.deviatoricstress.xx
    sigma_xx_a=averaging(md, sigma_xx, 0)
    sigma_yy=md.results.deviatoricstress.yy
    sigma_yy_a=averaging(md, sigma_yy, 0)
    longi=slope(md, 2*H*sigma_xx_a+H*sigma_yy_a)[0]
    longi_a=averaging(md, longi, 0)*-1
    
    mask=np.zeros(len(shear_a))
    thk=(md.results.TransientSolution[time_instance[i]].Thickness>10)[:,0]
    mask_Vx=np.zeros(len(Vx))
    mask_Vy=np.zeros(len(Vy))
    mask_Vx[thk]=Vx[thk,0]
    mask_Vy[thk]=Vy[thk,0]

    if i < 4:
        mask[thk]=(shear_a[thk])/1000
        vmin=-250
        vmax=250
        cmap=newcmp
    else:
        mask[thk]=(longi_a[thk])/1000
        vmin=-250
        vmax=250
        cmap=newcmp    
        ax22 = fig.add_subplot(gs[31:,5:16])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(ax22, cmap=newcmp, norm=norm, label='Lateral Shear Stress/Longitudinal Stress Gradient [kPa]', orientation='horizontal')

    plt.sca(ax[i])
    tripcolor(md.mesh.x, md.mesh.y, mask, vmin=vmin, vmax=vmax, cmap=cmap)

    xlim(32000,72000)
    ylim(10500,19000)
    yticks([12000,15000,18000], [2,5,8])
    xticks([35000,45000,55000,65000],[35,45,55,65])
    if i in [4,5]:
        xlabel('x [km]')
    if i < 4:
        ax[i].set_xticklabels([])
    if i in [1,3,5]:
        ax[i].yaxis.set_label_position("right")
        ax[i].yaxis.tick_right()
    ylabel('y [km]')
    
    plotcontour(md, md.results.TransientSolution[time_instance[i]].MaskGroundediceLevelset, levels=[0], colors='green')
    plotcontour(md, md.results.TransientSolution[time_instance[i]].Thickness, levels=[10], colors='black')
    
    






    triangles=mpl.tri.Triangulation(mod.mesh.x, mod.mesh.y, mod.mesh.elements-1)
    x=(mod.mesh.x[mod.mesh.elements[:,0]-1]+mod.mesh.x[mod.mesh.elements[:,1]-1]+mod.mesh.x[mod.mesh.elements[:,2]-1])/3
    y=(mod.mesh.y[mod.mesh.elements[:,0]-1]+mod.mesh.y[mod.mesh.elements[:,1]-1]+mod.mesh.y[mod.mesh.elements[:,2]-1])/3

    
    ds=drivingstress(mod, mod.results.TransientSolution[time_instance[i]].Surface, mod.results.TransientSolution[time_instance[i]].Thickness)[0]
    bd=basal_drag(mod, mod.results.TransientSolution[time_instance[i]].Thickness, mod.results.TransientSolution[time_instance[i]].Base, mod.results.TransientSolution[time_instance[i]].Vel)[0]
    ls=ds-bd

    tripcolor(triangles, (ls).flatten(), vmin=vmin, vmax=vmax, cmap=newcmp)
    
    xlim(30000,70000)
    ylim(11000,19000)
    plotcontour(mod, mod.results.TransientSolution[time_instance[i]].MaskGroundediceLevelset, levels=[0], colors='green')
    plotcontour(mod, mod.results.TransientSolution[time_instance[i]].Thickness, levels=[10], colors='black')
    
ax21 = fig.add_subplot(gs[30:,:])
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=newcmp, norm=norm, label='Longitudinal Stress', orientation='horizontal')


    

