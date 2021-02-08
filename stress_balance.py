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

def stresses(md, vx,vy):    
    #average speed over each element
    vx_el=vx[md.mesh.elements-1]
    vx_f=[np.mean(x) for x in vx_el]

    vy_el=vy[md.mesh.elements-1]
    vy_f=[np.mean(x) for x in vy_el]

    #Get flow angle (in rad)
    theta=[]
    for a in range(len(vx_f)):
        theta.append(math.atan2(vy_f[a], vx_f[a]))

    #calculate stresses
    md=mechanicalproperties(md, vx, vy)

    #Rotate reference frame
    sigma_along=np.zeros(md.mesh.numberofelements)
    sigma_across=np.zeros(md.mesh.numberofelements)
    sigma_shear=np.zeros(md.mesh.numberofelements)

    for t in range(md.mesh.numberofelements):
        R=[[cos(theta[t]), -sin(theta[t])], [sin(theta[t]), cos(theta[t])]]
        sigma=[[md.results.deviatoricstress.xx[t], md.results.deviatoricstress.xy[t]], [md.results.deviatoricstress.xy[t], md.results.deviatoricstress.yy[t]]]
        R_strich=np.linalg.inv(R)  
        sigma_new=R_strich.dot(sigma).dot(R)
        sigma_along[t]=sigma_new[0,0] #new xx component; use for basal perturbations
        sigma_across[t]=sigma_new[1,1] #new yy component
        sigma_shear[t]=sigma_new[0,1] #new xy component; use for lateral perturbations

    return sigma_along, sigma_across, sigma_shear


def height_above_flotation(md, thickness):
    flotation=np.squeeze(thickness)-md.geometry.bed*((md.materials.rho_water/md.materials.rho_ice)/md.materials.rho_water)*-1
    return flotation

def find_GL(md, t):
    ground_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset>0)
    float_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset<0)
    thk_mask=np.squeeze(md.results.TransientSolution[t].Thickness>2)
    bed_mask=np.squeeze(md.geometry.bed<0)
    no_ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset>0)
    ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset<0)
    vel_mask=np.squeeze(md.results.TransientSolution[t].Vel<1000)
    categorized=np.zeros(md.mesh.numberofvertices)
    floating=np.logical_or(np.logical_and(float_mask, thk_mask), np.logical_and(no_ice_mask, bed_mask))
    grounded=np.logical_and(np.logical_and(np.logical_and(bed_mask, thk_mask),ground_mask), ice_mask)
    categorized[floating]=-1
    categorized[grounded]=1
    categorized[vel_mask]=0

    vs=[]
    for i, q in enumerate(md.mesh.edges-1):
        voi=q[0]
        coi=categorized[voi]
        voi2=q[1]
        coi2=categorized[voi2]
        if coi == 1.0 and coi2==-1.0:
            vs.append(voi)
        if coi == -1 and coi2==1:
            vs.append(voi2)

    return vs


def streamlines(md, vx, vy, **kwargs):
    X,Y=np.meshgrid(np.arange(0,85000,10), np.arange(0,20000,10))
    VX=np.squeeze(scipy.interpolate.griddata((md.mesh.x,md.mesh.y), vx, (X,Y)))
    VY=np.squeeze(scipy.interpolate.griddata((md.mesh.x,md.mesh.y), vy, (X,Y)))
    streamplot(X,Y,VX,VY, **kwargs)
    
top = cm.get_cmap('Blues_r', 128)
bottom = cm.get_cmap('Oranges', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

interest=[ '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*off.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*.nc', '*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*off.nc', '*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*off.nc']

time_instance=[120,2,100,120]

geompath=['embayments/', 'bottlenecks/','depressions/', 'bumps/']

mds=[]

vmin=-200000
vmax=200000
for i,q in enumerate(interest):
    ax=[]
    fig = plt.figure(1)
    gs = GridSpec(nrows=25, ncols=7)
    ax.append(fig.add_subplot(gs[3:12,0:3]))
    ax.append(fig.add_subplot(gs[3:12,4:7]))
    ax.append(fig.add_subplot(gs[13:22,0:3]))
    ax.append(fig.add_subplot(gs[13:22,4:7]))
    md=glue_runs_md('./Models/'+geompath[i]+q)
    triangles=mpl.tri.Triangulation(md.mesh.x, md.mesh.y, md.mesh.elements-1)
    sigma_along, sigma_across, sigma_shear = stresses(md, md.results.TransientSolution[time_instance[i]].Vx, md.results.TransientSolution[time_instance[i]].Vy)
    ds=drivingstress(md, md.results.TransientSolution[time_instance[i]].Surface, md.results.TransientSolution[time_instance[i]].Thickness)[0] #px, py, pmag
    bd_x=basal_drag(md, md.results.TransientSolution[time_instance[i]].Thickness, md.results.TransientSolution[time_instance[i]].Base, md.results.TransientSolution[time_instance[i]].Vx)[1]
    bd_y=basal_drag(md, md.results.TransientSolution[time_instance[i]].Thickness, md.results.TransientSolution[time_instance[i]].Base, md.results.TransientSolution[time_instance[i]].Vy)[1]
    bd=basal_drag(md, md.results.TransientSolution[time_instance[i]].Thickness, md.results.TransientSolution[time_instance[i]].Base, md.results.TransientSolution[time_instance[i]].Vel)[1]
    
    sigma_along_a = averaging(md, sigma_along, 0)
    sigma_across_a = averaging(md, sigma_across, 0)
    sigma_shear_a = averaging(md, sigma_shear, 0)
    ds_a=averaging(md, ds, 0)

    thk=(md.results.TransientSolution[time_instance[i]].Thickness>10)[:,0]
    mask_shear=np.zeros(len(sigma_shear_a))
    mask_shear[thk]=abs(sigma_shear_a[thk])
    mask_stress=np.zeros(len(sigma_shear_a))
    mask_stress[thk]=abs(sigma_shear_a[thk])
    mask_ds=np.zeros(len(ds_a))
    mask_ds[thk]=ds_a[thk]
    mask_bd=np.zeros(len(bd))
    mask_bd[thk]=bd_x[thk]

    grad_longi=mask_ds-mask_shear-mask_bd


    plt.sca(ax[i])
    tripcolor(md.mesh.x, md.mesh.y, grad_longi, vmin=vmin, vmax=vmax, cmap=newcmp)

    
    
    plt.sca(ax[0])
    tripcolor(md.mesh.x, md.mesh.y, (mask_stress), cmap=newcmp, vmin=vmin, vmax=vmax)
    plt.sca(ax[1])
    tripcolor(md.mesh.x, md.mesh.y, mask_ds, cmap=newcmp, vmin=vmin, vmax=vmax) 
    plt.sca(ax[2])
    tripcolor(md.mesh.x, md.mesh.y, mask_bd, cmap=newcmp, vmin=vmin, vmax=vmax)
    plt.sca(ax[3])
    tripcolor(md.mesh.x, md.mesh.y, grad_longi, vmin=vmin, vmax=vmax, cmap=newcmp)
    for a in range(len(ax)):
        plt.sca(ax[a])
        xlim(35000,70000)
        ylim(11000,19000)
        plotcontour(md, md.results.TransientSolution[time_instance[i]].MaskGroundediceLevelset, levels=[0], colors='green')
        plotcontour(md, md.results.TransientSolution[time_instance[i]].Thickness, levels=[10], colors='black')

    ax22 = fig.add_subplot(gs[23:,1:6])
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax22, cmap=newcmp, norm=norm, label='Stress (a: shear stress, b: driving, c:basal drag, d: long gradient)', orientation='horizontal')



haf_GL=[]
for r in range(len(md.results.TransientSolution)):
    GL=find_GL(md, r)
    haf=height_above_flotation(md, md.results.TransientSolution[r].Thickness)
    haf_GL.append(np.mean(haf[GL]))




ax=[]
fig = plt.figure(2)
gs = GridSpec(nrows=27, ncols=7)
ax.append(fig.add_subplot(gs[3:12,0:3]))
ax.append(fig.add_subplot(gs[3:12,4:7]))
ax.append(fig.add_subplot(gs[13:22,0:3]))
ax.append(fig.add_subplot(gs[13:22,4:7]))

for i,q in enumerate(interest):
    md=glue_runs_md('./Models/'+geompath[i]+q)
    Vx=md.results.TransientSolution[time_instance[i]].Vx
    Vy=md.results.TransientSolution[time_instance[i]].Vy
    stress=mechanicalproperties(md, md.results.TransientSolution[time_instance[i]].Vx, md.results.TransientSolution[time_instance[i]].Vy)
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

    if i < 2:
        mask[thk]=(shear_a[thk])/1000
        vmin=-250
        vmax=250
        cmap=newcmp
#        ax21 = fig.add_subplot(gs[0:2,1:6])
#        norm = colors.Normalize(vmin=vmin, vmax=vmax)
#        cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=newcmp, norm=norm, label='Lateral shear stress [kPa]', orientation='horizontal')
    else:
        mask[thk]=(longi_a[thk])/1000
        vmin=-250
        vmax=250
        cmap=newcmp    
        ax22 = fig.add_subplot(gs[25:,1:6])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(ax22, cmap=newcmp, norm=norm, label='Longitudinal Stress Gradient/Lateral Shear Stress [kPa]', orientation='horizontal')

    plt.sca(ax[i])
    tripcolor(md.mesh.x, md.mesh.y, mask, vmin=vmin, vmax=vmax, cmap=cmap)
    if i == 0:
        seed_points=np.array([[40000,40000,40000,40000,40000],[13000,14000,15000,16000,17000]])
        streamlines(md, mask_Vx, mask_Vy, start_points=seed_points.T, color='grey')
    if i == 1:
        seed_points=np.array([[40000,40000,40000,40000,40000],[12500,14000,15000,16000,17500]])
        streamlines(md, mask_Vx, mask_Vy, start_points=seed_points.T, color='grey')
    xlim(37000,72000)
    ylim(11000,19000)
    yticks([12000,15000,18000], [2,5,8])
    xticks([40000,50000,60000,70000],[40,50,60,70])
    if i > 1:
        xlabel('x [km]')
    if i < 2:
        ax[i].set_xticklabels([])
    if i in [1,3]:
        ax[i].yaxis.set_label_position("right")
        ax[i].yaxis.tick_right()
    ylabel('y [km]')
    
    plotcontour(md, md.results.TransientSolution[time_instance[i]].MaskGroundediceLevelset, levels=[0], colors='green')
    plotcontour(md, md.results.TransientSolution[time_instance[i]].Thickness, levels=[10], colors='black')
