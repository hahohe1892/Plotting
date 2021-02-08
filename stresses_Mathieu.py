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

#md=loadmodel('./Models/embayments/GeomProj_SpinUp_load_SG_spcvx50_NL-450_FrM1200_FlMreal180_FC40_FT200_TS0.01_OF100_hmin350_BuH0_BuP0_BuS0_ByH450_ByP55000_ByS20000_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT825_inH0_smb55_smbPos30000_funnel300_FullMelt.nc')

#md=loadmodel('./Models/depressions/GeomProj_extenddomain_SG_spcvx50_NL-450_FrM1000_FlMreal150_FC40_FT300_TS0.01_OF100_hmin350_BuH-240_BuP55000_BuS20000_ByH0_ByP0_ByS0_Stallo_IF60000_MS1000000_MSF200000_xdim85000_accT1225_inH0_smb55_smbPos30000_funnel300_FullMelt_noCutoff.nc')

interest=[ '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH1350*ByP55000*ByS20000*off.nc', '*FrM1200*FlMreal180*BuH0*BuP0*BuS0*ByH-900*ByP55000*ByS20000*.nc', '*FrM800*FlMreal120*BuH-360*BuP55000*BuS20000*ByH0*ByP0*ByS0*off.nc', '*FrM1200*FlMreal180*BuH240*BuP55000*BuS20000*ByH0*ByP0*ByS0*off.nc']

time_instance=[120,2,100,120]

geompath=['embayments/', 'bottlenecks/','depressions/', 'bumps/']

mds=[]

top = cm.get_cmap('Blues_r', 128)
bottom = cm.get_cmap('Oranges', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

ax=[]
fig = plt.figure(2)
gs = GridSpec(nrows=25, ncols=7)
ax.append(fig.add_subplot(gs[3:12,0:3]))
ax.append(fig.add_subplot(gs[3:12,4:7]))
ax.append(fig.add_subplot(gs[13:22,0:3]))
ax.append(fig.add_subplot(gs[13:22,4:7]))
countr=0
countc=0

for i,q in enumerate(interest):
    md=glue_runs_md('./Models/'+geompath[i]+q)
#for i,md in enumerate(mds):
    plt.sca(ax[i])    
    triangles=mpl.tri.Triangulation(md.mesh.x, md.mesh.y, md.mesh.elements-1)
    x=(md.mesh.x[md.mesh.elements[:,0]-1]+md.mesh.x[md.mesh.elements[:,1]-1]+md.mesh.x[md.mesh.elements[:,2]-1])/3
    y=(md.mesh.y[md.mesh.elements[:,0]-1]+md.mesh.y[md.mesh.elements[:,1]-1]+md.mesh.y[md.mesh.elements[:,2]-1])/3

    sigma_along, sigma_across, sigma_shear = stresses(md, md.results.TransientSolution[time_instance[i]].Vx, md.results.TransientSolution[time_instance[i]].Vy)
    sigma_along_a = averaging(md, sigma_along, 0)
    sigma_across_a = averaging(md, sigma_across, 0)
    sigma_shear_a = averaging(md, sigma_shear, 0) 
    #tripcolor(triangles, abs((sigma_shear).flatten()), vmin=vmin, vmax=vmax, cmap='Blues') #abs vor sigma_along
    mask=np.zeros(len(sigma_along_a))
    thk=(md.results.TransientSolution[time_instance[i]].Thickness>10)[:,0]
    if i < 2:
        mask[thk]=(sigma_across_a[thk])
        #mask[thk]=md.results.TransientSolution[time_instance[i]].Vel.flatten()[thk]
        vmin=-250000
        vmax=250000
        cmap=newcmp#'Blues'
#        tricontour(md.mesh.x,md.mesh.y,abs(mask), levels=[100000], colors=['orange'])
        ax21 = fig.add_subplot(gs[0:2,1:6])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=cm.get_cmap('Blues'), norm=norm, label='abs(Lateral shear stress)', orientation='horizontal')
    else:
        mask[thk]=(sigma_across_a[thk])
        vmin=-250000
        vmax=250000
        cmap=newcmp    
        #tricontour(md.mesh.x,md.mesh.y,abs(sigma_shear_a), levels=[25000,50000,100000,180000,200000], colors=['brown','yellow','pink', 'orange', 'white'])
        ax22 = fig.add_subplot(gs[23:,1:6])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=newcmp, norm=norm, label='Longitudinal Stress', orientation='horizontal')
    plt.sca(ax[i])
    tripcolor(md.mesh.x, md.mesh.y, mask, vmin=vmin, vmax=vmax, cmap=cmap) 
    xlim(35000,70000)
    ylim(11000,19000)
    plotcontour(md, md.results.TransientSolution[time_instance[i]].MaskGroundediceLevelset, levels=[0], colors='green')
    plotcontour(md, md.results.TransientSolution[time_instance[i]].Thickness, levels=[10], colors='black')

    mds.append(md)


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
    flotation=thickness-md.geometry.bed*((md.materials.rho_water/md.materials.rho_ice)/md.materials.rho_water)*-1
    return flotation

def find_GL(md, time):
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



    #average speed over each element
    vx_el=md.results.TransientSolution[time_instance[i]].Vx[md.mesh.elements-1]
    vx=[np.mean(x) for x in vx_el]

    vy_el=md.results.TransientSolution[time_instance[i]].Vy[md.mesh.elements-1]
    vy=[np.mean(x) for x in vy_el]

    #Get flow angle (in rad)
    theta=[]
    for a in range(len(vx)):
        theta.append(math.atan2(vy[a], vx[a]))

    #calculate stresses
    md=mechanicalproperties(md, md.results.TransientSolution[time_instance[i]].Vx, md.results.TransientSolution[time_instance[i]].Vy)

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
