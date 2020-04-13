import glob
from matplotlib.backends.backend_pdf import PdfPages
from loadmodel import *
import numpy as np
from plotmodel import *
from GetAreas import *
from scipy.optimize import curve_fit
import gc

def generateEXP(params, x_dim, y_dim):
    e=open('./exp_file2.exp', 'w')
    
    x_coordinates=[0, int(params['bay_pos1']-(params['bay_spread1'])/2),int(params['bay_pos1']-(params['bay_spread1'])/2), int(params['bay_pos1']+(params['bay_spread1'])/2), int(params['bay_pos1']+(params['bay_spread1'])/2), x_dim, x_dim, int(params['bay_pos2']+(params['bay_spread2'])/2), int(params['bay_pos2']+(params['bay_spread2'])/2), int(params['bay_pos2']-(params['bay_spread2'])/2), int(params['bay_pos2']-(params['bay_spread2'])/2), 0, 0]
    y_coordinates=[10000, 10000, 10000-max(params['bay_height1'],0), 10000-max(params['bay_height1'],0), 10000, 10000, 10000+y_dim, 10000+y_dim, 10000+y_dim+max(params['bay_height2'],0), 10000+y_dim+max(params['bay_height2'],0), 10000+y_dim, 10000+y_dim, 10000]

    #x_coordinates=[0,x_dim, x_dim, 0,0]
    #y_coordinates=[10000-max(params['bay_height1'],0), 10000-max(params['bay_height1'],0), 10000+y_dim+max(params['bay_height1'],0), 10000+y_dim+max(params['bay_height1'],0), 10000-max(params['bay_height1'],0)]
    
    p=1
    tuple_list=[]
    for i in range(0, len(x_coordinates)):
        check=(x_coordinates[i], y_coordinates[i])
        if (check in tuple_list)==True:
            continue
        else:
            tuple_list.append(check)
            p+=1

    e.write('## Name:Square\n## Icon:0\n# Points Count  Value\n{} 1.\n# X pos Y pos\n'.format(p))

    tuple_list=[]
    for i in range(0, len(x_coordinates)):
        check=(x_coordinates[i], y_coordinates[i])
        if (check in tuple_list)==True:
            continue
        else:
            tuple_list.append(check)
            e.write('{} {} \n'.format(x_coordinates[i], y_coordinates[i]))
    e.write('0 {}'.format(10000))
    e.close()

    
 
def writelog(params, run_name, x_dim, y_dim, dc, gap_halfwidth, slope, steepness, step, clusterident, output_file, status):
    f=open(output_file, 'a')
    fill_string='%s,'*28+'%s\n'
    f.write(fill_string%(run_name,' ', str(params['spcvx']), str(params['null_level']), str(params['final_time']), str(params['timestepping']), str(params['output_frequency']), str(params['hmin']), str(params['frontal_melt']), str(params['floating_melt']),str(params['friction']),str(params['bump_height']), str(params['bump_pos']), str(params['bump_spread']), str(params['bay_height1']), str(params['bay_pos1']), str(params['bay_spread1']), str(params['bay_height2']), str(params['bay_pos2']), str(params['bay_spread2']), str(x_dim),str(y_dim), str(dc), str(gap_halfwidth), str(slope), str(round(params['steepness'],3)), str(step), clusterident, status))
    f.close()



def standardvalues():
    y_dim=10000
    x_dim=100000
    slope=-2./1000
    dc=2000
    gap_halfwidth=2500
    step=100

    return y_dim, x_dim, slope, dc, gap_halfwidth, step

def plotPDF(pattern, retrieve, name):
    all_model=glob.glob('./Models/'+pattern)
    pattern=pattern[1:].split('_')[0]
    pp=PdfPages('runplots_%s.pdf'%(name))
    for i in range(0,len(all_model)):
        md=loadmodel(all_model[i])
        r=np.linspace(0,len(md.results.TransientSolution),10)
        r = [ int(x) for x in r ]
        on_off={'on': md.miscellaneous.name, 'off': all_model[i]}
        if max(r)<=10:
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Vel, 'data', md.results.TransientSolution[-1].Vel, 'colorbar', 'off', 'title',  on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1]+'_Vel'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Thickness, 'data', md.results.TransientSolution[-1].Thickness, 'colorbar', 'off', 'title', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1] +'_Thk'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].IceMaskNodeActivation, 'data', md.results.TransientSolution[-1].IceMaskNodeActivation, 'colorbar', 'off', 'title', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1] +'_IceMask'))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].MaskGroundediceLevelset, 'data', md.results.TransientSolution[-1].MaskGroundediceLevelset,'log',10,'log',10,'log',10,'log',10, 'log',10, 'log',10,'log',10,'log',10,'log', 10, 'colorbar', 'off', 'title', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1] +'_GroundMask'))
        else:
            
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[r[0]].Vel,'data', md.results.TransientSolution[r[1]].Vel, 'data',md.results.TransientSolution[r[3]].Vel,'data',md.results.TransientSolution[r[4]].Vel, 'data',md.results.TransientSolution[r[5]].Vel, 'data',md.results.TransientSolution[r[6]].Vel,'data',md.results.TransientSolution[r[7]].Vel, 'data',md.results.TransientSolution[r[8]].Vel, 'data',md.results.TransientSolution[-1].Vel, 'colorbar', 'off', 'title#2',  on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1]+'\n_Vel(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].Thickness, 'data',md.results.TransientSolution[r[1]].Thickness, 'data',md.results.TransientSolution[r[3]].Thickness,'data',md.results.TransientSolution[r[4]].Thickness, 'data',md.results.TransientSolution[r[5]].Thickness, 'data',md.results.TransientSolution[r[6]].Thickness, 'data',md.results.TransientSolution[r[7]].Thickness, 'data',md.results.TransientSolution[r[8]].Thickness,'data',md.results.TransientSolution[-1].Thickness, 'colorbar', 'off', 'title#2', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1]+'\n_Thk(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[1]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[3]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[4]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[5]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[6]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[7]].IceMaskNodeActivation, 'data',md.results.TransientSolution[r[8]].IceMaskNodeActivation,'data',md.results.TransientSolution[-1].IceMaskNodeActivation, 'colorbar', 'off', 'title#2', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1]+'\n_IceMask(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))
            pp.savefig(plotmodel(md, 'data', md.results.TransientSolution[0].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[1]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[3]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[4]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[5]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[6]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[7]].MaskGroundediceLevelset, 'data',md.results.TransientSolution[r[8]].MaskGroundediceLevelset,'data',md.results.TransientSolution[-2].MaskGroundediceLevelset,'log',10,'log',10,'log',10, 'log',10, 'log',10,'log',10,'log',10,'log', 10,'log',10, 'colorbar', 'off', 'title#2', on_off[retrieve].split('hmin')[0]+'\n'+on_off[retrieve].split('hmin')[1]+'\n_GroundMask(%s%s%s%s%s%s%s)'%(r[1], r[3], r[4],r[5], r[6], r[7],r[8])))

        fig=plt.figure(10)
        for t in range(0, len(md.results.TransientSolution)):
            along(md,md.results.TransientSolution[t].Surface)
            along(md, md.results.TransientSolution[t].Base)
            #along(md, md.results.TransientSolution[t].Vel)
        pp.savefig(plt.figure(10))
        plt.close(10)

        mod=md
        ice_Volume=[]
        floating_area=[]
        grounded_area=[]
        calving=[]
        fV=[]
        mval=[]
        mvel=[]
        mthk=[]
        mvy=[]
        meanvy=[]
        sursteep=[]
        for q in range(0, (len(mod.results.TransientSolution))):
            ice_Volume.append(mod.results.TransientSolution[q].IceVolume)
            grounded_area.append(mod.results.TransientSolution[q].GroundedArea)
            floating_area.append(mod.results.TransientSolution[q].FloatingArea)
            #calving.append(mod.results.TransientSolution[q].TotalCalvingFluxLevelset)
            fV.append(mod.results.TransientSolution[q].IceVolume-mod.results.TransientSolution[q].IceVolumeAboveFloatation)
            mval.append(np.max(mod.mesh.x[np.where(mod.results.TransientSolution[q].Thickness>40)[0]]))
            #mvel.append(np.max(mod.results.TransientSolution[q].Vel))
            mthk.append(np.max(mod.results.TransientSolution[q].Thickness))
            mvy.append(np.max(abs(mod.results.TransientSolution[q].Vy)))
            meanvy.append(np.mean(mod.results.TransientSolution[q].Vy))
        for x in range(4,(len(mod.results.TransientSolution))):
            calving.append(mod.results.TransientSolution[x].TotalCalvingFluxLevelset)
            mvel.append(np.max(mod.results.TransientSolution[x].Vel))
            sursteep.append(np.max(mod.mesh.x[np.where(np.logical_and(mod.results.TransientSolution[x].Surface>1500, mod.results.TransientSolution[x].Thickness>5))[0]]))

                           
            
        fig=plt.figure(1, dpi=20)
        plt.plot(ice_Volume, label=i)
        fig.suptitle('Ice Volume')
        fig=plt.figure(4, dpi=20)
        plt.plot(ice_Volume, label=i)
        fig.suptitle('Ice Volume')
        fig=plt.figure(2, dpi=20)
        plt.plot(grounded_area, label=i)
        fig.suptitle('Grounded Area')
        fig=plt.figure(3, dpi=20)
        plt.plot(floating_area, label=i)
        fig.suptitle('Floating Area')
        fig=plt.figure(5, dpi=20)
        plt.plot(calving, label=i)
        fig.suptitle('Calving Flux')
        fig=plt.figure(6, dpi=20)
        plt.plot(mval, label=i)
        fig.suptitle('frontal position (mval)')
        fig=plt.figure(7, dpi=20)
        plt.plot(mvel, label=i)
        fig.suptitle('maximum vel')
        fig=plt.figure(8, dpi=20)
        plt.plot(mthk, label=i)
        fig.suptitle('maximum thickness')
        fig=plt.figure(9,dpi=20)
        plt.plot(meanvy, label=i)
        fig.suptitle('Mean Vy')

        
    plt.figure(1)
    plt.legend()
    plt.figure(2)
    plt.legend()
    plt.figure(3)
    plt.legend()
    plt.figure(4)
    plt.legend()
    plt.figure(5)
    plt.legend()
    plt.figure(6)
    plt.legend()
    plt.figure(7)
    plt.legend()
    plt.figure(8)
    plt.legend()
    plt.figure(9)
    plt.legend()
    
    pp.savefig(plt.figure(1))
    pp.savefig(plt.figure(2))
    pp.savefig(plt.figure(3))
    pp.savefig(plt.figure(4))
    pp.savefig(plt.figure(5))
    pp.savefig(plt.figure(6))
    pp.savefig(plt.figure(7))
    pp.savefig(plt.figure(8))
    pp.savefig(plt.figure(9))
    plt.close('all')
    pp.close()





def scale(data):
    output=(data-np.mean(data))/np.std(data)

    return output


def across(md, parameter, limit1, limit2):
    across=np.where(np.logical_and(md.mesh.x<limit2, md.mesh.x>limit1))
    array=np.array((md.mesh.y[across], np.squeeze(parameter[across])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    plt.plot(array[0],array[1])

def along(md, parameter, limit1=15001, limit2=14999):
    along=np.where(np.logical_and(md.mesh.y<limit1, md.mesh.y>limit2))
    array=np.array((md.mesh.x[along], np.squeeze(parameter[along])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    plt.plot(array[0], array[1])

def shape(md,limit1=-50, limit2=50):
    shape=np.where(np.logical_and(md.mesh.y>15000, np.logical_and(md.geometry.bed<limit2, md.geometry.bed>limit1)))
    array=np.array((md.mesh.x[shape], np.squeeze(md.mesh.y[shape])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    plt.plot(array[0], array[1])    




def getcontrib(md):
    areas=GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
    




def adddis(md, params):
    step=50
    bay_begin1=int(params['bay_pos1']-params['bay_spread1']/2)
    bay_end1=int(params['bay_pos1']+params['bay_spread1']/2)
    if params['bay_pos1']!=params['bay_pos2']:
        bay_begin2=int(params['bay_pos2']-params['bay_spread2']/2)
        bay_end2=int(params['bay_pos2']+params['bay_spread2']/2)
    else:
        bay_begin2=0
        bay_end2=0
    bump_begin=int(params['bump_pos']-params['bump_spread']/2)
    bump_end=int(params['bump_pos']+params['bump_spread']/2)

    if params['bay_height1']>0:
        iterator=list(range(bay_begin1, bay_end1,step))
        if params['bay_height2']>0 and bay_end2!=0:
            iterator.append(list(range(bay_begin2, bay_end2,step)))
        for i in iterator:
            uplift_area=np.where(np.logical_and(md.mesh.x<i, md.mesh.x>i-step))[0]
            central_area=np.where(np.logical_and(md.mesh.y<17000, md.mesh.y>13000))
            reference=np.intersect1d(central_area, uplift_area)
            for q in uplift_area:
                if md.geometry.surface[q]<np.mean(md.geometry.surface[reference]):
                    md.geometry.surface[q]=np.mean(md.geometry.surface[reference])
                    md.geometry.base[q]=md.geometry.bed[q]
                    md.geometry.thickness[q]=md.geometry.surface[q]-md.geometry.base[q]
                    md.mask.ice_levelset[q]=-1
                else:
                    continue

    if params['bay_height1']<0:
        pos=np.where(np.logical_and(md.mesh.x<bay_end1, md.mesh.x>bay_begin1))
        md.geometry.thickness[pos]=md.geometry.surface[pos]-md.geometry.base[pos]

    if params['bay_height2']<0 and bay_end2!=0:
        pos2=np.where(np.logical_and(md.mesh.x<bay_end2, md.mesh.x>bay_begin2))
        md.geometry.thickness[pos2]=md.geometry.surface[pos2]-md.geometry.base[pos2]

    if params['bump_pos']!=0:       
        pos3=np.where(np.logical_and(md.mesh.x<bump_end, md.mesh.x>bump_begin))
        md.geometry.thickness[pos3]=md.geometry.surface[pos3]-md.geometry.bed[pos3]
        
    return md
    
def addbotbump(md, params):
    bay_begin1=int(params['bay_pos1']-params['bay_spread1']/2)
    bay_end1=int(params['bay_pos1']+params['bay_spread1']/2)
    if params['bay_pos1']!=params['bay_pos2']:
        bay_begin2=int(params['bay_pos2']-params['bay_spread2']/2)
        bay_end2=int(params['bay_pos2']+params['bay_spread2']/2)
    else:
        bay_begin2=0
        bay_end2=0
    pos=np.where(np.logical_and(md.mesh.x<bay_end1, md.mesh.x>bay_begin1))
    md.geometry.thickness[pos]=md.geometry.surface[pos]-md.geometry.base[pos]
    if bay_begin2!=0:
        pos2=np.where(np.logical_and(md.mesh.x<bay_end2, md.mesh.x>bay_begin2))
        md.geometry.thickness[pos2]=md.geometry.surface[pos2]-md.geometry.base[pos2]

    bump_begin=int(params['bump_pos']-params['bump_spread']/2)
    bump_end=int(params['bump_pos']+params['bump_spread']/2)
    pos3=np.where(np.logical_and(md.mesh.x<bump_end, md.mesh.x>bump_begin))
    md.geometry.thickness[pos3]=md.geometry.surface[pos3]-md.geometry.base[pos3]

    return md


##fit curve
def test_func(x,a,b,c):
    return -a*np.exp(-b*x)+c

#mval=[]
def forecast_mval(mval, testbegin=15, testend=30, extrapolationtime=150):
    x=np.array(list(range(testbegin,testend)))
    y=mval[testbegin:testend]
    p0=[30000, 0.09, 35000]
    popt, pcov = curve_fit(test_func, x, y,p0)
    x2=np.array(list(range(0,extrapolationtime)))
    plot(test_func(x2, *popt))

from exportVTK import *
def vtk(md, run_name):
    vtk_name=run_name
    exportVTK(vtk_name, md, 'geometry', 'mesh','mask',  'stressbalance', 'materials', 'friction', 'calving')


def wetted(md, limit1, limit2):
    across=np.where(np.logical_and(np.logical_and(md.mesh.x<limit2, md.mesh.x>limit1), md.geometry.bed<0))
    array=np.array((md.mesh.y[across], np.squeeze(md.geometry.bed[across])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    increment=[]
    centerval=[]
    for i in range(0, len(array[0])-1):
        increment.append(array[0][i+1]-array[0][i])
        centerval.append((array[1][i+1]+array[1][i])/2)
    wetted_area=np.sum(np.array(centerval)*np.array(increment))

    return wetted_area

def wetted_per(md, limit1, limit2):
    across=np.where(np.logical_and(np.logical_and(md.mesh.x<limit2, md.mesh.x>limit1), md.geometry.bed<0))
    array=np.array((md.mesh.y[across], np.squeeze(md.geometry.bed[across])))
    ind=np.argsort(array[0])
    array=array[:,ind]
    increment=[]
    for i in range(0, len(array[0])-1):
        increment.append(array[0][i+1]-array[0][i])
    wp=np.sum(increment)
    return wp

def findfluxgate(md,params):
    bump_height=np.linspace(0,-300,10)
    bump_pos=30000
    bump_spread=np.linspace(5000,40000,10)
    bump_skew=0
    bay_height1=np.linspace(3000, -1000,10)
    bay_pos1=30000
    bay_spread1=np.linspace(5000,40000,10)
    bay_skew1=0
    #bay_height2=np.linspace(3000,-1000,3)
    #bay_pos2=30000
    #bay_spread2=np.linspace(5000,40000,3)
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    start_icefront=params['start_icefront']
    slab_thickness=params['slab_thickness']
    steepness=1./300
    null_level=-500
    x_dim=55000
    params['x_dim']=55000
    params['bay_height1']=4000
    params['bay_height2']=4000
    params['bay_pos1']=30000
    params['bay_pos2']=30000
    params['bump_pos']=30000
    params['bay_spread1']=40000
    params['bay_spread2']=40000
    params['bump_spread']=40000
    generateEXP(params, x_dim, y_dim) 
    md=bamg(model(), 'domain','./exp_file2.exp', 'hmax',100, 'hmin', 100)
    #result=np.ones((len(bump_height)*len(bay_height1),3))
    result=[]
    bumps=[]
    bays=[]
    good_pars=[]
    diflist=[]
    wetstep=200
    for i in range(0, len(bump_height)):
        for q in range(0, len(bay_height1)):
            #for u in range(0, len(bump_skew)):
            for r in range(0, len(bump_spread)):
                for z in range(0, len(bay_spread1)):
                    if bay_spread1[z]/bay_height1[q]<5:
                        continue
                        #for y in range(0, len(bay_skew1)):
                    md.geometry.bed=Geometry(md, y_dim+20000, x_dim, slope, bump_spread[r],  bump_height[i], bump_pos, bump_skew, steepness, gap_halfwidth, dc, bay_spread1[z], bay_height1[q], bay_pos1,bay_skew1, bay_spread1[z], bay_height1[q], bay_pos1, bay_skew1, step)
                    old_mesh_elements=md.mesh.elements
                    old_mesh_x=md.mesh.x
                    old_mesh_y=md.mesh.y
                    old_mesh_geometry=md.geometry.bed

                    h=np.nan*np.ones(md.mesh.numberofvertices)
                    h[np.where(np.logical_and(np.logical_and(md.mesh.y<(17800+bay_spread1[z]), md.mesh.y>(12200-bay_spread1[z])), np.logical_and(md.mesh.x<52000, md.mesh.x>8000)))]=200

                    md=bamg(md, 'field', old_mesh_geometry, 'hmax', 2000, 'hmin', params['hmin'], 'gradation', 1.5, 'hVertices', h)
                    md.geometry.bed=InterpFromMeshToMesh2d(old_mesh_elements, old_mesh_x, old_mesh_y, old_mesh_geometry, md.mesh.x, md.mesh.y)[0][:,0]+null_level

                    bumps.append(bump_height[i])
                    bays.append(bay_height1[q])
                    result.append(wetted(md, 30000,30500))

                    #control_gate=[i > -2100000 and i<-1900000 for i in result]
                    wetarea=[]
                    for t in range(int(bump_pos-(max(bay_spread1[z], bump_spread[r])/2)-1000), int(bump_pos+(max(bay_spread1[z],bump_spread[r])/2)+1000),wetstep):
                        wetarea.append(wetted(md, t, t+wetstep))
                    reference=(mean(wetarea[0:5])-mean(wetarea[(len(wetarea)-5):(len(wetarea))]))/-len(wetarea)*np.array(list(range(0, len(wetarea))))+mean(wetarea[0:5])
                    percentdif=(max(wetarea)-min(wetarea))/wetarea[0]
                    referencedif=mean(abs((wetarea-reference)))#)/len(wetarea)
                    diflist.append(percentdif)
                    if referencedif< 300000 and bump_spread[r]!=0 and bump_height[i]!=0 and bay_height1[q]!=0 and bay_spread1[z]!=0 and bay_spread1[z]/bay_height1[q]>4: #wetted(md, 30000,30500)<-1800000 and wetted(md, 30000,30500)>-2200000:#percentdif>-0.5:# and wetted(md, 30000,30500)<-1800000 and wetted(md, 30000,30500)>-2200000:
                        good_pars.append([i,q,r,z,percentdif,referencedif, bump_height[i], bay_height1[q], bump_spread[r], bay_spread1[z]])
                    gc.collect()

                    if wetted(md, 30000, 30500) >-2200000 and wetted(md, 30000,30500)<-1800000:
                        good_pars.append([i, q, q])

                    bump_height[i], bump_spread[r], bay_spread1[z], bay_height1[q]
    ind=[i for i, x in enumerate(control_gate) if x]
    
def findgateman(BuH, BuS, ByH, ByS):
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    x_dim=85000
    params={}
    params['bay_height1']=params['bay_height2']=bay_height1=bay_height2=ByH
    params['bay_spread1']=params['bay_spread2']=bay_spread1=bay_spread2=ByS
    params['bump_height']=bump_height=BuH
    params['bump_spread']=bump_spread=BuS
    params['bump_pos']=bump_pos=50000
    params['bay_pos1']=params['bay_pos2']=bay_pos1=bay_pos2=50000
    bay_skew1=0
    bump_skew=0
    generateEXP(params, x_dim, y_dim)
    md=bamg(model(), 'domain','./exp_file2.exp', 'hmax',100, 'hmin', 100)
    steepness=1./300
    null_level=-450
    acc_pos=30000
    funnel=300
    md.geometry.bed=Geometry(md, y_dim+20000, x_dim, slope, bump_spread,  bump_height, bump_pos, bump_skew, steepness, gap_halfwidth, dc, bay_spread1, bay_height1, bay_pos1,bay_skew1, bay_spread1, bay_height1, bay_pos1, bay_skew1, step, acc_pos, funnel)

    wetstep=150
    p_list=[]
    for r in range(0,int(np.max(md.mesh.x)),wetstep):
        p=wetted_per(md, r, r+wetstep)
        p_list.append(p)

def wets(x,bump_height, bump_spread, bay_height1, bay_spread1):
    y_dim, x_dim, slope, dc, gap_halfwidth, step = standardvalues()
    bay_pos1=55000
    bump_pos=55000
    acc_pos=30000
    funnel=300
    null_level=-450
    steepness=1./300
    bay=0
    bump=0
    acc=0
    y=np.linspace(0,20000,20000)
    if (x > bay_pos1-bay_spread1/2 and x < bay_pos1+bay_spread1/2):
        bay=-np.sin((x-bay_spread1/4-bay_pos1)*2*np.pi/bay_spread1)*.5*bay_height1+.5*bay_height1
    if (x > bump_pos-bump_spread/2 and x < bump_pos+bump_spread/2):
        bump=-np.sin((x-bump_spread/4-bump_pos)*2*np.pi/bump_spread)*.5*bump_height+.5*bump_height
    if x<acc_pos:
        acc=((x-acc_pos)/funnel)**2


    Bx=null_level+slope*x
    By=(((dc-bump)/(1+np.exp(steepness*((y-5000)-y_dim/2+gap_halfwidth+bay+acc)))+(dc-bump)/(1+np.exp(-steepness*((y-5000)-y_dim/2-gap_halfwidth-bay-acc))))+Bx+bump)

    q=np.where(By<0)[0]
    y1=(((dc-bump)/(1+np.exp(steepness*((q-5000)-y_dim/2+gap_halfwidth+bay+acc)))+(dc-bump)/(1+np.exp(-steepness*((q-5000)-y_dim/2-gap_halfwidth-bay-acc))))+Bx+bump)
    r=q+1
    y2=(((dc-bump)/(1+np.exp(steepness*((r-5000)-y_dim/2+gap_halfwidth+bay+acc)))+(dc-bump)/(1+np.exp(-steepness*((r-5000)-y_dim/2-gap_halfwidth-bay-acc))))+Bx+bump)
    dy=y1-y2
    ds=(np.sqrt(1+(dy**2)))
    wetted_per=np.sum(ds)

    wetted_area=np.sum(By[np.where(By<0)])*-1

    return wetted_per, wetted_area, By, ds


