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

  

#### define functions for retreat parameters 

  

def snapshot_birdseye(md, par, val): 

    r=np.linspace(0,len(md.results.TransientSolution),10) 

    r = [ int(x) for x in r ] 

    plotmodel(md, 'data', md.results.TransientSolution[r[0]].__dict__[par],'contourlevels', [val], 'data', md.results.TransientSolution[r[1]].__dict__[par], 'contourlevels',[val],'data',md.results.TransientSolution[r[3]].__dict__[par],'contourlevels',[val],'data',md.results.TransientSolution[r[4]].__dict__[par], 'contourlevels', [val],'data',md.results.TransientSolution[r[5]].__dict__[par], 'contourlevels',[val],'data',md.results.TransientSolution[r[6]].__dict__[par],'contourlevels',[val],'data',md.results.TransientSolution[r[7]].__dict__[par],'contourlevels',[val], 'data',md.results.TransientSolution[r[8]].__dict__[par], 'contourlevels',[val],'data',md.results.TransientSolution[-1].__dict__[par],'contourlevels',[val])   

  

  

def ssop(md, *args): #scalar standard output pars 

    dic={} 

    for arg in args: 

        par_list=[] 

        for i in range(0, len(md.results.TransientSolution)): 

            if hasattr(md.results.TransientSolution[i],arg): 

                par_list.extend(md.results.TransientSolution[i].__dict__[arg]) 

                dic[arg]=par_list 

    return dic 

  

def msop(md, *args): #max standard output pars 

    dic={} 

    for arg in args:  

        par_list=[] 

        for i in range(0, len(md.results.TransientSolution)): 

            if hasattr(md.results.TransientSolution[i], arg): 

                par_list.append(np.max(abs(md.results.TransientSolution[i].__dict__[arg]))) 

                dic['max '+arg]=par_list 

    return dic 

  

  

def GLvalf(md, depression='no'): 

    GLval=[] 
    central_lowGL=np.where(np.logical_and(md.mesh.y<15500, md.mesh.y>14500))
    for q in range(0, len(md.results.TransientSolution)):
        try:
            if depression=='no':
                GL_criteria=np.where(md.results.TransientSolution[q].MaskGroundediceLevelset>0)
                GLval.append(np.max(md.mesh.x[np.intersect1d(central_lowGL, GL_criteria)]))
            else:
                GL_criteria=np.where(md.results.TransientSolution[q].MaskGroundediceLevelset>0)
                GLval.append(np.max(md.mesh.x[np.intersect1d(central_lowGL, GL_criteria)]))
        except ValueError:
            GLval.append(10000)
    return GLval 

  

def mvalf(md): 

    mval=[] 
    central_low=np.where(np.logical_and(md.mesh.y<16500, md.mesh.y>13500))
    for q in range(0, len(md.results.TransientSolution)): 
        thick_criteria=np.where(md.results.TransientSolution[q].Thickness>10)
        try:
            mval.append(np.max(md.mesh.x[np.intersect1d(central_low,thick_criteria)]))
        except ValueError:
            mval.append(10000)
    return mval 

def findGL(md,t):
    ground_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset>0)
    float_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset<0)
    thk_mask=np.squeeze(md.results.TransientSolution[t].Thickness>2)
    bed_mask=np.squeeze(md.geometry.bed<0)
    no_ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset>0)
    ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset<0)
    floating=np.logical_or(np.logical_and(float_mask, thk_mask), np.logical_and(no_ice_mask, bed_mask))
    grounded=np.logical_and(np.logical_and(np.logical_and(bed_mask, thk_mask),ground_mask), ice_mask)
    categorized=np.zeros(md.mesh.numberofvertices)
    categorized[floating]=-1
    categorized[grounded]=1
    
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
    
def GLvel(md, GLval): 

    central_lowGL=np.where(np.logical_and(md.mesh.y<15500, md.mesh.y>14500)) 

    GLvel=[] 

    GLvelx=[] 

    for i in range(len(md.results.TransientSolution)): 

        GL_area=np.where(np.logical_and(md.mesh.x<(GLval[i]+100), md.mesh.x>(GLval[i]-100))) 

        vel_area=np.intersect1d(central_lowGL, GL_area) 

        GLvel.append(np.mean(md.results.TransientSolution[i].Vel[vel_area])) 

        GLvelx.append(np.mean(md.results.TransientSolution[i].Vx[vel_area])) 

    return GLvel, GLvelx 


def GLvel_mean(md):
    mean_vels=[]
    med_vels=[]
    for i in range(len(md.results.TransientSolution)):
        mean_vels.append(np.mean(md.results.TransientSolution[i].Vel[findGL(md, i)]))
        med_vels.append(np.median(md.results.TransientSolution[i].Vel[findGL(md, i)]))

    return mean_vels, med_vels

def deltaval(val): 

    delta=[] 

    for i in range(len(val)-1): 

        delta.append(val[i+1]-val[i]) 

    return delta 

     

def getcolors(n, palette): 

    try: 

        colors = cm.__dict__[palette](np.linspace(0,1,n)) 

    except KeyError: 

        colors=[str(palette)]*n 

    return colors 

  

  

def along_evol(md, palette, bar='yes', title='', limit=-1, *args, **kwargs): 

    colors=getcolors(len(md.results.TransientSolution[:limit]), palette) 

    for arg in args: 

        for i in range(0, len(md.results.TransientSolution[:limit])): 

            along(md, md.results.TransientSolution[i].__dict__[arg], color=colors[i], **kwargs) 

    if bar=='yes': 

        sm = plt.cm.ScalarMappable(cmap=palette, norm=plt.Normalize(vmin=0, vmax=len(md.results.TransientSolution[:limit]))) 

        plt.colorbar(sm) 

    if title=='m': 

        plt.title(input('Set title: ')) 

    elif title=='no': 

        print('No title set') 

    else: 

        plt.title(args) 

  

def across_evol(md, palette, limit1, limit2, *args): 

    colors=getcolors(len(md.results.TransientSolution), palette) 

    for arg in args: 

        for i in range(0, len(md.results.TransientSolution)): 

            across(md,md.results.TransientSolution[i].__dict__[arg], limit1,limit2, color=colors[i]) 

  

def val_evol(palette, title, bar='no', *args, **kwargs): 

    for arg in args: 

        colorrange=len(arg) 

        colors=getcolors(colorrange, palette) 

        for i in range(0,len(arg)): 

            plt.scatter(range(0, len(arg))[i],arg[i], color=colors[i], **kwargs)

        lc='lightblue'
        plt.plot(arg, color=lc) 

    if bar=='yes': 

        sm = plt.cm.ScalarMappable(cmap=palette, norm=plt.Normalize(vmin=0, vmax=colorrange)) 

        plt.colorbar(sm) 

    plt.title(title) 

    plt.grid(True) 

  

def scatter_pos(GL, palette, markdot, *args, **kwargs): 

    colors=getcolors(len(GL), palette) 

    for arg in args: 

        for i in range(len(arg)): 

            if i in markdot: 

                plt.scatter(GL[i], arg[i], color='red', **kwargs) 

            else: 

                plt.scatter(GL[i], arg[i], color=colors[i], **kwargs) 

    grid(True) 

  

def getbar(palette, arg, **kwargs): 

    sm = plt.cm.ScalarMappable(cmap=palette, norm=plt.Normalize(vmin=0, vmax=len(arg))) 

    plt.colorbar(sm, **kwargs) 

    return sm 

def greybox(begin=45000, end=65000): 

    axvspan(begin,end, color='gainsboro')#, alpha=0.15) 

  

def flux_station(md, limit1, limit2):    ### decompose into vel and thk; need to integrate vel as well. 

    flux=[] 

    for i in range(len(md.results.TransientSolution)):  

        thk_int=integration(md, md.results.TransientSolution[i].Thickness,limit1,limit2)  

        vel=integration(md, md.results.TransientSolution[i].Vx, limit1, limit2) 

        width=np.max(md.mesh.y)-np.min(md.mesh.y) 

        flux.append((np.array(vel)/width)*(np.array(thk_int))) 

#        flux.append(thk_int) 

    return flux 

  

  

#### define functions for fjord geometry parameters 

  

def integration(md, par, limit1, limit2): 

    across=np.where(np.logical_and(md.mesh.x<limit2, md.mesh.x>=limit1)) 

    array=np.array((md.mesh.y[across], np.squeeze(par[across]))) 

    ind=np.argsort(array[0]) 

    array=array[:,ind] 

    increment=[] 

    centerval=[] 

    for i in range(0, len(array[0])-1): 

        increment.append(array[0][i+1]-array[0][i]) 

        centerval.append(((array[1][i])+array[1][i+1])/2) 

    wetted_area=np.sum(np.array(centerval)*np.array(increment)) 

    return wetted_area 

  

def along(md, parameter,limit1=15050, limit2=14950, **kwargs): 
    if hasattr(md.results, 'TransientSolution'): 
        if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 
            lens=[] 
            for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 
                lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0])) 
    #along=np.where(np.logical_and(np.logical_and(md.mesh.y<=limit1, md.mesh.y>=limit2),np.logical_or(parameter>1, parameter<-1))) 
    along=np.where(np.logical_and(md.mesh.y<=limit1, md.mesh.y>=limit2)) 
    try: 
        array=np.array((md.mesh.x[along], np.squeeze(parameter[along]))) 
    except IndexError: 
        right_index=int(np.nonzero(lens==np.intersect1d(lens,len(parameter)))[0]) 
        along=np.where(np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][1]<=limit1, md.results.TransientSolution[-1].catch_mesh[right_index][1]>=limit2)) 
        array=np.array((md.results.TransientSolution[-1].catch_mesh[right_index][0][along], np.squeeze(parameter[along]))) 
    ind=np.argsort(array[0])  
    array=array[:,ind] 
    if array[1][0]<0 and max(array[1])>-1:
        try:
            array=array[:,0:np.where(array[1]>-2)[0][1]]
        except:
            array=array[:,0:np.where(array[1]>-2)[0][0]]
    if array[1][0]>0: 
        try: 
            array=array[:,0:np.where(array[1]<2)[0][1]] 
        except: 
            print('zero value not found') 
    plt.plot(array[0], array[1],**kwargs)
    
    #plt.grid(axis='x') 

def along_vel(md, parameter,limit1=15050, limit2=14950, **kwargs):
    if hasattr(md.results, 'TransientSolution'): 
        if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 
            lens=[] 
            for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 
                lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0])) 
    along=np.where(np.logical_and(md.mesh.y<=limit1, md.mesh.y>=limit2))
    try:
        array=np.array((md.mesh.x[along], np.squeeze(parameter[along])))
    except IndexError:
        right_index=int(np.nonzero(lens==np.intersect1d(lens,len(parameter)))[0]) 
        along=np.where(np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][1]<=limit1, md.results.TransientSolution[-1].catch_mesh[right_index][1]>=limit2)) 
        array=np.array((md.results.TransientSolution[-1].catch_mesh[right_index][0][along], np.squeeze(parameter[along])))
        
    ind=np.argsort(array[0])  
    array=array[:,ind]
    #try:
    #    array=array[:,0:int(np.where(array[1]==np.max(array[1][np.where(array[0]>70000)[0]]))[0])]
    #except:
    if 'x_area' in kwargs:
        x_area=kwargs.get('x_area')
        array=array[:,0:int(np.where(array[1]==np.max(array[1][np.where(array[0]>x_area)]))[0])-2]
        del kwargs['x_area']
    else:
        array=array[:,0:int(np.where(array[1]==np.max(array[1]))[0])-2]
#    array=array[:,0:int(np.where(np.logical_and(array[1]==np.max(array[1]), array[0]>70000))[0])]
    plt.plot(array[0], array[1],**kwargs) 

    

def plotcontour(md, par, **kwargs): 
    if hasattr(md.results, 'TransientSolution'): 
        if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 
            lens=[] 
            for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 
                lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0]))
    try:
        tricontour(md.mesh.x, md.mesh.y, np.squeeze(par), **kwargs)
    except:
        right_index=int(np.nonzero(lens==np.intersect1d(lens,len(par)))[0])
        tricontour(md.results.TransientSolution[-1].catch_mesh[right_index][0], md.results.TransientSolution[-1].catch_mesh[right_index][1], np.squeeze(par), **kwargs)
 

  

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

    wetted_area=np.sum(np.array(centerval)*np.array(increment))*-1. 

    return wetted_area 

     

def along_wetted(md, wetstep): 

    wetarea=[[],[]] 

    for i in range(0, int(max(md.mesh.x)), wetstep): 

        if i<20000: 

            wetarea[0].append(wetted(md, i, i+800)) 

        else: 

            wetarea[0].append(wetted(md, i, i+wetstep)) 

    wetarea[1]=np.linspace(0,85000, len(wetarea[0])) 

    return wetarea 

  

def dalong_wetted(v_wetarea, wetstep):  ## dx is wetstep 

    dP=[[],[]] 

    dP[0]=np.array(list(reversed(np.array(deltaval(list(reversed(v_wetarea))))/wetstep))) 

    dP[1]=np.linspace(0,85000, len(dP[0])) 

    return dP 

  

def along_waterdepth(md, depthstep=350): 

    depth=[[],[]] 

    for i in range(0, int(max(md.mesh.x)), depthstep): 

        if i<20000: 

            depth[0].append(np.min(md.geometry.bed[np.where(np.logical_and(md.mesh.x>i, md.mesh.x<i+700))])) 

        else: 

            depth[0].append(np.min(md.geometry.bed[np.where(np.logical_and(md.mesh.x>i, md.mesh.x<i+depthstep))])) 

    depth[1]=np.linspace(0,85000, len(depth[0])) 

    return depth 

  

def along_width2(md, parameter, widthstep=1000, height=500): 

    width=[[],[]] 

    for i in range(0, int(max(md.mesh.x)), widthstep): 

        width[0].append(width_pos(md, parameter, i, widthstep,height)*-1) 

    width[1]=np.linspace(0,85000, len(width[0])) 

    return width 

  

def along_width(md, par, **kwargs): 

    cont=tricontour(md.mesh.x, md.mesh.y, par, **kwargs) 

    lines=cont.__dict__['allsegs'][0] 

    if len(lines)==2: 

        upper=np.transpose(lines[0])[:,:len(lines[1])] 

        lower=np.transpose(list(reversed(lines[1])))[:,:len(lines[0])] 

        width=[(np.array(upper[1])-np.array(lower[1])), upper[0]]                

    else: 

        x_coords=np.transpose(lines)[0] 

        x_coords=x_coords[:int(len(x_coords)/2)] 

        w=np.transpose(lines)[1] 

        width=[w[:int(ceil(len(w)/2))]-w[int(len(w)/2):len(w)], x_coords] 

        width[0]=width[0][:len(width[1])] 

        width[1]=width[1][:len(width[0])] 

    return width 

  

def value_pos(md, GL, wetstep): 

    value_pos=[] 

    for i in GL: 

        value_pos.append(wetted(md, i, i+wetstep)) 

    return value_pos 

  

def width_pos(md, parameter, GL, widthstep, height): 

    across=np.where(np.logical_and(md.mesh.x>GL, md.mesh.x<GL+widthstep)) 

    array=np.array((md.mesh.y[across], np.squeeze(parameter[across]))) 

    ind=np.argsort(array[0]) 

    array=array[:,ind] 

    can_vals=array[:,np.where(np.logical_and(array[1]<height+200, array[1]>height))] 

    lo_vals=[can_vals[0][0][np.nonzero(can_vals[0]<15000)[1]], can_vals[1][0][np.nonzero(can_vals[0]<15000)[1]]] 

    hi_vals=[can_vals[0][0][np.nonzero(can_vals[0]>15000)[1]], can_vals[1][0][np.nonzero(can_vals[0]>15000)[1]]] 

    width=np.mean(lo_vals[0])-np.mean(hi_vals[0]) 

    return width 

  

def slopebox(md): 

    bedslope=slope(md, md.geometry.bed) 

    beg_n=np.min(bedslope[0][np.where(bedslope[1]<0)]) 

    end_n=np.max(bedslope[0][np.where(bedslope[1]<0)])  

    beg_p=np.min(bedslope[0][np.where(bedslope[1]>0)]) 

    end_p=np.max(bedslope[0][np.where(bedslope[1]>0)]) 

    axvspan(beg_n,end_n, color='red', alpha=0.15) 

    axvspan(beg_p,end_p, color='green', alpha=0.15) 

  

def common(GL, indic,AOI,**kwargs): ## get the values for fjord characteristics at the position of the GL or mval in order to be able to correlate them with retreat characteristics 

    dic={} 

    inds_dic={} 

    for arg in kwargs: 

        close=[] 

        inds_list=[] 

        for i in GL: 

            x_coord=min(indic[arg][1], key=lambda x:abs(x-i)) 

            inds_pot_add=np.nonzero(GL==i)[0] 

            inds_add_bool=[y not in inds_list for y in inds_pot_add] 

            if i < 65000 and i>45000 and AOI==True: 

                inds_add=inds_pot_add[inds_add_bool] 

                insert_len=len(inds_add) 

                close.extend([x_coord]*insert_len) 

                inds_list.extend(inds_add) 

            elif AOI == False: 

                close.append(x_coord)

            elif isinstance(AOI, list) and i < AOI[0] and i > AOI[1]:
                
                inds_add=inds_pot_add[inds_add_bool] 

                insert_len=len(inds_add) 

                close.extend([x_coord]*insert_len) 

                inds_list.extend(inds_add)                 

            else: 

                continue 

        common=[] 

        for x in close: 

            inds=np.nonzero(indic[arg][1]==x)[0][0] 

            common.append(indic[arg][0][inds]) 

        if AOI == False: 

            inds_list=list(range(len(GL))) 

        dic[arg]=common 

        inds_dic[arg]=inds_list 

    return dic, inds_dic 

  

def fit_dP(dP): 

    dP_adj=dP 

    aoi=np.where(np.logical_and(dP[1]<65000, dP[1]>45000)) 

    dat=dP[0][aoi] 

    x=dP[1][aoi] 

    z=np.polyfit(x, dat, 4) 

    p=np.poly1d(z) 

    dP_adj=np.array(dP_adj) 

    dP_adj[0][aoi]=p(x) 

    dP_adj[1][aoi]=x 

    return dP_adj 

  

def plot3d(md, par, **kwargs): 

#    fig=plt.figure() 

    ax = plt.gca(projection='3d') 

    ax.plot_trisurf(md.mesh.x, md.mesh.y, par, linewidth=0.2, antialiased=True, **kwargs) 

  

def wa_from_name(mod): 

    newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in mod) 

    lon = np.array([float(i) for i in newstr.split()]) 

    if lon[0]<0 or lon[3]>0: 

        wa=np.max(findgateman(lon[0], lon[2],lon[3],lon[5])[0]) 

    if lon[0]>0 or lon[3]<0: 

        wa=np.min(findgateman(lon[0], lon[2],lon[3],lon[5])[0]) 

    if lon[0]==0 and lon[3]==0: 

        wa=np.mean(findgateman(lon[0], lon[2],lon[3],lon[5])[0]) 

    return wa 

  

def runmean(timeseries, window): 

    runmean=np.convolve(timeseries, np.ones((window,))/window, mode='valid') 

    return runmean 

def runsum(timeseries, window): 

    s=[] 

    for i in range(0,len(timeseries)-window): 

        s.append(np.sum(timeseries[i:i+window])) 

    return s 

def GLgate(md, par, GL, buff): 

    if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 

        lens=[] 

        for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 

            lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0])) 

    if len(par)==len(md.mesh.x): 

        GLarea=np.where(np.logical_and(np.logical_and(md.mesh.x>(GL-buff), md.mesh.x<(GL+50)), np.logical_and(md.mesh.y<15500, md.mesh.y>14500))) 

        GLthk=np.mean(par[GLarea]) 

        GLbed=np.mean(md.geometry.bed[GLarea]) 

        across=np.where(np.logical_and(np.logical_and(md.mesh.x>(GL-buff), md.mesh.x<(GL+50)), md.geometry.bed<GLbed+GLthk)) 

        array=np.array((md.mesh.y[across], np.squeeze(par[across]))) 

    elif hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 

        right_index=int(np.nonzero(lens==np.intersect1d(lens,len(par)))[0]) 

        GLarea=np.where(np.logical_and(np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][0]>(GL-buff), md.results.TransientSolution[-1].catch_mesh[right_index][0]<(GL+50)), np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][1]<17000, md.results.TransientSolution[-1].catch_mesh[right_index][1]>13000))) 

        GLthk=np.mean(par[GLarea]) 

        GLbed=np.mean(md.results.TransientSolution[-1].catch_bed[right_index][GLarea]) 

        across=np.where(np.logical_and(np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][0]>(GL-buff), md.results.TransientSolution[-1].catch_mesh[right_index][0]<(GL+50)), md.results.TransientSolution[-1].catch_bed[right_index]<GLbed+GLthk)) 

        array=np.array((md.results.TransientSolution[-1].catch_mesh[right_index][1][across], np.squeeze(par[across]))) 

    else: 

        print('No matching mesh found') 

    ind=np.argsort(array[0]) 

    array=array[:,ind] 

    increment=[] 

    centerval=[] 

    for i in range(0, len(array[0])-1): 

        increment.append(array[0][i+1]-array[0][i]) 

        centerval.append((array[1][i+1]+array[1][i])/2) 

    wetted_area=np.sum(np.array(centerval)*np.array(increment)) 

    return wetted_area, GLthk 

  

def GLgate_new(md, par, GL, buff, pattern): 
    if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 
        lens=[] 
        for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 
            lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0]))
    if len(par)==len(md.mesh.x): 
        GLarea=np.where(np.logical_and(np.logical_and(md.mesh.x>(GL-buff), md.mesh.x<(GL+50)), np.logical_and(md.mesh.y<15500, md.mesh.y>14500))) 
        GLthk=np.mean(par[GLarea]) 
        GLbed=np.mean(md.geometry.bed[GLarea])
    elif hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 
        right_index=int(np.nonzero(lens==np.intersect1d(lens,len(par)))[0])
        GLarea=np.where(np.logical_and(np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][0]>(GL-buff), md.results.TransientSolution[-1].catch_mesh[right_index][0]<(GL+50)), np.logical_and(md.results.TransientSolution[-1].catch_mesh[right_index][1]<17000, md.results.TransientSolution[-1].catch_mesh[right_index][1]>13000))) 

        GLthk=np.mean(par[GLarea]) 

        GLbed=np.mean(md.results.TransientSolution[-1].catch_bed[right_index][GLarea]) 
    per, area = wets(GL, pattern[0], pattern[1], pattern[2], pattern[3], GLthk=GLthk+GLbed)[0:2]
    
    return per, area, GLthk 

  

         

def along_GLgate(md, par,buff): 

    wetarea=[[],[]] 

    GLthk=[[],[]] 

    for i in range(0+buff, int(max(md.mesh.x)), buff): 

        if i<20000: 

            wetarea[0].append(GLgate(md, par, i, buff+200)[0]) 

            GLthk[0].append(GLgate(md, par, i, buff+200)[1]) 

        else: 

            wetarea[0].append(GLgate(md, par,i ,buff)[0]) 

            GLthk[0].append(GLgate(md, par, i, buff)[1]) 

    wetarea[1]=np.linspace(0,85000, len(wetarea[0])) 

    GLthk[1]=np.linspace(0,85000, len(wetarea[0])) 

    return wetarea, GLthk 

  

def getallpars(md, cut='all', depression='no'): 

    mod=deepcopy(md) 

    if cut=='all': 

        mod.results.TransientSolution=md.results.TransientSolution[:] 

    elif isinstance(cut, int): 

        mod.results.TransientSolution=md.results.TransientSolution[cut:] 

    elif isinstance(cut, tuple): 

        mod.results.TransientSolution=md.results.TransientSolution[cut[0]:cut[1]]        

    ssop_vals=ssop(mod, 'IceVolumeAboveFloatation',  'GroundedArea','FloatingArea', 'IcefrontMassFluxLevelset', 'GroundinglineMassFlux','IceVolume','TotalCalvingFluxLevelset', 'TotalFloatingBmb') 

    msop_vals=msop(mod, 'Vel', 'Thickness', 'Vy', 'Vx') 

    #msop_vals=msop(md, 'Vel', 'Vy', 'Vx') 

  

    all_values={**ssop_vals, **msop_vals} 

  

    all_values['GLval']=GLvalf(mod, depression=depression) 

    all_values['mval']=mvalf(mod) 

    all_values['GLvel']=GLvel(mod, all_values['GLval'])[0] 

    #all_values['GLvelx']=GLvel(md, all_values['GLval'])[1] 

    #all_values['Frvel']=GLvel(md, all_values['mval']) 

    all_values['Flux']=flux_station(mod, 20000,20300) 

  

    all_values['dGL']=deltaval(all_values['GLval']) 

    all_values['dGL_smooth']=runmean(all_values['dGL'],5) 

    all_values['dFr']=deltaval(all_values['mval']) 

  

    all_values['shelf_length']=np.array(all_values['mval'])-np.array(all_values['GLval']) 

  

    all_values_fr=['FloatingArea', 'IcefrontMassFluxLevelset', 'max Vel','max Vx', 'mval','dFr', 'shelf_length']#, 'TotalCalvingFluxLevelset'] 

  

    return all_values, all_values_fr 

def cutallpars(all_values, keyword_list, cut):
    inds=np.where(np.array(all_values['GLval'])>cut)[0]
    parlist=[np.array(all_values[x])[inds[:-1]] for x in keyword_list]
    dic={}
    for i,q in enumerate(parlist):
        dic[keyword_list[i]]=q
    return dic
    
    

    #### create fjord geometry parameters 

def ffj_chars(md, AOI, **kwargs): 

    fj_chars={}

    fj_chars['x']=[np.linspace(0,85000,851), np.linspace(0,85000,851)]

    fj_chars['P']=along_wetted(md, 600) 

    dP=dalong_wetted(fj_chars['P'][0], 600) 

    fj_chars['dP']=fit_dP(dP) 

    ddP=dalong_wetted(fj_chars['dP'][0],600) 

    fj_chars['ddP']=fit_dP(ddP) 

    fj_chars['D']=along_waterdepth(md, depthstep=350) 

    fig = plt.figure(1000) 

    fj_chars['W']=along_width(md, md.geometry.bed, levels=[500]) 

    plt.close(1000) 

    wpr=0 

    war=0 

    dwpr=0 

    dwar=0 

    if 'pattern' in kwargs: 

        pattern=kwargs.get('pattern') 

        if 'asy' in kwargs: 

            asy=kwargs.get('asy') 

        else: 

            asy='no' 

        wp=[[],np.linspace(0,85000,851)] 

        wa=[[],np.linspace(0,85000,851)] 

        if 'GL' in kwargs: 

            GL=kwargs.get('GL') 

            wpr=[[],GL] 

            war=[[], GL] 

            for q in GL:

                if AOI == False: 

                    per, area=wets(q, pattern[0], pattern[1], pattern[2], pattern[3], asy=asy)[0:2] 

                    wpr[0].append(per) 

                    war[0].append(area)

                elif AOI == True and q < 65000 and q >45000: 

                    per, area=wets(q, pattern[0], pattern[1], pattern[2], pattern[3], asy=asy)[0:2] 

                    wpr[0].append(per) 

                    war[0].append(area)

                elif isinstance(AOI, list) and q < AOI[0] and q> AOI[1]:
                    
                    per, area=wets_long(q, pattern[0], pattern[1], pattern[2], pattern[3], pattern[4], pattern[5], asy=asy)[0:2] 

                    wpr[0].append(per) 

                    war[0].append(area)

                else: 

                    continue 

            dwpr=[np.array(deltaval(wpr[0]))*-1, wpr[1]] 

            dwar=[np.array(deltaval(wpr[0]))*-1, wpr[1]]

        for q in np.linspace(0,85000,851): 
            if isinstance(AOI, list):
                per, area=wets_long(q, pattern[0], pattern[1], pattern[2], pattern[3], pattern[4], pattern[5],asy=asy)[0:2]
            else:
                per, area=wets(q, pattern[0], pattern[1], pattern[2], pattern[3], asy=asy)[0:2] 

            wp[0].append(per) 

            wa[0].append(area) 

        fj_chars['WPer']=wp 

        fj_chars['dWPer']=[np.array(deltaval(wp[0]))*-1, wp[1]] 

        fj_chars['WA']=wa 

        fj_chars['dWA']=[np.array(deltaval(wa[0]))*-1, wa[1]] 

         

    return fj_chars, wpr, war, dwpr, dwar 

  

    #### combine fjord and retreat parameters 

  

  

def glue_runs(pattern, manual_mods='no', depression='no'):
    
    if manual_mods == 'no':
        mods=glob.glob(pattern) 

        mods.sort(key=os.path.getmtime) 
    else:
        mods=manual_mods        

    all_mods={} 

    all_pars_mods={} 

    combs={} 

    for i in mods: 

        all_mods[i]=loadmodel(i) 

        all_pars_mods[i+'_all_values']=getallpars(all_mods[i], depression=depression)[0] 

    for t in all_pars_mods[mods[0]+'_all_values']: 

        for q in all_pars_mods: 

            if t in combs: 

                combs[t].extend(all_pars_mods[q][t]) 

            else: 

                combs[t]=[] 

                combs[t].extend(all_pars_mods[q][t]) 

    combs['dGL']=deltaval(combs['GLval']) 

    combs['dFr']=deltaval(combs['mval']) 

    combs['dGL_smooth']=runmean(combs['dGL'],5) 

    return combs 

            

def glue_runs_md(pattern, manual_mods='no'): 

    if manual_mods == 'no':
        mods=glob.glob(pattern) 

        mods.sort(key=os.path.getmtime) 
    else:
        mods=manual_mods
        
    all_mods={} 

    all_pars_mods={} 

    combs={} 

    catch_mesh=[] 

    catch_bed=[] 

    md=loadmodel(mods[0]) 

    for i in mods[1:]: 

        all_mods[i]=loadmodel(i) 

        if list(all_mods[i].mesh.x)!=list(md.mesh.x): 

            catch_mesh.append([all_mods[i].mesh.x, all_mods[i].mesh.y]) 

            catch_bed.append(all_mods[i].geometry.bed) 

        md.results.TransientSolution.extend(all_mods[i].results.TransientSolution) 

    md.results.TransientSolution[-1].__dict__['catch_mesh']=catch_mesh 

    md.results.TransientSolution[-1].__dict__['catch_bed']=catch_bed 

    return md 

  

def cont_all(md, front_or_GL, intervall, bar, **kwargs): 

    plotcontour(md, md.geometry.bed, levels=[0], colors='black', **kwargs)
    if 'cut' in kwargs: 
        cut=kwargs.get('cut')
    else:
        cut=0

    if hasattr(md.results.TransientSolution[-1], 'catch_mesh'): 

        print('Will have to catch mesh') 

        lens=[] 

        for q in range(len(md.results.TransientSolution[-1].catch_mesh)): 

            lens.append(len(md.results.TransientSolution[-1].catch_mesh[q][0])) 

    print('Now plotting...') 

    if 'Front' in front_or_GL: 

        colors=getcolors(len(md.results.TransientSolution)-cut, 'autumn')  

        for i in range(0,len(md.results.TransientSolution)-cut, intervall):  

            try:  

                plotcontour(md, md.results.TransientSolution[i].MaskIceLevelset, levels=[0], colors=colors[i].reshape(-1,4), **kwargs)  

            except ValueError: 

                right_index=int(np.nonzero(lens==np.intersect1d(lens,len(md.results.TransientSolution[i].MaskIceLevelset)))[0]) 

                tricontour(md.results.TransientSolution[-1].catch_mesh[right_index][0], md.results.TransientSolution[-1].catch_mesh[right_index][1], np.squeeze(md.results.TransientSolution[i].MaskIceLevelset), levels=[0], colors=colors[i].reshape(-1,4), **kwargs) 

        if bar=='yes': 

            getbar('autumn', md.results.TransientSolution, orientation='horizontal') 

    if 'GL' in front_or_GL: 

        colors=getcolors(len(md.results.TransientSolution)-cut, 'viridis')  

        for i in range(0,len(md.results.TransientSolution)-cut,intervall):  

            try:  

                plotcontour(md, md.results.TransientSolution[i].MaskGroundediceLevelset, levels=[0], colors=colors[i].reshape(-1,4), **kwargs)  

            except ValueError: 

                right_index=int(np.nonzero(lens==np.intersect1d(lens,len(md.results.TransientSolution[i].MaskIceLevelset)))[0]) 

                tricontour(md.results.TransientSolution[-1].catch_mesh[right_index][0], md.results.TransientSolution[-1].catch_mesh[right_index][1], np.squeeze(md.results.TransientSolution[i].MaskGroundediceLevelset), levels=[0], colors=colors[i].reshape(-1,4), **kwargs) 

        if bar=='yes': 

            getbar('winter', md.results.TransientSolution, orientation='horizontal')         

  

  

def get_fjord(mod, all_values, AOI, **kwargs): 
    if 'pattern' in kwargs: 
        pattern=kwargs.get('pattern') 
        if 'asy' in kwargs: 
            asy=kwargs.get('asy') 
        else: 
            asy='no' 
        if 'GL' in kwargs: 
            GL=kwargs.get('GL') 
            fj_chars, wpr, war, dwpr, dwar =ffj_chars(mod, AOI, pattern=pattern, GL=GL, asy=asy)
        else: 
            fj_chars, wpr, war, dwpr, dwar =ffj_chars(mod, AOI, pattern=pattern, asy=asy) 
    else: 
        fj_chars=ffj_chars(mod, AOI)[0]
    rfj_chars_GL, inds_dic_GL=common(all_values['GLval'], fj_chars, AOI,  **fj_chars)     
    #rfj_chars_mval, inds_dic_mval=common(all_values['mval'], fj_chars, AOI, **fj_chars) 
    rfj_chars_mval, inds_dic_mval=0,0 
    if 'pattern' and 'GL' in kwargs:
        rfj_chars_GL['WPer','WA', 'dWPer', 'dWA']= wpr, war, dwpr, dwar 
    rfj_chars_GL['Pthk']=[] 
    rfj_chars_GL['GLthk']=[] 
    rfj_chars_GL['WPthk']=[] 
    inds_dic_GL['WPthk']=[] 
    inds_dic_GL['Pthk']=[] 
    inds_dic_GL['GLthk']=[]  
    for i in range(len(all_values['GLval'])):
        if AOI == True or isinstance(AOI, list):
            if AOI==True:
                uplim=65000
                lowlim=45000
            elif isinstance(AOI, list):
                uplim=AOI[0]
                lowlim=AOI[1]
            if all_values['GLval'][i] < uplim and all_values['GLval'][i]>lowlim: 
                inds_dic_GL['Pthk'].append(i) 
                inds_dic_GL['GLthk'].append(i) 
                inds_dic_GL['WPthk'].append(i) 
                if 'pattern' in kwargs: 
                    rfj_chars_GL['WPthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[0]) 
                    rfj_chars_GL['Pthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[1]) 
                    rfj_chars_GL['GLthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[2]) 
                else: 
                    rfj_chars_GL['Pthk'].append(GLgate(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200)[0]) 
                    rfj_chars_GL['GLthk'].append(GLgate(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200)[1]) 
            else: 
                continue 
        else: 
            if 'pattern' in kwargs: 
                rfj_chars_GL['WPthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[0]) 
                rfj_chars_GL['Pthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[1]) 
                rfj_chars_GL['GLthk'].append(GLgate_new(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200, pattern)[2]) 
                inds_dic_GL['Pthk']=list(range(len(all_values['GLval']))) 
                inds_dic_GL['GLthk']=list(range(len(all_values['GLval']))) 
                inds_dic_GL['WPthk']=list(range(len(all_values['GLval']))) 
                
            else: 
                rfj_chars_GL['Pthk'].append(GLgate(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200)[0]) 
                rfj_chars_GL['GLthk'].append(GLgate(mod, mod.results.TransientSolution[i].Thickness, all_values['GLval'][i], 200)[1]) 
                inds_dic_GL['Pthk']=list(range(len(all_values['GLval']))) 
                inds_dic_GL['GLthk']=list(range(len(all_values['GLval']))) 

    return fj_chars, rfj_chars_GL, rfj_chars_mval, inds_dic_GL, inds_dic_mval 

 
