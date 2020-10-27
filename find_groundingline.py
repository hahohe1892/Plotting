import six
import sys
sys.modules['sklearn.externals.six'] = six
import mlrose
from scipy.spatial import distance

max_vels=[]
min_coords=[]
for t in range(len(md.results.TransientSolution)):
    geofence_y=np.logical_or(md.mesh.y>-2260400, md.mesh.y<-2284900)
    geofence_x=np.logical_or(md.mesh.x<-251300, md.mesh.x>-123300)
    ground_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset>0)
    float_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset<0)
    thk_mask=np.squeeze(md.results.TransientSolution[t].Thickness>2)
    bed_mask=np.squeeze(md.geometry.bed<0)
    no_ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset>0)
    ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset<0)
    vel_mask=np.squeeze(md.results.TransientSolution[t].Vel<2000)
    categorized=np.zeros(md.mesh.numberofvertices)
    floating=np.logical_or(np.logical_and(float_mask, thk_mask), np.logical_and(no_ice_mask, bed_mask))
    grounded=np.logical_and(np.logical_and(np.logical_and(bed_mask, thk_mask),ground_mask), ice_mask)
    categorized[floating]=-1
    categorized[grounded]=1
    categorized[geofence_x]=0
    categorized[geofence_y]=0
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

    vels=[]
    min_xs=[]
    for p in vs:
        vels.append(md.results.TransientSolution[t].Vel[p])
        min_xs.append(md.mesh.x[p])
    max_vels.append(max(vels))
    if len(min_xs)%2==0:
        min_xs=np.sort(min_xs)[1:]
    min_ind=np.where(md.mesh.x==np.median(min_xs))
    min_coords.append([float(md.mesh.x[min_ind]), float(md.mesh.y[min_ind])])


f=open('./Jakobshavn/Flowline_points.txt', 'r')

x=[]
y=[]
for line in f:
    c=line.split(';')
    xs=c[3].split(',')
    ys=c[2].split(',')
    try:
        int_x=int(xs[0])
    except:
        continue
    dez_x=float('0.'+xs[1])
    int_y=int(ys[0])
    dez_y=float('0.'+ys[1])
    if int_x==0:
        continue
    else:
        x.append(int_x+dez_x)
        y.append(int_y+dez_y)

f.close()

close_inds=[]
for k,r in enumerate(min_coords):
    dis=[]
    for s in range(len(x)-1):
        dis.append(distance.euclidean(r, (x[s], y[s])))
    close=int(np.where(dis==np.min(dis))[0])
    close_inds.append(close)

    
    


es=[]
e=open('./exp_file.exp', 'w')
e.write('## Name:Square\n## Icon:0\n# Points Count  Value\n{} 1.\n# X pos Y pos\n'.format(len(vs)))
for v in vs:
    e.write('{} {} \n'.format(md.mesh.x[v], md.mesh.y[v]))
    #e.write('{} {} \n'.format(es[v][0], es[v][1]))
#    es.append((md.mesh.x[v], md.mesh.y[v]))
e.close()
                                                                                        
        

plotmodel(md, 'data', categorized, 'expdisp', 'exp_file.exp', 'data', md.results.TransientSolution[t].Vel, 'data', md.geometry.bed, 'data', md.results.TransientSolution[t].MaskGroundediceLevelset)

for i,element in enumerate(md.mesh.elements-1):
    adj=md.mesh.elementconnectivity[i]-1
    for q in adj:
        if md.results.TransientSolution[0].MaskGroundediceLevelset[q]<0:
            print('grounded')
        else:
            print('floating')

#md.stressbalance.requested_outputs=['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', 'MassFlux1']
#md.outputdefinition.definitions = [massfluxatgate('name', 'MassFlux1', 'profilename', './Jakobshavn/MF.exp', 'definitionstring', 'Outputdefinition1')]

# this can make a node based (len(f)==md.mesh.numberofnodes) field to an element (len(f)==md.mesh.numberofelements) based on
t=np.zeros(md.mesh.numberofelements)
for r in range(len(categorized)):
    for u,z in enumerate(md.mesh.elements-1):
        for q in z:
            if q==r:
                t[u]=categorized[r]
