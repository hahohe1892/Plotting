import six
import sys
sys.modules['sklearn.externals.six'] = six
import mlrose


max_vels=[]
for t in range(len(md.results.TransientSolution)):
    ground_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset>0)
    float_mask=np.squeeze(md.results.TransientSolution[t].MaskGroundediceLevelset<0)
    thk_mask=np.squeeze(md.results.TransientSolution[t].Thickness>2)
    bed_mask=np.squeeze(md.geometry.bed<0)
    no_ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset>0)
    ice_mask=np.squeeze(md.results.TransientSolution[t].MaskIceLevelset<0)
    categorized=np.zeros(md.mesh.numberofvertices)
    floating=np.logical_or(np.logical_and(float_mask, thk_mask), np.logical_and(no_ice_mask, bed_mask))
    grounded=np.logical_and(np.logical_and(np.logical_and(bed_mask, thk_mask),ground_mask), ice_mask)
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

    vels=[]
    for p in vs:
        vels.append(md.results.TransientSolution[t].Vel[p])

    max_vels.append(max(vels))

es=[]
e=open('./exp_file.exp', 'w')
e.write('## Name:Square\n## Icon:0\n# Points Count  Value\n{} 1.\n# X pos Y pos\n'.format(len(vs)))
for v in cs:
    #e.write('{} {} \n'.format(md.mesh.x[v], md.mesh.y[v]))
    e.write('{} {} \n'.format(es[v][0], es[v][1]))
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
