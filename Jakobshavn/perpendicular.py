import numpy as np
import rasterio
from rasterio.plot import show
import scipy.signal
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

profiles=[]
profiles_bash=[]
for i in range(len(x)-1):
    x1=x[i]
    y1=y[i]
    x2=x[i+1]
    y2=y[i+1]

    line_x=x2-x1
    line_y=y2-y1
    
    vec=np.array((line_x, line_y))

    if 0 in vec:
        profile_x=np.array([x[i], x[i]])
        profile_y=np.array([y[i]+10000, y[i]-10000])
        profiles.append([profile_x, profile_y])
        profiles_bash.append(profile_x[0])
        profiles_bash.append(profile_x[1])
        profiles_bash.append(profile_y[0])
        profiles_bash.append(profile_y[1])
        continue
    
    xg, yg= np.array((x1,y1))+vec/2
    
    mo=(line_y/line_x)

    m=-1/mo

    k=y[i]-x[i]*m

    profile_x=np.array([x[i]-10000/m, x[i]+10000/m])
    profile_y=profile_x*m+k
    #print(np.sqrt((profile_x[1]-profile_x[0])**2+(profile_y[1]-profile_y[0])**2))
    profiles.append([profile_x, profile_y])

    profiles_bash.append(profile_x[0])
    profiles_bash.append(profile_x[1])
    profiles_bash.append(profile_y[0])
    profiles_bash.append(profile_y[1])
    #profiles.append((tuple(profile_x), tuple(profile_y)))


f=open('Profiles_bash.txt','w')
for element in profiles_bash:
    f.write('{}\n'.format(element))
f.close()


f=open('Profiles.txt','w')
for i in range(np.shape(profiles)[1]):
    f.write('{} {}\n'.format(profiles[0][i], profiles[1][i]))
f.close()

for q in profiles:
    plot(q[0], q[1])



f=open('./Jakobshavn/all_profiles.txt', 'r')
all_profiles=np.empty((np.shape(profiles)[0],1000)) 
all_profiles[:]=np.nan
prof_count=0
i=0
max_ind=0
for line in f:
    if i==0:
        ind_count=0
    l=line.split(' ')[-1].strip()
    if l == 'new_profile':
        prof_count+=1
        ind_count=0
        continue
    else:
        all_profiles[prof_count,ind_count]=int(l)
        ind_count+=1
    i+=1
    if ind_count> max_ind:
        max_ind=ind_count
f.close()

    



fp = r'./Jakobshavn/BedMachine_clip.tif'
topo = rasterio.open(fp)

loc=810
fig, ax = plt.subplots(2)
show(topo)
plot(x,y, linewidth=5)
plot(profiles[loc][0], profiles[loc][1]) 
ax[0].plot(all_profiles[loc])

S=[]
for a in all_profiles:
    n, m = 1,1
    na_len=np.count_nonzero(~np.isnan(a))
    ha_len=int(na_len/2)
    w=a[ha_len]*150
    while a[ha_len+m]<0 and (ha_len+m)*150<=20000:
        w+=a[ha_len+m]*150
        m+=1
    while a[ha_len-n]<0 and (ha_len-n)*150>=0:
        w+=a[ha_len-n]*150
        n+=1
    S.append(w*-1)


dS=[S[s+1]-S[s] for s in range(len(S)-1)]
dS.append(dS[-1])

S_fil=scipy.signal.savgol_filter(S, window_length=51, polyorder=3)
dS_fil=[S_fil[s+1]-S_fil[s] for s in range(len(S_fil)-1)]
dS_fil.append(dS_fil[-1])
dS_filfil=scipy.signal.savgol_filter(dS_fil, window_length=31, polyorder=5)

fig, ax = plt.subplots(4)#, sharex=True)
show(topo)
ax[3].plot(np.array(x)[close_inds], np.array(y)[close_inds], 'rd')
ax[0].plot(S)
ax[0].plot(S_fil)
ax[0].plot(np.array(list(range(len(S_fil))))[close_inds],np.array(S_fil)[close_inds], 'rd')
ax[1].plot(dS_fil)
ax[1].plot(np.array(list(range(len(dS_fil))))[close_inds],np.array(dS_fil)[close_inds], 'rd')
ax[2].plot(np.array(list(range(len(dS_fil))))[close_inds], (np.squeeze(GLMF)[:len(close_inds)]/np.array(S_fil[close_inds])), 'rd')

plot(np.array(dS_fil)[close_inds],(np.squeeze(GLMF)[:len(close_inds)]/np.array(S_fil[close_inds])), 'rd')

#inds_mask=np.where(np.logical_and(np.array(close_inds)>850, np.array(close_inds)<1050))
inds_mask=np.where(np.logical_and(np.array(x)[close_inds]>-177000, np.array(x)[close_inds]<-155000))
close_inds_mask=np.array(close_inds)[inds_mask]
plot(np.array(dS_fil)[close_inds_mask],(np.squeeze(GLMF)[inds_mask]/np.array(S_fil[close_inds_mask])), 'rd')

lw=2
fig = plt.figure()
gs = fig.add_gridspec(6,9)
ax0=fig.add_subplot(gs[2:4,0:6])
ax1=fig.add_subplot(gs[2:4,6:9])
ax2=fig.add_subplot(gs[4:6,6:9])
ax3=fig.add_subplot(gs[0:2,0:6])
ax4=fig.add_subplot(gs[0:2,6:9])
ax5=fig.add_subplot(gs[4:6,0:6])
plt.sca(ax3)
show(topo, norm=divnorm, cmap=terrain_map)
plot(x,y, linewidth=3, color='lightgreen', linestyle='--')
xlim(-255000, -135000)
ylim(-2290000,-2250000)
plt.sca(ax4)
show(topo, cmap=terrain_map, norm=divnorm)
xlim(-177000, -155000)
ylim(-2284000,-2274000)
ax4.scatter(np.array(x)[close_inds_mask], np.array(y)[close_inds_mask], color='orange')
#ax0.plot(np.array(x)[np.min(close_inds_mask):np.max(close_inds_mask)],S[np.min(close_inds_mask):np.max(close_inds_mask)], alpha=0.3)
#ax0.plot(np.array(x)[np.min(close_inds_mask):np.max(close_inds_mask)],S_fil[np.min(close_inds_mask):np.max(close_inds_mask)])
#ax0.plot(np.array(x)[close_inds_mask],np.array(S_fil)[close_inds_mask], 'go')
ax0.plot(x[:-1],S, color='lightgreen', linewidth=lw)
ax0.plot(x[:-1],S_fil, color='darkgreen', linewidth=lw)
ax1.scatter(np.array(x)[close_inds_mask],np.array(dS_fil)[close_inds_mask], color='orange')
ax1.plot(np.array(x)[np.min(close_inds_mask):np.max(close_inds_mask)],dS_fil[np.min(close_inds_mask):np.max(close_inds_mask)], linewidth=lw)
ax2.scatter(np.array(x)[close_inds_mask],(np.squeeze(GLMF)[inds_mask]/np.array(S_fil[close_inds_mask])), color='orange')
ax5.plot(x[:-1],dS_fil, color='lightgreen', linewidth=lw)

ax1.set_xticklabels('')
ax3.set_xticklabels('')
plt.sca(ax0)
ylabel('$\it{\mathregular{S}}$ [m\u00b2]', fontsize=fs)
xlabel('Easting')
pyplot.locator_params(axis='x', nbins=3)
plt.sca(ax1)
ylabel('$\it{\mathregular{dS}}$ [m\u00b2/100 m]', fontsize=fs)
pyplot.locator_params(axis='y', nbins=3)
plt.sca(ax2)
ylabel('$\it{\mathregular{Q_{GL}/S}}$ [km\u00b3/yr]', fontsize=fs)
xlabel('Easting')
pyplot.locator_params(axis='y', nbins=3)
pyplot.locator_params(axis='x', nbins=4)
plt.sca(ax3)
pyplot.locator_params(axis='x', nbins=3)
pyplot.locator_params(axis='y', nbins=3)
plt.sca(ax4)
pyplot.locator_params(axis='x', nbins=3)
pyplot.locator_params(axis='y', nbins=3)
ax1.yaxis.tick_right()
ax2.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
ax2.yaxis.set_label_position("right")
ax0.xaxis.set_label_position("top")
ax0.xaxis.tick_top()

ax21 = fig.add_subplot(gs[4:6,1:2])
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=terrain_map, norm=divnorm, orientation='vertical')
pyplot.locator_params(axis='y', nbins=4)
ax21.yaxis.set_label_position("left")
ax21.yaxis.tick_left()

fit=np.polyfit(np.array(x)[close_inds_mask],(np.squeeze(GLMF)[inds_mask]/(np.array(S_fil[close_inds_mask])/1e6)),5)
p=np.poly1d(fit)
o=np.linspace(np.min(np.array(x)[close_inds_mask]), np.max(np.array(x)[close_inds_mask]))
ax2.plot(o, p(o), linewidth=lw)

dGL_criterion=[True]
dGL_criterion.extend([close_inds_mask[r]!=close_inds_mask[r+1] for r in range(len(close_inds_mask)-1)])
dGL_criterion_inv=[not l for l in dGL_criterion]

colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 256))
colors_land = plt.cm.Greys(np.linspace(0.25, 1, 256))
all_colors = np.vstack((colors_undersea, colors_land))
terrain_map = colors.LinearSegmentedColormap.from_list('terrain_map',
    all_colors)
divnorm = colors.TwoSlopeNorm(vmin=-1500., vcenter=0, vmax=500)

