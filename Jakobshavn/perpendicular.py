import numpy as np
import rasterio
from rasterio.plot import show

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
for i in range(shape(profiles)[1]):
    f.write('{} {}\n'.format(profiles[0][i], profiles[1][i]))
f.close()

for q in profiles:
    plot(q[0], q[1])



f=open('./Jakobshavn/all_profiles.txt', 'r')
all_profiles=np.empty((shape(profiles)[0],1000)) 
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

loc=146
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

S_fil=scipy.signal.savgol_filter(S, window_length=31, polyorder=5)
dS_fil=[S_fil[s+1]-S_fil[s] for s in range(len(S_fil)-1)]
dS_fil.append(dS_fil[-1])


fig, ax = plt.subplots(3)
ax[1].plot(S)
ax[1].plot(S_fil)
ax[1].plot(np.array(list(range(len(S_fil))))[close_inds],np.array(S_fil)[close_inds], 'rd')
xlim(0, len(S_fil))
ax[0].plot(dS_fil)
ax[0].plot(np.array(list(range(len(dS_fil))))[close_inds],np.array(dS_fil)[close_inds], 'rd')
xlim(0, len(S_fil))
ax[2].plot(np.array(list(range(len(dS_fil))))[close_inds], (GLMF/np.array(S_fil[close_inds]))[0], 'rd')


plot(np.array(dS_fil)[close_inds_mask],(GLMF/np.array(S_fil[close_inds_mask]))[0], 'rd')

close_inds_mask=np.array(close_inds)[np.array(close_inds)>800]
