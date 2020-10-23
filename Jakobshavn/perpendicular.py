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
        continue
    xg, yg= np.array((x1,y1))+vec/2
    
    mo=(line_y/line_x)

    if mo > 0:
        m=mo-1
    if mo < 0:
        m=1-mo

    k=yg-xg*m

    profile_x=np.array([xg-10000, xg+10000])
    profile_y=profile_x*m+k

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
all_profiles=np.empty((shape(profiles)[1],402)) 
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

loc=1000 
fig, ax = plt.subplots(2)
show(topo)
plot(x,y, linewidth=5)
plot(profiles[loc][0], profiles[loc][1]) 
ax[0].plot(all_profiles[loc])

