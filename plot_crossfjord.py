from generate_exp2 import *
fs=26
ls=22
lw=5
plot(wets(50000,0,0,0,0)[2], color='black', linewidth=lw)
plot(wets(50000,240,20000,0,0)[2], color='Gray', linewidth=lw)
plot(wets(50000,-360, 20000,0,0)[2], color='Green', linewidth=lw)
plot(wets(50000,0,0,-900,20000)[2], color='orange', linewidth=lw)
plot(wets(50000,0,0,1350, 20000)[2], color='blue', linewidth=lw)
xlim(50000,150000)
xlabel('y-coordinates [km]', fontsize=fs)
ylabel('z [m]', fontsize=fs)
xticks([50000,75000,100000,125000,150000], [0, 2.5, 5, 7.5, 10])
tick_params(axis="x", labelsize=ls)
tick_params(axis="y", labelsize=ls)
