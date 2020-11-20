from matplotlib import colors

fs=15
lw=2


colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 256))
colors_land = plt.cm.Greys(np.linspace(0.25, 1, 256))
all_colors = np.vstack((colors_undersea, colors_land))
terrain_map = colors.LinearSegmentedColormap.from_list('terrain_map',
    all_colors)
divnorm = colors.TwoSlopeNorm(vmin=-1500., vcenter=0, vmax=500)

fig = plt.figure()
gs = fig.add_gridspec(8,9)
ax2=fig.add_subplot(gs[5:8,0:6])
ax0=fig.add_subplot(gs[0:3,0:6], sharex=ax2)
#ax1=fig.add_subplot(gs[3:5,0:6], sharex=ax2)
ax5=fig.add_subplot(gs[5:8,6:9])
ax3=fig.add_subplot(gs[0:2,6:9], sharex=ax5)
ax4=fig.add_subplot(gs[2:5,6:9], sharex=ax5)

s_col='grey'
ds_col='black'

plt.sca(ax0)
show(topo, norm=divnorm, cmap=terrain_map)
plot(x,y, linewidth=3, color=ds_col, linestyle='--')
xlim(-255000, -135000)
ylim(-2290000,-2250000)
#ax0.set_xticklabels('')

ax2.plot(x[:-1],dS_fil, color=ds_col, linewidth=3, alpha=0.15)
ax2.set_ylabel('$\it{\mathregular{dS}}$ [m\u00b2/100 m]', fontsize=fs)
ax2.set_xlabel('Easting')

ax1 = ax2.twinx()
ax1.plot(x[:-1],S_fil/1e6, color='black', linewidth=lw, linestyle='--')
ax1.set_ylabel('$\it{\mathregular{S}}$ [km\u00b2]', fontsize=fs)


plt.sca(ax3)
show(topo, cmap=terrain_map, norm=divnorm)
xlim(-177000, -155000)
ylim(-2284000,-2274000)
ax3.scatter(np.array(x)[close_inds_mask], np.array(y)[close_inds_mask], color='orange')
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
#ax3.set_xticklabels('')

ax4.plot(np.array(x)[:-1],dS_fil, linewidth=4, color=ds_col, alpha=0.15)
ax4.scatter(np.array(x)[close_inds_mask],np.array(dS_fil)[close_inds_mask], color='orange')
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
#ax4.set_xticklabels('')
ax4.set_ylabel('$\it{\mathregular{dS}}$ [m\u00b2/100 m]', fontsize=fs)
ax4.set_ylim([-26000,50000])
pyplot.locator_params(axis='y', nbins=3)

ax5.scatter(np.array(x)[close_inds_mask][dGL_criterion],(np.squeeze(GLMF)[inds_mask]/(np.array(S_fil[close_inds_mask])/1e6))[dGL_criterion], color='orange')
ax5.scatter(np.array(x)[close_inds_mask][dGL_criterion_inv],(np.squeeze(GLMF)[inds_mask]/(np.array(S_fil[close_inds_mask])/1e6))[dGL_criterion_inv], color='blue')
ax5.plot(o, p(o), linewidth=lw, linestyle=':')
ax5.yaxis.tick_right()
ax5.yaxis.set_label_position("right")
ax5.set_ylabel('$\it{\mathregular{Q_{GL}/S}}$ [km\u00b3/yr]', fontsize=fs)
ax5.set_xlabel('Easting')
pyplot.locator_params(axis='y', nbins=3)

ax21 = fig.add_subplot(gs[3:4,0:6])
cb1 = mpl.colorbar.ColorbarBase(ax21, cmap=terrain_map, norm=divnorm, orientation='horizontal')
ax21.xaxis.set_label_position("bottom")
ax21.xaxis.tick_bottom()
xlabel('Elevation [m]', fontsize=fs)
pyplot.locator_params(axis='x', nbins=5)
