# %%
import numpy as np
import matplotlib.pyplot as plt
import os
from flyvbjerg_petersen_std_err import fp_stderr

# %%
plt.rc('text', usetex=True)
#plt.rc('ps', usedistiller='xpdf')
plt.rc('font', family='serif')
plt.rc('axes', labelsize='xx-large')
plt.rc('xtick', labelsize='xx-large')
plt.rc('ytick', labelsize='xx-large')
plt.rc('axes', titlesize='xx-large')
plt.rc('legend', fontsize='large')

# %%
def load_data(files):
    data = [np.loadtxt(f) for f in files]
    return np.array(data)

# %%
dchi = 0.05
dbins = np.arange(-1.0, 1.0+dchi, dchi)
dbins_center = (dbins[1:] + dbins[:-1]) / 2

cbins = np.arange(0, 20, 1)

# %%
mean_flag = False

# %%
if mean_flag == True:
    prefix_hist = 'mean_trivec_hist'
    prefix_avg = 'mean_trivec_time_avg'
elif mean_flag == False:
    prefix_hist = 'trivec_hist'
    prefix_avg = 'trivec_time_avg'


# %%
# Load data
results_dir = '/mnt/hdb/Research_Data/enantiomers/results/water/10ns/'

chi_w = np.loadtxt(os.path.join(results_dir, prefix_hist+'.dat'))
norm_chi_w = np.sum(chi_w, axis=0)*dchi

chains_w = np.loadtxt(os.path.join(results_dir, 'chain_hist.dat'))
norm_chain_w = np.sum(chains_w, axis=0)

chi_w_avg = np.loadtxt(os.path.join(results_dir, prefix_avg+'.dat'))

# %%
def plot_chi_water(chi_w, bins_center, norm, results_dir, output_file, mean=False):

    fig, axes = plt.subplots(1, 1, figsize=(4.5, 4))
    ax = axes

    ax.plot(bins_center, chi_w[:, 0] / norm[0], label='$o_1$')
    ax.plot(bins_center, chi_w[:, 1] / norm[1], label='$o_2$')
    ax.plot(bins_center, chi_w[:, 2] / norm[2], label='$o_3$')
    ax.plot(bins_center, chi_w[:, 3] / norm[3], label='$o_4$')
    if mean == True:
        ax.plot(bins_center, chi_w[:, 4] / norm[4], label='all')

    ax.set_aspect(2/0.6)
    plt.ylim(0, 0.6)
    plt.xlim(-1, 1)
    if mean == True:
        plt.xlabel(r'$\overline{\chi}$')
    else:
        plt.xlabel(r'$\chi$')
    plt.ylabel('Probability')
    #plt.title('pure water')
    plt.legend(loc='lower center')
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_file))

def plot_chain_water(chains_w, bins_center, norm, results_dir, output_file):
    
    fig, axes = plt.subplots(2, 2, figsize=(6, 6))
    ax = axes.flatten()

    #plt.plot(bins_center, chains_w[:, 0] / norm[0], label='$p_1$')
    #plt.plot(bins_center, chains_w[:, 1] / norm[1], label='$p_2$')
    #plt.plot(bins_center, chains_w[:, 2] / norm[2], label='$p_3$')
    #plt.plot(bins_center, chains_w[:, 3] / norm[3], label='$p_4$')

    ax[0].bar(bins_center, chains_w[:, 0] / norm[0], label='$p_1$')
    ax[1].bar(bins_center, chains_w[:, 1] / norm[1], label='$p_2$')
    ax[2].bar(bins_center, chains_w[:, 2] / norm[2], label='$p_3$')
    ax[3].bar(bins_center, chains_w[:, 3] / norm[3], label='$p_4$')

    #plt.ylim(0, 5)
    for i, a in enumerate(ax):
        a.set_xlim(0, 20)
        if i >= 2:
            a.set_xlabel(r'distinct chains at $o_i$')
        if i == 0 or i == 2:
            a.set_ylabel('Probability')
        a.grid()
        
    ax[0].set_title(r'water, $o_1$')
    ax[1].set_title(r'water, $o_2$')
    ax[2].set_title(r'water, $o_3$')
    ax[3].set_title(r'water, $o_4$')
    #plt.legend()
    
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, output_file))


# %%
plot_chi_water(chi_w, dbins_center, norm_chi_w, results_dir, prefix_hist+'.pdf', mean=mean_flag)

# %%
plot_chain_water(chains_w, cbins, norm_chain_w, results_dir, 'chain_hist.pdf')


# %%
files = ['p1_trivec.npy', 'p2_trivec.npy', 'p3_trivec.npy', 'p4_trivec.npy']
files_avg = ['p1_mean_trivec.npy', 'p2_mean_trivec.npy', 'p3_mean_trivec.npy', 'p4_mean_trivec.npy']
p = ['p1', 'p2', 'p3', 'p4']

# %%
enantiomer = 'alanine'
rdirectory = '/mnt/hdb/Research_Data/enantiomers/results/'+enantiomer+'/r/280K/1000ns/'
sdirectory = '/mnt/hdb/Research_Data/enantiomers/results/'+enantiomer+'/s/280K/1000ns/'
bins = np.linspace(-1, 1, 201)
bins2= np.linspace(-1, 1, 51)
mid_bins = (bins[1:] + bins[:-1]) / 2

# %%
i=3
file = files[i]
file_avg = files_avg[i]
pi = p[i]
    
rdata = np.load(os.path.join(rdirectory, file))
sdata = np.load(os.path.join(sdirectory, file))

ravgdata = np.load(os.path.join(rdirectory, file_avg))
savgdata = np.load(os.path.join(sdirectory, file_avg))

#rclean = rdata[:,0][~np.isnan(rdata[:,0])]
rclean = rdata[~np.isnan(rdata)]
ravgclean = ravgdata[:,0][~np.isnan(ravgdata[:,0])]

#sclean = sdata[:,0][~np.isnan(sdata[:,0])]
sclean = sdata[~np.isnan(sdata)]
savgclean = savgdata[:,0][~np.isnan(savgdata[:,0])]

# %%
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
ax = axes.flatten()

for a in ax:
    a.grid()

rcolor = '#750014'
scolor = '#8b959e'

#rhist, _ = np.histogram(rclean, bins=bins, density=True)
#shist, _ = np.histogram(sclean, bins=bins, density=True)

ax[0].hist(rdata, edgecolor=rcolor, facecolor='None', bins=bins2, density=True,lw=1.5, label='R')
#ax[2].hist(rdata, edgecolor=rcolor, histtype='step', bins=bins, density=True,lw=1.5, label='R', log=True)

ax[0].hist(sdata, edgecolor=scolor, facecolor='None', bins=bins2, density=True,lw=1.5, label='S')
#ax[2].hist(sdata, edgecolor=scolor, histtype='step', bins=bins, density=True,lw=1.5, label='S', log=True)

ax[1].hist(ravgclean, edgecolor=rcolor, facecolor='None', bins=bins, density=True,lw=1.5, label='R')
#ax[3].hist(ravgclean, edgecolor=rcolor, facecolor='None', bins=bins, density=True,lw=1.5, label='R', log=True)

ax[1].hist(savgclean, edgecolor=scolor, facecolor='None', bins=bins, density=True,lw=1.5, label='S')
#ax[3].hist(savgclean, edgecolor=scolor, facecolor='None', bins=bins, density=True,lw=1.5, label='S', log=True)

ax[0].set_xlim(-1.0, 1.0)
ax[0].set_ylim(0.3, 0.6)
ax[1].set_xlim(-0.3, 0.3)
#ax[3].set_xlim(-0.3, 0.3)
#ax[2].set_xlim(-1.0, 1.0)

ax[0].set_xlabel(r'$\chi$')
ax[1].set_xlabel(r'$\overline{\chi}$')

ax[0].set_ylabel('Probability')
#ax[2].set_ylabel('Probability')

ax[0].legend()
#title = enantiomer+r', $o_1$, $0.0 < d < 3.5$'
#fig.suptitle(title, fontsize='xx-large')
ax[0].set_title(r'$o_4$')
ax[1].set_title(r'$o_4$')

plt.tight_layout()
plt.savefig(enantiomer+'_'+pi+'_chi_1000ns_280K.pdf')


# %%
rall = np.array([])
ravgall = np.array([])
sall = np.array([])
savgall = np.array([])
for i, file in enumerate(files):
#file = files[i]
    file_avg = files_avg[i]
    pi = p[i]
        
    rdata = np.load(os.path.join(rdirectory, file))
    sdata = np.load(os.path.join(sdirectory, file))

    ravgdata = np.load(os.path.join(rdirectory, file_avg))
    savgdata = np.load(os.path.join(sdirectory, file_avg))

    #rclean = rdata[:,0][~np.isnan(rdata[:,0])]
    rclean = rdata[~np.isnan(rdata)]
    ravgclean = ravgdata[:,0][~np.isnan(ravgdata[:,0])]

    rall = np.concatenate((rall, rclean))
    ravgall = np.concatenate((ravgall, ravgclean))
    

    #sclean = sdata[:,0][~np.isnan(sdata[:,0])]
    sclean = sdata[~np.isnan(sdata)]
    savgclean = savgdata[:,0][~np.isnan(savgdata[:,0])]

    sall = np.concatenate((sall, sclean))
    savgall = np.concatenate((savgall, savgclean))

# %%
rmean = np.nanmean(rall)
rerr = fp_stderr(rall)
print(rmean, rerr)

smean = np.nanmean(sall)
serr = fp_stderr(sall)
print(smean, serr)

rmean = np.nanmean(ravgall)
rerr = fp_stderr(ravgall)
print(rmean, rerr)

smean = np.nanmean(savgall)
serr = fp_stderr(savgall)
print(smean, serr)

# %%
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
ax = axes.flatten()

for a in ax:
    a.grid()

rcolor = '#750014'
scolor = '#8b959e'

#rhist, _ = np.histogram(rclean, bins=bins, density=True)
#shist, _ = np.histogram(sclean, bins=bins, density=True)

ax[0].hist(rall, edgecolor=rcolor, facecolor='None', bins=bins2, density=True,lw=1.5, label='R')
#ax[2].hist(rdata, edgecolor=rcolor, histtype='step', bins=bins, density=True,lw=1.5, label='R', log=True)

ax[0].hist(sall, edgecolor=scolor, facecolor='None', bins=bins2, density=True,lw=1.5, label='S')
#ax[2].hist(sdata, edgecolor=scolor, histtype='step', bins=bins, density=True,lw=1.5, label='S', log=True)

ax[1].hist(ravgall, edgecolor=rcolor, facecolor='None', bins=bins, density=True,lw=1.5, label='R')
#ax[3].hist(ravgclean, edgecolor=rcolor, facecolor='None', bins=bins, density=True,lw=1.5, label='R', log=True)

ax[1].hist(savgall, edgecolor=scolor, facecolor='None', bins=bins, density=True,lw=1.5, label='S')
#ax[3].hist(savgclean, edgecolor=scolor, facecolor='None', bins=bins, density=True,lw=1.5, label='S', log=True)

ax[0].set_xlim(-1.0, 1.0)
ax[0].set_ylim(0.3, 0.6)
ax[1].set_xlim(-0.3, 0.3)
#ax[3].set_xlim(-0.3, 0.3)
#ax[2].set_xlim(-1.0, 1.0)

ax[0].set_xlabel(r'$\chi$')
ax[1].set_xlabel(r'$\overline{\chi}$')

ax[0].set_ylabel('Probability')
#ax[2].set_ylabel('Probability')

ax[0].legend()
#title = enantiomer+r', $o_1$, $0.0 < d < 3.5$'
#fig.suptitle(title, fontsize='xx-large')
ax[0].set_title(r'all')
ax[1].set_title(r'all')

plt.tight_layout()
plt.savefig(enantiomer+'_all_chi_1000ns_280K.pdf')


# %%
