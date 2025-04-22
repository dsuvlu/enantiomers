# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
from natsort import natsorted
import os


# %%
def plot_chi_water(chi_w, bins_center, norm, results_dir, output_file, mean=False):

    fig = plt.figure(figsize=(4, 4))

    plt.plot(bins_center, chi_w[:, 0] / norm[0], label='$p_1$')
    plt.plot(bins_center, chi_w[:, 1] / norm[1], label='$p_2$')
    plt.plot(bins_center, chi_w[:, 2] / norm[2], label='$p_3$')
    plt.plot(bins_center, chi_w[:, 3] / norm[3], label='$p_4$')
    if mean == True:
        plt.plot(bins_center, chi_w[:, 4] / norm[4], label='all')

    #plt.ylim(0, 5)
    plt.xlim(-1, 1)
    if mean == True:
        plt.xlabel(r'$\overline{\chi}$')
    else:
        plt.xlabel(r'$\chi$')
    plt.ylabel('Probability')
    plt.title('water')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_file))

def plot_chain_water(chains_w, bins_center, norm, results_dir, output_file):
    
    fig, axes = plt.subplots(2, 2, figsize=(5, 5))
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
            a.set_xlabel(r'chain participation at $p_i$')
        if i == 0 or i == 2:
            a.set_ylabel('Probability')
        a.grid()
        
    ax[0].set_title('water, p1')
    ax[1].set_title('water, p2')
    ax[2].set_title('water, p3')
    ax[3].set_title('water, p4')
    #plt.legend()
    
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, output_file))

# %%
def plot_chi_enantiomers(chi_r, norm_r, chi_r_avg, chi_s, norm_s, chi_s_avg, chi_w, norm_w, 
                     bins_center, pos, results_dir, output_title, output_file, mean=False):
    
    if mean == True:
        height = 12
        n = 5
    else:
        height = 10
        n = 4
    
    fig, axes = plt.subplots(n, 2, figsize=(5, height))
    ax = axes.flatten()

    R = np.arange(0, n*2, 2)
    S = np.arange(1, n*2, 2)

    for p in range(n):
        ax[R[p]].plot(bins_center, chi_w[:, p] / norm_w[p], color='black', label='water')
        ax[S[p]].plot(bins_center, chi_w[:, p] / norm_w[p], color='black', label='water')

        d = 1
        ax[R[p]].plot(bins_center, chi_r[p,:,d] / norm_r[p,d], color = 'C1', label='$3.5 < d < 7.0$')
        ax[S[p]].plot(bins_center, chi_s[p,:,d] / norm_s[p,d], color = 'C1', label='$3.5 < d < 7.0$')

        d = 0
        ax[R[p]].plot(bins_center, chi_r[p,:,d] / norm_r[p,d], color = 'C0', label='$0.0 < d < 3.5$')
        ax[S[p]].plot(bins_center, chi_s[p,:,d] / norm_s[p,d], color = 'C0', label='$0.0 < d < 3.5$')

        ax[R[p]].text(0.3, 0.2, 'avg: {:.6f}'.format(chi_r_avg[p, 0]), transform=ax[R[p]].transAxes, fontsize=8, color='C0')
        ax[S[p]].text(0.3, 0.2, 'avg: {:.6f}'.format(chi_s_avg[p, 0]), transform=ax[S[p]].transAxes, fontsize=8, color='C0')

        ax[R[p]].text(0.3, 0.1, 'avg: {:.6f}'.format(chi_r_avg[p, 1]), transform=ax[R[p]].transAxes, fontsize=8, color='C1')
        ax[S[p]].text(0.3, 0.1, 'avg: {:.6f}'.format(chi_s_avg[p, 1]), transform=ax[S[p]].transAxes, fontsize=8, color='C1')

        r_title = 'R-' + output_title + ', ' + pos[p]
        s_title = 'S-' + output_title + ', ' + pos[p]
        ax[R[p]].set_title(r_title)
        ax[S[p]].set_title(s_title)

    bottom = S[-1:]-1
    for i, a in enumerate(ax):
        #a.set_ylim(0, 5)
        
        if mean == True:
            if i >= bottom: a.set_xlabel(r'$\overline{\chi}$')
            a.set_xlim(-0.5, 0.5)
        else:
            if i >= bottom: a.set_xlabel(r'$\chi$')
            a.set_xlim(-1.0, 1.0)

        if i in R:
            a.set_ylabel('Probability')
        a.grid()

    ax[1].legend(loc='center')
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_file))

def plot_chains_enantiomers(chains_r, norm_r, chains_s, norm_s, chains_w, norm_w, 
                     bins, pos, results_dir, output_title, output_file, mean=False):
    
    if mean == True:
        height = 12
        n = 5
    else:
        height = 10
        n = 4
    
    fig, axes = plt.subplots(n, 2, figsize=(5, height))
    ax = axes.flatten()

    R = np.arange(0, n*2, 2)
    S = np.arange(1, n*2, 2)

    for p in range(n):
        ax[R[p]].plot(bins, chains_w[:, p] / norm_w[p], color='black', label='water')
        ax[S[p]].plot(bins, chains_w[:, p] / norm_w[p], color='black', label='water')

        d = 1
        ax[R[p]].plot(bins, chains_r[p,:,d] / norm_r[p,d], color = 'C1', label='$3.5 < d < 7.0$')
        ax[S[p]].plot(bins, chains_s[p,:,d] / norm_s[p,d], color = 'C1', label='$3.5 < d < 7.0$')

        d = 0
        ax[R[p]].plot(bins, chains_r[p,:,d] / norm_r[p,d], color = 'C0', label='$0.0 < d < 3.5$')
        ax[S[p]].plot(bins, chains_s[p,:,d] / norm_s[p,d], color = 'C0', label='$0.0 < d < 3.5$')

        r_title = 'R-' + output_title + ', ' + pos[p]
        s_title = 'S-' + output_title + ', ' + pos[p]
        ax[R[p]].set_title(r_title)
        ax[S[p]].set_title(s_title)

    bottom = S[-1:]-1
    for i, a in enumerate(ax):
        #a.set_ylim(0, 5)
        
        if mean == True:
            if i >= bottom: a.set_xlabel(r'$\langle \chi \rangle$')
            a.set_xlim(-0.5, 0.5)
        else:
            if i >= bottom: a.set_xlabel(r'chain participation at $p_i$')
            a.set_xlim(0, 20)

        if i in R:
            a.set_ylabel('Probability')
        a.grid()

    ax[1].legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_file))


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
plot_chi_water(chi_w, dbins_center, norm_chi_w, results_dir, prefix_hist+'.png', mean=mean_flag)
plot_chain_water(chains_w, cbins, norm_chain_w, results_dir, 'chain_hist.png')


# %%
results_dir = ['/mnt/hdb/Research_Data/enantiomers/results/alanine/',
                '/mnt/hdb/Research_Data/enantiomers/results/2-butanol/',
                '/mnt/hdb/Research_Data/enantiomers/results/lactic_acid/']

pos = ['p1', 'p2', 'p3', 'p4', 'all']

output_title = ['alanine', '2-butanol', 'lactic_acid']

direc = [0, 1, 2]

# %%
for d in direc:
    r_files = natsorted(glob.glob(os.path.join(results_dir[d], 'r/1000ns/'+prefix_hist+'*.dat')))
    s_files = natsorted(glob.glob(os.path.join(results_dir[d], 's/1000ns/'+prefix_hist+'*.dat')))
    if mean_flag == True:
        r_files = r_files[1:] + r_files[:1]
        s_files = s_files[1:] + s_files[:1]

    chi_r = load_data(r_files)
    chi_s = load_data(s_files)

    r_files = natsorted(glob.glob(os.path.join(results_dir[d], 'r/1000ns/'+prefix_avg+'*.dat')))
    s_files = natsorted(glob.glob(os.path.join(results_dir[d], 's/1000ns/'+prefix_avg+'*.dat')))
    if mean_flag == True:
        r_files = r_files[1:] + r_files[:1]
        s_files = s_files[1:] + s_files[:1]

    chi_r_avg = load_data(r_files)
    chi_s_avg = load_data(s_files)

    r_files = natsorted(glob.glob(os.path.join(results_dir[d], 'r/1000ns/chain_hist_*.dat')))
    s_files = natsorted(glob.glob(os.path.join(results_dir[d], 's/1000ns/chain_hist_*.dat')))

    chains_r = load_data(r_files)
    chains_s = load_data(s_files)

    if mean_flag == True:
        n = 5
    else:
        n = 4

    norm_r = np.array([np.sum(chi_r[i,:,:], axis=0)*dchi for i in range(n)])
    norm_s = np.array([np.sum(chi_s[i,:,:], axis=0)*dchi for i in range(n)])

    plot_chi_enantiomers(chi_r, norm_r, chi_r_avg, chi_s, norm_s, chi_s_avg, chi_w, norm_chi_w, 
                     dbins_center, pos, results_dir[d], output_title[d], 
                     output_file=prefix_hist+'.png', mean=mean_flag)
    
    norm_r = np.array([np.sum(chains_r[i,:,:], axis=0) for i in range(4)])
    norm_s = np.array([np.sum(chains_s[i,:,:], axis=0) for i in range(4)])
    
    plot_chains_enantiomers(chains_r, norm_r, chains_s, norm_s, chains_w, norm_chain_w, 
                     cbins, pos, results_dir[d], output_title[d], 
                     output_file='chain_hist.png', mean=False)

# %%
