# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
from natsort import natsorted
import os


# %%
def plot_water(w_chi_avg, bins_center, norm, results_dir, output_file):
    fig = plt.figure(figsize=(4, 4))

    plt.plot(bins_center, w_chi_avg[:, 0] / norm[0], label='$p_1$')
    plt.plot(bins_center, w_chi_avg[:, 1] / norm[1], label='$p_2$')
    plt.plot(bins_center, w_chi_avg[:, 2] / norm[2], label='$p_3$')
    plt.plot(bins_center, w_chi_avg[:, 3] / norm[3], label='$p_4$')
    plt.plot(bins_center, w_chi_avg[:, 4] / norm[4], label='all')

    plt.ylim(0, 5)
    plt.xlim(-1, 1)
    plt.xlabel(r'$\chi$')
    plt.ylabel('Probability')
    plt.title('water')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_file))

# %%
def plot_enantiomers(er_chi_avg, norm_er, es_chi_avg, norm_es, 
                    w_chi_avg, norm_w, bins_center, pos, results_dir, output):
    fig, axes = plt.subplots(5, 2, figsize=(5, 12))
    ax = axes.flatten()

    R = np.arange(0, 10, 2)
    S = np.arange(1, 10, 2)

    for p in range(5):
        ax[R[p]].plot(bins_center, w_chi_avg[:, p] / norm_w[p], color='black', label='water')
        ax[S[p]].plot(bins_center, w_chi_avg[:, p] / norm_w[p], color='black', label='water')

        d = 1
        ax[R[p]].plot(bins_center, er_chi_avg[p,:,d] / norm_er[p,d], color = 'C1', label='$3.0 < d < 6.0$')
        ax[S[p]].plot(bins_center, es_chi_avg[p,:,d] / norm_es[p,d], color = 'C1', label='$3.0 < d < 6.0$')

        d = 0
        ax[R[p]].plot(bins_center, er_chi_avg[p,:,d] / norm_er[p,d], color = 'C0', label='$0.0 < d < 3.0$')
        ax[S[p]].plot(bins_center, es_chi_avg[p,:,d] / norm_es[p,d], color = 'C0', label='$0.0 < d < 3.0$')

        r_title = 'R-' + output + ', ' + pos[p]
        s_title = 'S-' + output + ', ' + pos[p]
        ax[R[p]].set_title(r_title)
        ax[S[p]].set_title(s_title)


    for i, a in enumerate(ax):
        #a.set_ylim(0, 5)
        a.set_xlim(-0.5, 0.5)
        if i >= 8:
            a.set_xlabel(r'$\chi$')
        if i in R:
            a.set_ylabel('Probability')
        a.grid()

    ax[1].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, 'hist_mean.png'))

# %%
def load_data(files):
    data = [np.loadtxt(f) for f in files]
    return np.array(data)

# %%
dchi = 0.05
bins = np.arange(-1.0, 1.0+dchi, dchi)
bins_center = (bins[1:] + bins[:-1]) / 2

# %%
# Load data
results_dir = '/home/dsuvlu/git/enantiomers/results/water/'

w_chi_avg = np.loadtxt(os.path.join(results_dir, 'hist_mean.dat'))
norm_w = np.sum(w_chi_avg, axis=0)*dchi

# %%
plot_water(w_chi_avg, bins_center, norm_w, results_dir, 'hist_mean.png')


# %%
results_dir = ['/home/dsuvlu/git/enantiomers/results/alanine/',
                '/home/dsuvlu/git/enantiomers/results/2-butanol/',
                '/home/dsuvlu/git/enantiomers/results/lactic_acid/']

pos = ['p1', 'p2', 'p3', 'p4', 'all']

output = ['alanine', '2-butanol', 'lactic_acid']

# %%
r_files = natsorted(glob.glob(os.path.join(results_dir[0], 'r/hist_mean_*.dat')))
r_files = r_files[1:] + r_files[:1]
s_files = natsorted(glob.glob(os.path.join(results_dir[0], 's/hist_mean_*.dat')))
s_files = s_files[1:] + s_files[:1]

# %%
er_chi_avg = load_data(r_files)
es_chi_avg = load_data(s_files)

# %%
norm_er = np.array([np.sum(er_chi_avg[i,:,:], axis=0)*dchi for i in range(5)])
norm_es = np.array([np.sum(es_chi_avg[i,:,:], axis=0)*dchi for i in range(5)])

# %%
for d in range(3):
    plot_enantiomers(er_chi_avg, norm_er, es_chi_avg, norm_es, 
                 w_chi_avg, norm_w, bins_center, pos, results_dir[d], output[d])

# %%
