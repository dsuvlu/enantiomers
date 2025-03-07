# %%
import numpy as np
import matplotlib.pyplot as plt
from flyvbjerg_petersen_std_err import fp_stderr
import os
from scipy.stats import skew, kurtosis
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp
from scipy.stats import mannwhitneyu
import scipy.stats as stats

# %%
files = ['p1_trivec.npy', 'p2_trivec.npy', 'p3_trivec.npy', 'p4_trivec.npy']
p = ['p1', 'p2', 'p3', 'p4']


# %%
rdirectory = '/home/dsuvlu/git/enantiomers/results/alanine/r/1000ns/'
sdirectory = '/home/dsuvlu/git/enantiomers/results/alanine/s/1000ns/'
bins = np.linspace(-1, 1, 201)

# %%
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
ax = axes.flatten()
for a, file in enumerate(files):
    
    rdata = np.load(os.path.join(rdirectory, file))
    sdata = np.load(os.path.join(sdirectory, file))

    rclean = rdata[:,0][~np.isnan(rdata[:,0])]
    rmean = np.nanmean(rclean)
    rerr = fp_stderr(rclean)
    print(rmean, rerr)

    sclean = sdata[:,0][~np.isnan(sdata[:,0])]
    smean = np.nanmean(sclean)
    serr = fp_stderr(sclean)
    print(smean, serr)

    ks_stat, p_value = ks_2samp(sclean, rclean)
    print("K-S test statistic:", ks_stat, "p-value:", p_value)

    ad_stat, ad_val, ad_pval = anderson_ksamp([sclean, rclean])
    print("Anderson-Darling test statistic:", ad_stat, "p-value:", ad_pval)

    u_stat, p_value = mannwhitneyu(sclean, rclean)
    print("Mann-Whitney U test statistic:", u_stat, "p-value:", p_value)

    plt.rcParams.update({'font.size': 10})

    ax[a].grid(zorder=-1000)

    ax[a].hist(rdata[:,0], edgecolor='C0', facecolor='None', bins=bins, density=True,lw=1.5, label='R', log=True)
    ax[a].hist(sdata[:,0], edgecolor='C1', facecolor='None', bins=bins, density=True,lw=1.5, label='S', log=True)

    ax[a].text(0.01, 0.9, f"Mean: {rmean:.2e} ± {rerr:.2e}", transform=ax[a].transAxes, color='C0')
    ax[a].text(0.01, 0.85, f"Mean: {smean:.2e} ± {serr:.2e}", transform=ax[a].transAxes, color='C1')

    ax[a].text(0.01, 0.8, f"Std: {np.std(rclean):.2e}", transform=ax[a].transAxes, color='C0')
    ax[a].text(0.01, 0.75, f"Std: {np.std(sclean):.2e}", transform=ax[a].transAxes, color='C1')

    ax[a].text(0.01, 0.7, f"Skewness: {skew(rclean):.2e}", transform=ax[a].transAxes, color='C0')
    ax[a].text(0.01, 0.65, f"Skewness: {skew(sclean):.2e}", transform=ax[a].transAxes, color='C1')

    ax[a].text(0.01, 0.6, f"Kurtosis: {kurtosis(rclean):.2e}", transform=ax[a].transAxes, color='C0')
    ax[a].text(0.01, 0.55, f"Kurtosis: {kurtosis(sclean):.2e}", transform=ax[a].transAxes, color='C1')

    ax[a].text(0.65, 0.85, f"K-S test statistic: {ks_stat:.2e}", transform=ax[a].transAxes)
    ax[a].text(0.65, 0.80, f"P-value: {p_value:.2e}", transform=ax[a].transAxes)

    ax[a].text(0.65, 0.75, f"A-D test statistic: {ad_stat:.2e}", transform=ax[a].transAxes)
    ax[a].text(0.65, 0.70, f"P-value: {ad_pval:.2e}", transform=ax[a].transAxes)

    ax[a].text(0.65, 0.65, f"M-W U test statistic: {u_stat:.2e}", transform=ax[a].transAxes)
    ax[a].text(0.65, 0.60, f"P-value: {p_value:.2e}", transform=ax[a].transAxes)

    ax[a].set_xlim(-1.0, 1.0)
    ax[a].set_xlabel(r'$\overline{\chi}$')
    ax[a].set_ylabel('Probability')
    ax[a].legend()
    title = 'alanine, ' + p[a] + ', $0.0 < d < 3.5$'
    ax[a].set_title(title)

plt.tight_layout()
plt.savefig('alanine.png')



# %%
stats.describe(sclean)

# %%
stats.describe(rclean)

# %%
# %%
for a, file in enumerate(files):
    
    rdata = np.load(os.path.join(rdirectory, file))
    sdata = np.load(os.path.join(sdirectory, file))

    rclean = rdata[:,0][~np.isnan(rdata[:,0])]
    rmean = np.nanmean(rclean)
    rerr = fp_stderr(rclean)
    print(rmean, rerr)

    sclean = sdata[:,0][~np.isnan(sdata[:,0])]
    smean = np.nanmean(sclean)
    serr = fp_stderr(sclean)
    print(smean, serr)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))

    plt.subplot(1, 2, 1)
    stats.probplot(sclean, dist="norm", plot=plt)
    plt.title(r's enantiomer, $0.0 < d < 3.5$')

    plt.subplot(1, 2, 2)
    stats.probplot(rclean, dist="norm", plot=plt)
    plt.title(r'r enantiomer, $0.0 < d')

    plt.tight_layout()
    plt.show()

# %%
# For two-sample Q-Q plot, you might compare quantiles directly:
quantiles1 = np.percentile(sclean, np.linspace(0, 100, 100))
quantiles2 = np.percentile(rclean, np.linspace(0, 100, 100))
plt.plot(quantiles1, quantiles2, 'o')
plt.xlabel("Data 1 Quantiles")
plt.ylabel("Data 2 Quantiles")
plt.title("Two-sample Q-Q Plot")
plt.show()



binc = (bins[1:] + bins[:-1])/2
rhist, _ = np.histogram(rclean, bins=bins, density=True)
shist, _ = np.histogram(sclean, bins=bins, density=True)

# %%
print("Data s mean:", np.mean(sclean), "Data r mean:", np.mean(rclean))
print("Data s variance:", np.var(sclean), "Data r variance:", np.var(rclean))
print("Data s skewness:", skew(sclean), "Data r skewness:", skew(rclean))
print("Data s kurtosis:", kurtosis(sclean), "Data r kurtosis:", kurtosis(rclean))
