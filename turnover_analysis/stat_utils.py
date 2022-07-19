import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kstest



def compute_corrcoef_std(max_sample_size, R=5000):
    """
    Computing standard deviations of null correlation models at given sample size
    """
    ns = np.arange(2,max_sample_size)
    stds = dict()
    for n in ns:
        rhos = []
        for _ in range(R):
            r,p = spearmanr(np.random.rand(n), np.random.rand(n))
            rhos.append(r)
        stds[n] = np.std(rhos)
    return stds


def compute_rescaled_corr(corr_samples1, corr_samples2, is_self, stds, temp_shift=0):
    k = list(corr_samples1.keys())[0]
    n_m1, n_m2 = len(corr_samples1[k]), len(corr_samples2[k])
    rhos = np.zeros((len(corr_samples1), n_m1, n_m2))
    pvals = np.zeros((len(corr_samples1), n_m1, n_m2))
    
    # Computing the re-scaled correlation coefficients
    n_sample = 0
    for sample in corr_samples1:
        for i1, m1 in enumerate(corr_samples1[sample]):
            for i2, m2 in enumerate(corr_samples2[sample]):
                if not is_self or i1>i2:
                    r,p = spearmanr(m1, m2)
                    if temp_shift >= 0:
                        r,p = spearmanr(m1[:len(m1)-temp_shift], m2[temp_shift:])
                    else:
                        r,p = spearmanr(m1[abs(temp_shift):], m2[:len(m1)-abs(temp_shift)])
                    rhos[n_sample][i1][i2] = r/stds[len(m1)]
                    pvals[n_sample][i1][i2] = p/stds[len(m1)]
        n_sample += 1
        
    return rhos, pvals
        
        
def compute_KS_stat(rho_matrix, is_self):
    n_samp, n_m1, n_m2 = rho_matrix.shape
    
    kss = np.zeros((n_m1, n_m2))
    kspval = np.zeros((n_m1, n_m2))
    av_rho = np.zeros((n_m1, n_m2))
    for i1 in range(n_m1):
        for i2 in range(n_m2):
            if not is_self or i1>i2:
                ks, ksp = kstest(rho_matrix[:,i1,i2], 'norm')
                kss[i1, i2] = ks
                kspval[i1, i2] = ksp
                av_rho[i1, i2] = rho_matrix[:,i1, i2].mean()
    return kss, kspval, av_rho


def plot_correlation_heatmap(ax1, ax2, kspval, av_rho, m1_labels, m2_labels, is_self):
    if is_self:
        kspval_frame = pd.DataFrame(np.log10(kspval[1:,:-1]), index=m1_labels[1:], columns=m2_labels[:-1])
        avrho_frame = pd.DataFrame(av_rho[1:,:-1], index=m1_labels[1:], columns=m2_labels[:-1])
        mask = np.triu(np.ones_like(kspval_frame, dtype=bool),1)
    else:
        kspval_frame = pd.DataFrame(np.log10(kspval), index=m1_labels, columns=m2_labels)
        avrho_frame = pd.DataFrame(av_rho, index=m1_labels, columns=m2_labels)
        mask = np.zeros_like(kspval_frame, dtype=bool)

    cmap1 = sns.color_palette("coolwarm_r", as_cmap=True)
    cmap2 = sns.diverging_palette(150, 295, s=50, as_cmap=True)

    sns.heatmap(kspval_frame, ax=ax1, mask=mask, cmap=cmap1, vmax=0, vmin=-4, center=-2,
                linewidths=.5, annot=True, cbar_kws={'shrink':.8, 'label':'Log10 K-S p-value'})

    sns.heatmap(avrho_frame, ax=ax2, mask=mask, cmap=cmap2, vmax=1, vmin=-1, center=0, yticklabels=[],
                linewidths=.5, annot=True, cbar_kws={'shrink':.8, 'label':'Average rescaled corr'})

    plt.tight_layout()
    return ax1, ax2


def _convert_lognormal_params(mu, cov):
    """Helper function to convert lognormal means and covariances to normal
    means and covariances.
    """
    sigma = np.diag(cov)
    norm_sigma = np.log(1 + (sigma / (mu ** 2)))
    norm_mu = np.log(mu) - 0.5 * norm_sigma
    i, j = np.indices(cov.shape)
    norm_cov = np.log(1 + (cov[i, j] / np.exp(norm_mu[i] + norm_mu[j] + 0.5 * (norm_sigma[i] + norm_sigma[j]))))
    return norm_mu, norm_cov


def multivariate_lognormal(mu, cov, size=1, clip=True, spectrum_min=1e-12):
    """
    Draw random samples from a multivariate lognormal distribution.
    The multivariate lognormal is a generalization of the one-dimensional
    lognormal distribution to higher dimensions.  Such a distribution
    is specified by its mean and covariance matrix.  These parameters
    are analogous to the mean (average or "center") and
    variance (standard deviation, or "width," squared) of the
    one-dimensional lognormal distribution.
    Parameters
    ----------
    mean : 1-D array_like, of length N
        Mean of the N-dimensional distribution.
    cov : 2-D array_like, of shape (N, N)
        Covariance matrix of the distribution. It must be symmetric and
        positive-semidefinite for proper sampling.
    size : int or tuple of ints, optional
        Given a shape of, for example, ``(m,n,k)``, ``m*n*k`` samples are
        generated, and packed in an `m`-by-`n`-by-`k` arrangement.  Because
        each sample is `N`-dimensional, the output shape is ``(m,n,k,N)``.
        If no shape is specified, a single (`N`-D) sample is returned.
    clip: boolean, optional
        Clip the spectrum of the normal covariance matrix in order to ensure
        that the conversion from lognormal to normal does not violate the
        non-negative semi-definiteness of the covariance matrix.
    spectrum_min: float, optional
        The clipping value for the eigenvalues of the covariance matrix for the
        associated normal distribution.
    Returns
    -------
    out : ndarray
        The drawn samples, of shape *size*, if that was provided.  If not,
        the shape is ``(N,)``.
        In other words, each entry ``out[i,j,...,:]`` is an N-dimensional
        value drawn from the distribution.
    Notes
    -----
    The samples for a multivariate lognormal distribution are drawn by computing
    the mean vector and covariance matrix for the associated normal distribution.
    The normal samples are then exponentiated in order to obtain the multivariate
    lognormal samples.
    Examples
    --------
    >>> import numpy as np
    >>> N = 10
    >>> mu = np.arange(1, N+1)
    >>> rho = 0.5
    >>> sigmas = np.arange(1, N+1)
    >>> cov = np.diag(sigmas)
    >>> multivariate_lognormal(mu, cov, size=100)
    """
    norm_mu, norm_cov = _convert_lognormal_params(mu, cov)
    if clip is True:
        eigen_vals, eigen_vec = np.linalg.eig(norm_cov)
        eigen_vals[eigen_vals <= 0] = spectrum_min
        norm_cov = np.dot(eigen_vec, eigen_vals[:, np.newaxis] * eigen_vec.T)
    norm_samples = np.random.multivariate_normal(norm_mu, norm_cov, size=size)
    return np.exp(norm_samples)