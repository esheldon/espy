def extract_mcmc_stats(data, weights=None):
    """
    Extract mean and covariance from mcmc trials

    parameters
    ----------
    trials: array
        (nsteps, npars) array.
    weights: optional
        Optional weights for each trial

    Returns
    -------
    means, covs

    means: array
        Mean for each parameter
    covs: array
        Covariance for parameters (npar, npar)
    """
    if weights is not None:
        return _extract_weighted_stats(data, weights)
    else:
        return _extract_stats(data)


def _extract_stats(data):
    import numpy as np

    ntrials = data.shape[0]
    npar = data.shape[1]

    means = np.zeros(npar, dtype="f8")
    cov = np.zeros((npar, npar), dtype="f8")

    for i in range(npar):
        means[i] = data[:, i].mean()

    num = ntrials

    for i in range(npar):
        idiff = data[:, i] - means[i]
        for j in range(i, npar):
            if i == j:
                jdiff = idiff
            else:
                jdiff = data[:, j] - means[j]

            cov[i, j] = (idiff * jdiff).sum() / (num - 1)

            if i != j:
                cov[j, i] = cov[i, j]

    return means, cov


def _extract_weighted_stats(data, weights):
    import numpy as np

    if weights.size != data.shape[0]:
        raise ValueError("weights not same size as data")

    npar = data.shape[1]

    wsum = weights.sum()

    means = np.zeros(npar, dtype="f8")
    cov = np.zeros((npar, npar), dtype="f8")

    for i in range(npar):
        dsum = (data[:, i] * weights).sum()
        means[i] = dsum / wsum

    for i in range(npar):
        idiff = data[:, i] - means[i]
        for j in range(i, npar):
            if i == j:
                jdiff = idiff
            else:
                jdiff = data[:, j] - means[j]

            wvar = (weights * idiff * jdiff).sum() / wsum
            cov[i, j] = wvar

            if i != j:
                cov[j, i] = cov[i, j]

    return means, cov
