def bootstrap(data, nrand, weights=None, rng=None):
    """
    Bootstrap the data to get the error on the mean

    Parameters
    ----------
    data: array
        An array of data.  Can be one dimensional or of shape (npoints, ndim)
    nrand: int
        Number of times to sample the data, e.g. 100.
    weight: array, optional
        Optional weights for the data
    rng: np.random.default_rng
        Default random number generator.  If not sent, one will
        be created, seeded by entropy from the OS

    Returns
    -------
    err:  The mean and uncertainty on the mean.
    """

    import numpy as np

    if weights is None:
        weights = np.ones(data.shape)

    if weights.shape != data.shape:
        raise ValueError(
            f'weights shape {weights.shape} does not match '
            f'data shape {data.shape}'
        )

    if rng is None:
        rng = np.random.default_rng()

    if len(data.shape) > 1:
        nbin = data.shape[1]
        vals = np.zeros((nrand, nbin))
    else:
        vals = np.zeros(nrand)

    for i in range(nrand):
        ind = rng.integers(0, data.shape[0], size=data.shape[0])
        dsum = (data[ind] * weights[ind]).sum(axis=0)
        wsum = weights[ind].sum(axis=0)
        vals[i] = dsum / wsum

    return vals.std(axis=0)


def bootstrap_sums(sums, wsums, nrand, rng=None):
    """
    Bootstrap the data to get the error on the mean

    Parameters
    ----------
    sums: array
        An array of data representing sums for the numerator
    wsums: array
        An array of data representing sums for the denominator
    nrand: int
        Number of times to sample the data, e.g. 100.
    rng: np.random.default_rng
        Default random number generator.  If not sent, one will
        be created, seeded by entropy from the OS

    Returns
    -------
    err:  The mean and uncertainty on the mean.
    """

    import numpy as np

    if rng is None:
        rng = np.random.default_rng()

    if wsums.shape != sums.shape:
        raise ValueError(
            f'wsums shape {wsums.shape} does not match '
            f'sumsshape {sums.shape}'
        )

    if len(sums.shape) > 1:
        nbin = sums.shape[1]
        vals = np.zeros((nrand, nbin))
    else:
        vals = np.zeros(nrand)

    for i in range(nrand):
        ind = rng.integers(0, sums.shape[0], size=sums.shape[0])
        vals[i] = sums[ind].sum(axis=0) / wsums[ind].sum(axis=0)

    return vals.std(axis=0)
