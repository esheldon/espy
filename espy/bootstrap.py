def bootstrap(data, nrand, weights=None, rng=None):
    """
    Bootstrap the data to get the error on the mean

    Parameters
    ----------
    data: array
        An array of data
    nrand: int
        Number of times to sample the data, e.g. 100.
    weight: array, optional
        Optional weights for the data
    rng: np.random.default_rng
        Default random number generator.  If not sent, one will
        be created, seeded by entropy from the OS

    Returns
    -------
    mean, err:  The mean and uncertainty on the mean.
    """

    import numpy as np

    if weights is None:
        weights = np.ones(data.size)

    vals = np.zeros(nrand)
    for i in range(nrand):
        ind = rng.integers(0, data.size, size=data.size)
        vals[i] = (data[ind] * weights[ind]).sum() / weights[ind].sum()

    return vals.mean(), vals.std()
