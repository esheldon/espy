def test_bootstrap():
    """
    test against gaussian errors
    """
    import numpy as np
    from espy.bootstrap import bootstrap

    rng = np.random.default_rng(812)

    sigmas = np.array([2, 3, 4])

    num = 10_000

    data = rng.normal(
        scale=sigmas,
        size=(num, sigmas.size),
    )

    nrand = 100
    errs = bootstrap(data=data, nrand=nrand, rng=rng)

    expected_errs = sigmas / np.sqrt(num)

    fracdiff = errs / expected_errs - 1
    print('fracdiff:', fracdiff)
    assert np.all(fracdiff < 0.1)


def test_bootstrap_sums():
    """
    test against gaussian errors
    """
    import numpy as np
    from espy.bootstrap import bootstrap_sums

    rng = np.random.default_rng(1998)

    sigmas = np.array([2, 3, 4])
    nchunks = 1000

    nbins = sigmas.size

    nperbin_mean = 1000

    sums = np.zeros((nchunks, nbins))
    wsums = np.zeros((nchunks, nbins), dtype='i8')
    for i in range(nchunks):
        for ibin, sigma in enumerate(sigmas):

            num = rng.poisson(nperbin_mean)
            data = rng.normal(size=num, scale=sigma)

            sums[i, ibin] = data.sum()
            wsums[i, ibin] = num

    nrand = 100
    errs = bootstrap_sums(sums=sums, wsums=wsums, nrand=nrand, rng=rng)

    expected_errs = sigmas / np.sqrt(wsums.sum(axis=0))

    fracdiff = errs / expected_errs - 1
    print('fracdiff:', fracdiff)
    assert np.all(fracdiff < 0.1)


if __name__ == '__main__':
    test_bootstrap()
    test_bootstrap_sums()
