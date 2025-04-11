ERR_TOL = 0.1
NRAND = 1000


def test_bootstrap_1d():
    """
    test against gaussian errors
    """
    import numpy as np
    from espy.bootstrap import bootstrap

    rng = np.random.default_rng(991)

    num = 10_000

    for i in range(10):
        sigma = 2.0

        data = rng.normal(scale=sigma, size=num)

        mn, err = bootstrap(data=data, nrand=NRAND, rng=rng)

        expected_err = sigma / np.sqrt(num)

        fracdiff = err / expected_err - 1
        print('fracdiff:', fracdiff)
        assert fracdiff < 0.1


def test_bootstrap_dims():
    """
    test against gaussian errors
    """
    import numpy as np
    from espy.bootstrap import bootstrap

    rng = np.random.default_rng(812)

    sigmas = np.array([2, 3, 4])

    num = 10_000

    for i in range(10):
        data = rng.normal(scale=sigmas, size=(num, sigmas.size))

        mns, errs = bootstrap(data=data, nrand=NRAND, rng=rng)

        expected_errs = sigmas / np.sqrt(num)

        fracdiff = errs / expected_errs - 1
        print('expected_errs:', expected_errs)
        print('bootstrap errs:', errs)
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

    for i in range(10):
        sums = np.zeros((nchunks, nbins))
        wsums = np.zeros((nchunks, nbins), dtype='i8')
        for i in range(nchunks):
            for ibin, sigma in enumerate(sigmas):

                num = rng.poisson(nperbin_mean)
                data = rng.normal(size=num, scale=sigma)

                sums[i, ibin] = data.sum()
                wsums[i, ibin] = num

        mns, errs = bootstrap_sums(
            sums=sums, wsums=wsums, nrand=NRAND, rng=rng,
        )

        expected_errs = sigmas / np.sqrt(wsums.sum(axis=0))

        fracdiff = errs / expected_errs - 1
        print('fracdiff:', fracdiff)
        assert np.all(fracdiff < 0.1)


if __name__ == '__main__':
    test_bootstrap_1d()
    test_bootstrap_dims()
    test_bootstrap_sums()
