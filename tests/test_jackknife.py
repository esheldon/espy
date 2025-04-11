ERR_TOL = 0.1


def test_jackknife():
    """
    test against gaussian errors
    """
    import numpy as np
    from espy.jackknife import jackknife

    rng = np.random.default_rng(991)

    num = 10_000

    for i in range(10):
        sigma = 2.0

        data = rng.normal(scale=sigma, size=num)

        mn, err = jackknife(data=data)

        expected_err = sigma / np.sqrt(num)

        fracdiff = err / expected_err - 1
        print('fracdiff:', fracdiff)
        assert fracdiff < 0.1


if __name__ == '__main__':
    test_jackknife()
