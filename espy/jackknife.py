"""
code for jackknifing 1 dimensional data
"""
import numpy as np
from numba import njit


def jackknife(
    *, data, weights=None, chunksize=1, err_err=False,
):
    if weights is not None:
        assert data.size == weights.size, "weights and data must be same size"
        mn, err = _wjackknife(data, weights, chunksize)
    else:
        mn, err = _jackknife(data, chunksize)

    if err_err:
        nchunks = data.size // chunksize
        errors = np.zeros(nchunks)
        logic = np.ones(data.size, dtype=bool)

        for i in range(nchunks):
            logic[:] = True

            start = i * chunksize
            end = (i + 1) * chunksize

            logic[start:end] = False

            if weights is not None:
                imn, ierr = _wjackknife(data[logic], weights[logic], chunksize)
            else:
                imn, ierr = _jackknife(data[logic], chunksize)
            errors[i] = ierr

        fac = (nchunks - 1) / nchunks

        err_err_cov = fac * ((err - errors)**2).sum()
        err_err = np.sqrt(err_err_cov)
        return mn, err, err_err
    else:
        return mn, err


def jackknife_ratio(
    *, data1, data2, weights1=None, weights2=None, chunksize=1,
):
    assert data1.size == data2.size, "data1 and data2 must be same size"

    if weights1 is not None or weights2 is not None:
        assert weights2 is not None and weights1 is not None, (
            "send both weights1 and weights2"
        )
        assert data1.size == weights1.size, (
            "data and weights must be same size"
        )
        assert weights1.size == weights2.size, (
            "data and weights must be same size"
        )

        return _wjackknife_ratio(data1, weights1, data2, weights2, chunksize)
    else:
        return _jackknife_ratio(data1, data2, chunksize)


@njit
def _jackknife(data, chunksize):

    nchunks = data.size//chunksize

    num = data.size

    dsum = data.sum()
    mn = dsum/data.size

    mns = np.zeros(nchunks)

    for i in range(nchunks):
        beg = i*chunksize
        end = (i+1)*chunksize

        subsum = 0.0
        for j in range(beg, end):
            subsum += data[j]

        tsum = dsum - subsum
        mns[i] = tsum/(num-chunksize)

    fac = (nchunks-1)/float(nchunks)
    var = fac*(((mns - mn)**2).sum())

    err = np.sqrt(var)
    return mn, err


@njit
def _wjackknife(data, weights, chunksize):

    nchunks = data.size//chunksize

    dsum = (data * weights).sum()
    wsum = weights.sum()

    mn = dsum/wsum

    mns = np.zeros(nchunks)

    for i in range(nchunks):

        beg = i*chunksize
        end = (i+1)*chunksize

        subsum = 0.0
        wsubsum = 0.0
        for j in range(beg, end):
            subsum += data[j] * weights[j]
            wsubsum += weights[j]

        tsum = dsum - subsum
        twsum = wsum - wsubsum

        mns[i] = tsum/twsum

    fac = (nchunks-1)/float(nchunks)
    var = fac*(((mns - mn)**2).sum())

    err = np.sqrt(var)
    return mn, err


@njit
def _jackknife_ratio(data1, data2, chunksize):

    nchunks = data1.size//chunksize

    num = data1.size

    dsum1 = data1.sum()
    dsum2 = data2.sum()

    mn1 = dsum1/num
    mn2 = dsum2/num

    rat = mn1/mn2

    rats = np.zeros(nchunks)

    for i in range(nchunks):

        beg = i*chunksize
        end = (i+1)*chunksize

        subsum1 = 0.0
        subsum2 = 0.0
        for j in range(beg, end):
            subsum1 += data1[j]
            subsum2 += data2[j]

        tdsum1 = dsum1 - subsum1
        tdsum2 = dsum2 - subsum2

        tmn1 = tdsum1/(num - chunksize)
        tmn2 = tdsum2/(num - chunksize)

        rats[i] = tmn1/tmn2

    fac = (nchunks-1)/float(nchunks)
    rat_var = fac*(((rat - rats)**2).sum())

    rat_err = np.sqrt(rat_var)
    return rat, rat_err


@njit
def _wjackknife_ratio(data1, weights1, data2, weights2, chunksize):

    nchunks = data1.size//chunksize

    wsum1 = weights1.sum()
    wsum2 = weights2.sum()

    dsum1 = (data1 * weights1).sum()
    dsum2 = (data2 * weights2).sum()

    mn1 = dsum1/wsum1
    mn2 = dsum2/wsum2

    rat = mn1/mn2

    rats = np.zeros(nchunks)

    for i in range(nchunks):

        subsum1 = 0.0
        subsum2 = 0.0
        wsubsum1 = 0.0
        wsubsum2 = 0.0

        beg = i*chunksize
        end = (i+1)*chunksize

        for j in range(beg, end):
            subsum1 += data1[j] * weights1[j]
            subsum2 += data2[j] * weights2[j]
            wsubsum1 += weights1[j]
            wsubsum2 += weights2[j]

        tdsum1 = dsum1 - subsum1
        tdsum2 = dsum2 - subsum2

        twsum1 = wsum1 - wsubsum1
        twsum2 = wsum2 - wsubsum2

        tmn1 = tdsum1/twsum1
        tmn2 = tdsum2/twsum2

        rats[i] = tmn1/tmn2

    fac = (nchunks-1)/float(nchunks)
    rat_var = fac*(((rat - rats)**2).sum())

    rat_err = np.sqrt(rat_var)
    return rat, rat_err


def test_jackknife():
    seed = 5123
    rng = np.random.RandomState(seed)

    data = rng.normal(size=100000)

    mn = data.mean()
    err = data.std()/np.sqrt(data.size)

    jmn, jerr = jackknife(data=data)
    jmn100, jerr100 = jackknife(data=data, chunksize=100)

    print('simple:', mn, err)
    print('jack 1:', jmn, jerr)
    print('jack 100:', jmn100, jerr100)


def test_jackknife_err_err():
    rng = np.random.RandomState()

    nper = 100000
    def make_data():
        return rng.normal(size=100000)

    _, err, err_err_predicted = jackknife(
        data=make_data(), chunksize=100, err_err=True,
    )

    ntrials = 100
    errors = np.zeros(ntrials)
    for i in range(ntrials):
        _, errors[i] = jackknife(data=make_data(), chunksize=100)

    err_err_measured = errors.std()
    print(f'err: {err:g}')
    print(f'err_err predicted: {err_err_predicted:g}')
    print(f'err_err measured: {err_err_measured:g}')
    print(f'err_err sqrt: {err / np.sqrt(2 * nper):g}')


def test_wjackknife():
    import esutil as eu
    seed = 800
    num = 50_000
    rng = np.random.RandomState(seed)

    mn_low = 0
    std_low = 1

    mn_high = 0.5
    std_high = 2

    data_low = rng.normal(size=num, loc=mn_low, scale=std_low)
    data_high = rng.normal(size=num, loc=mn_high, scale=std_high)

    err_low = data_low*0 + std_low
    err_high = data_high*0 + std_high

    data = np.hstack((data_low, data_high))
    err = np.hstack((err_low, err_high))

    # need to scramble, or else we will see the different samples
    # too visibly
    r = np.arange(data.size)
    rng.shuffle(r)
    data = data[r]
    err = err[r]

    weights = 1.0/err**2
    mn, err = eu.stat.wmom(data, weights, calcerr=True)

    jmn, jerr = jackknife(data=data, weights=weights)
    jmn100, jerr100 = jackknife(data=data, weights=weights, chunksize=100)

    print('wmom:', mn, err)
    print('jack 1:', jmn, jerr)
    print('jack 100:', jmn100, jerr100)


def test_jackknife_ratio():
    seed = 99
    num = 50_000
    rng = np.random.RandomState(seed)

    mn_low = 0
    std_low = 1

    mn_high = 0.5
    std_high = 2

    data_low = rng.normal(size=num, loc=mn_low, scale=std_low)
    data_high = rng.normal(size=num, loc=mn_high, scale=std_high)

    err_low = data_low*0 + std_low
    err_high = data_high*0 + std_high

    rat = data_low.mean()/data_high.mean()

    jrat, jrat_err = jackknife_ratio(
        data1=data_low, data2=data_high,
    )
    jrat100, jrat_err100 = jackknife_ratio(
        data1=data_low, data2=data_high, chunksize=100,
    )

    print('simple: %g' % rat)
    print('jack: %g +/- %g' % (jrat, jrat_err))
    print('jack100: %g +/- %g' % (jrat100, jrat_err100))

    weights_low = 1.0/err_low**2
    weights_high = 1.0/err_high**2

    wjrat, wjrat_err = jackknife_ratio(
        data1=data_low, data2=data_high,
        weights1=weights_low, weights2=weights_high,
    )
    wjrat100, wjrat_err100 = jackknife_ratio(
        data1=data_low, data2=data_high, chunksize=100,
        weights1=weights_low, weights2=weights_high,
    )

    print('wjack: %g +/- %g' % (wjrat, wjrat_err))
    print('wjack100: %g +/- %g' % (wjrat100, wjrat_err100))


if __name__ == '__main__':
    print('-'*70)
    print('jackknife')
    test_jackknife()

    print('-'*70)
    print('weighted jackknife')
    test_wjackknife()

    print('-'*70)
    print('jackknife ratio')
    test_jackknife_ratio()

    print('-'*70)
    print('jackknife err_err')
    test_jackknife_err_err()
