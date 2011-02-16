"""
Module:
    sdsspy.util

Functions:
    photoid(run,rerun,camcol,field,id)
    photoid(array with fields):
        Create a 65-bit index from the id info.
    photoid_extract(superid): Extract run,rerun,camcol,field,id from
        a super id created using the photoid() function.


    nmgy2mag(nmgy, ivar=None):
        Convert nano-maggies to log magnitude.  If ivar is sent, the tuple
        (mag, err) is returned.
    mag2nmgy(mag, magerr=None):
        Convert from magnitudes to nano-maggies.
    nmgy2lups(nmgy, ivar=None, band=None):
        Convert from nano-maggies to luptitudes, which are asinh based
        mags.

    make_cmodelflux(recarr):
        Combine dev and exp fluxes into a composite flux.  The basic
        formula is
            fc = (1-fracdev)*flux_exp + fracdev*flux_dev
    make_cmodelmag(recarr, doerr=True, dered=False, lups=False)
        Combine dev and exp fluxes into a composite flux, and then convert
        to mags.

    dered_fluxes(extinction, flux, ivar=None):
        Deredden fluxes using the 'extinction' given in mags by
        the SDSS pipeline.  It is slightly more precise to correct
        the fluxes before converting to mags.
 

"""
import numpy
from numpy import log10, log, sqrt, where

BANDS = {'u':0,
         'g':1,
         'r':2,
         'i':3,
         'z':4}

def photoid(*args, **keys):
    """
    Name:
        photoid
    Purpose:
        Convert run,rerun,camcol,field,id to a single 64-bit superid
    Usage:
        pid = photoid(run,rerun,camcol,field,id)
        pid = photoid(recarray)
    Inputs:
        run,rerun,camcol,field,id: SDSS id info.  Note you can enter
            a subset as well as long as they occur in the expected
            order and you start with run, e.g. just run,rerun,camcol

        OR

        recarray: An array with the fields run,rerun,camcol,field,id

    Keywords:
        old: If true, use the older set of exponents.  These
        old exponents were chosen poorly and should not be used
        for new ids.

    Outputs:
        A super id combining all the id info into a single 64-bit
        number.  This is done in such a way that you can read off
        the individual ids.  E.g. if id info is
            run: 1048
            rerun: 301
            camcol: 4
            field: 125
            id: 994
        Then the result is
            10480301401250994

        You can use the photoid_extract() function to extract
        the original id info.
    """

    def photoid_usage():
        raise ValueError("Send either run,rerun,camcol,field,id or "
                         "a recarray with those fields")

    nargs = len(args)
    if nargs == 0:
        photoid_usage()

    if nargs == 1:
        arg1 = args[0]
        if isinstance(arg1, numpy.ndarray):
            # see if this has the fields, if so call sphotoid
            names = arg1.dtype.names
            if names is not None:
                return sphotoid(arg1)
            else:
                photoid_usage()
        else:
            photoid_usage()


    # order: run,rerun,camcol,field,id
    #pvals = [15,12,11,6,0]
    # 6 for run
    # 4 for rerun
    # 1 for camcol
    # 4 for field
    # 4 for id
    old = keys.get('old',False)
    if not old:
        pvals = [13,9,8,4,0]
    else:
        pvals = [15,12,11,6,0]
    if nargs > 5:
        nargs=5

    # make ten an array scalar with exactly 64 bits
    ten = numpy.array(10, dtype='i8')
    superid = None
    for i in range(nargs):
        arg = numpy.array(args[i], dtype='i8', copy=False, ndmin=1)

        if superid is None:
            superid = numpy.zeros(arg.size,dtype='i8')

        p = pvals[i]

        superid += arg*ten**p

    return superid



def sphotoid(arr, old=False):
    """
    This just extracts id info from an array with fields.  There is not
    need to call this directly since photoid will call sphotoid if needed.
    """
    names = arr.dtype.names
    if names is None:
        raise ValueError("array must have fields")

    args = []

    # we allow subsets from the left
    for name in ['run','rerun','camcol','field','id']:
        if name in names:
            tarr = numpy.array(arr[name], dtype='i8', copy=False, ndmin=1)
            args.append(tarr)
        else:
            break

    if len(args) == 0:
        raise ValueError("the struct must contain at least 'run'")
    args = tuple(args)
    return photoid( *args, old=old )

def photoid_extract(photoid, old=False):
    """
    Name:
        photoid_extract
    Purpose:
        Extract run,rerun,camcol,field,id from a super id created using
        the photoid() function.
    Usage:
        run,rerun,camcol,field,id = photoid_extract(superid)
    Inputs:
        superid: A super id created using photoid()
    Outputs:
        The SDSS id info.
    Keywords:
        old: If true, assume the old exponents were used.
    """
    if not old:
        pvals = [13,9,8,4,0]
    else:
        pvals = [15,12,11,6,0]

    # make ten an array scalar with exactly 64 bits
    ten = numpy.array(10, dtype='i8')

    run = photoid/ten**pvals[0]

    rerun = photoid/ten**pvals[1] - run*ten**(pvals[0]-pvals[1])

    camcol = \
        photoid/ten**pvals[2] - \
        run*ten**(pvals[0]-pvals[2]) - \
        rerun*ten**(pvals[1]-pvals[2])

    field = \
        photoid/ten**pvals[3] - \
        run*ten**(pvals[0]-pvals[3]) - \
        rerun*ten**(pvals[1]-pvals[3]) - \
        camcol*ten**(pvals[2]-pvals[3])

    id =  \
        photoid/ten**pvals[4] - \
        run*ten**(pvals[0]-pvals[4]) - \
        rerun*ten**(pvals[1]-pvals[4]) - \
        camcol*ten**(pvals[2]-pvals[4]) - \
        field*ten**(pvals[3]-pvals[4])

    return run,rerun,camcol,field,id




def nmgy2mag(nmgy, ivar=None):
    """
    Name:
        nmgy2mag
    Purpose:
        Convert SDSS nanomaggies to a log10 magnitude.  Also convert
        the inverse variance to mag err if sent.  The basic formulat
        is 
            mag = 22.5-2.5*log_{10}(nanomaggies)
    Calling Sequence:
        mag = nmgy2mag(nmgy)
        mag,err = nmgy2mag(nmgy, ivar=ivar)
    Inputs:
        nmgy: SDSS nanomaggies.  The return value will have the same
            shape as this array.
    Keywords:
        ivar: The inverse variance.  Must have the same shape as nmgy.
            If ivar is sent, then a tuple (mag,err) is returned.

    Outputs:
        The magnitudes.  If ivar= is sent, then a tuple (mag,err)
        is returned.

    Notes:
        The nano-maggie values are clipped to be between 
            [0.001,1.e11]
        which corresponds to a mag range of 30 to -5
    """
    nmgy = numpy.array(nmgy, ndmin=1, copy=False)

    nmgy_clip = numpy.clip(nmgy,0.001,1.e11)

    mag = nmgy_clip.copy()
    mag[:] = 22.5-2.5*log10(nmgy_clip)

    if ivar is not None:

        ivar = numpy.array(ivar, ndmin=1, copy=False)
        if ivar.shape != nmgy.shape:
            raise ValueError("ivar must be same shape as input nmgy array")

        err = mag.copy()
        err[:] = 9999.0

        w=where( ivar > 0 )

        if w[0].size > 0:
            err[w] = sqrt(1.0/ivar[w])

            a = 2.5/log(10)
            err[w] *= a/nmgy_clip[w]

        return mag, err
    else:
        return mag

def mag2nmgy(mag, magerr=None):
    """
    Name:
        mag2nmgy
    Purpose:
        Convert from magnitudes to nano-maggies.  The basic formula
        is 
            mag = 22.5-2.5*log_{10}(nanomaggies)
        The mag error can optionally be sent, in which case the inverse
        variance of the nanomaggies is returned.
    Calling Sequence:
        nmgy = mag2nmgy(mag)
        nmgy,ivar = mag2nmgy(mag, magerr=magerr)
            
    """
    mag = numpy.array(mag, ndmin=1, copy=False)

    nmgy = 10.0**( (22.5-mag)/2.5 )
    if magerr is not None:
        ivar = nmgy.copy()
        ivar[:] = 0.0

        w = numpy.where( (nmgy > 0) & (magerr > 0) )

        if w[0].size > 0:
            a = 2.5/log(10)
            ivar[w] = ( a/nmgy[w]/magerr[w] )**2

        return nmgy, ivar
    else:
        return nmgy

def nmgy2lups(nmgy, ivar=None, band=None):
    """
    Name:
        nmgy2lups
    Purpose:
        Convert from nano-maggies to luptitudes, which are asinh based
        mags.  The default parameters for SDSS are used.
    Calling Sequence:
        lup = nmgy2lups(nmgy)
        lup,err = nmgy2lups(nmgy, ivar=ivar)
    Inputs:
        nmgy: SDSS nanomaggies.  Can either be a [5,Nobj] array or
            an array for a single band, in which case the band must
            be given.
    Keywords:
        ivar: The inverse variance.  Must have the same shape as nmgy.
            If ivar is sent, then a tuple (lup,luperr) is returned.

    Outputs:
        The luptitudes as asinh values.  If ivar= is sent, a tuple
        is returned (lup,luperr)
    """
    s = nmgy.shape
    if ivar is not None:
        sivar = ivar.shape
        if len(sivar) != len(s):
            raise ValueError("ivar and fluxes must be same shape")
        for i in xrange(len(s)):
            if sivar[i] != s[i]:
                raise ValueError("ivar and fluxes must be same shape")

    if len(s) == 2:
        if s[1] != 5:
            raise ValueError("Either enter a 1-d array or a (nobj, 5) array")
        nband = 5
        band=[0,1,2,3,4]
    else:
        if band is None:
            raise ValueError("For 1-d input, specify a band in [0,4]")
        nband = 1
        try:
            if len(band) != 1:
                raise ValueError("for 1-d input, enter a single band")
        except:
            band = [band]

    # make sure band values makes sense
    for b in band:
        if b not in [0,1,2,3,4]:
            raise ValueError("band must be in [0,4]")

    lups = numpy.array( nmgy, copy=True )
    lups[:] = -9999.0
    if ivar is not None:
        lups_err = numpy.array(ivar, copy=True)
        lups_err[:] = -9999.0

    for b in band:
        if nband == 1:
            lups[:] = _nmgy2lups_1band(nmgy, b)
            if ivar is not None:
                lups_err[:] = _ivar2luperr_1band(nmgy, ivar, b)
        else:
            lups[:,b] = _nmgy2lups_1band(nmgy[:,b], b)
            if ivar is not None:
                lups_err[:,b] = _ivar2luperr_1band(nmgy[:,b], ivar[:,b], b)

    if ivar is not None:
        return lups, lups_err
    else:
        return lups

_bvalues=[1.4, 0.9, 1.2, 1.8, 7.4]
_log10 = numpy.log(10.0)

ln10_min10 = -23.02585
def _nmgy2lups_1band(nmgy, band):
    b=_bvalues[band]
    lups = 2.5*(10.0-numpy.log10(b)) - 2.5*numpy.arcsinh(5.0*nmgy/b)/_log10
    return lups

def _ivar2luperr_1band(nmgy, ivar, band):
    b=_bvalues[band]
    lups_err = numpy.array(ivar, copy=True)
    lups_err[:] = -9999.0
    w,=numpy.where(ivar > 0.0)
    if w.size > 0:
        terr = 1.0/numpy.sqrt(ivar[w])
        lups_err[w] = 2.5*terr
        lups_err[w] /= 0.2*b*_log10*numpy.sqrt(1.0 + (5.0*nmgy[w]/b)**2 )
    return lups_err

def make_cmodelflux(objs, doivar=True):
    """
    Name:
        make_cmodelfllux
    Purpose:
        Combine dev and exp fluxes into a composite flux.  The basic
        formula is
            fc = (1-fracdev)*flux_exp + fracdev*flux_dev
    Calling Sequence:
        cmodel = make_cmodelflux(str, doiver=True)
    Inputs:
        str: A recarray with fields 'devflux','expflux' and one of 
            'fracdev' or 'fracpsf'.  fracpsf is allowed because
            the PHOTO pipeline now puts the value of fracdev into
            the old fracpsf field.
    Outputs:
        A [5,Nobj] array holding the composite model flux, or in
        the case that doivar=True a tuple (flux,ivar).
    """
    devflux = objs['devflux']
    expflux = objs['expflux']
    if 'fracdev' in objs.dtype.names:
        fracdev = objs['fracdev']
    else:
        fracdev = objs['fracpsf']

    if not doivar:
        flux = devflux*fracdev + expflux*(1.0-fracdev)
        return flux


    devflux_ivar = objs['devflux_ivar']
    expflux_ivar = objs['expflux_ivar']

    # get same shape, type
    flux = objs['devflux'].copy()
    ivar = objs['devflux'].copy()
    flux[:] = -9999.0
    ivar[:] = 0.0

    for band in [0,1,2,3,4]:
        bflux, bivar = _make_cmodelflux_1band(fracdev[:,band],
                                              devflux[:,band],
                                              devflux_ivar[:,band],
                                              expflux[:,band],
                                              expflux_ivar[:,band])
        flux[:,band] = bflux
        ivar[:,band] = bivar

    return flux, ivar

def _make_cmodelflux_1band(fracdev, dev, dev_ivar, exp, exp_ivar):
    fracexp = 1.0-fracdev

    flux = dev*fracdev + exp*fracexp

    ivar = exp.copy()
    ivar[:] = 0.0

    w,=where( (exp_ivar > 0) & (dev_ivar > 0) )
    if w.size > 0:
        # measurements are correlated... doing a weighted average
        err2 = fracdev[w]/dev_ivar[w] + fracexp[w]/exp_ivar[w]

        w2, = where( err2 > 0)
        if w2.size > 0:
            w=w[w2]
            ivar[w] = 1.0/err2[w2]

    return flux,ivar

def make_cmodelmag(objs, doerr=True, dered=False, lups=False):
    """
    Name:
        make_cmodelmag
    Purpose:
        Combine dev and exp fluxes into a composite flux.  The basic
        formula is
            fc = (1-fracdev)*flux_exp + fracdev*flux_dev
        Then convert the flux to mags using 
            mag = 22.5-2.5*log_{10}(fc)
    Calling Sequence:
        cmodelmag = make_cmodelmag(str, doerr=True, dered=False, lups=False)
    Inputs:
        str: A recarray with fields 'devflux','expflux' and one of 
            'fracdev' or 'fracpsf'.  fracpsf is allowed because
            the PHOTO pipeline now puts the value of fracdev into
            the old fracpsf field.

            If dred=True then the 'extinction' field must be present.
    Keywords:
        doerr: Return a tuple (mag,err)
        dered: Apply an extinction correction.
        lups: Return luptitudes instead of log mags.
    """

    if doerr:
        flux, ivar = make_cmodelflux(objs, doivar=True)

        if dered:
            flux, ivar = dered_fluxes(objs['extinction'], flux, ivar)
        if lups:
            mag, err = nmgy2lups(flux, ivar=ivar)
        else:
            mag, err = nmgy2mag(flux, ivar=ivar)
        return mag, err
    else:
        flux = make_cmodelflux(objs, doivar=False)

        if dered:
            flux = dered_fluxes(objs['extinction'], flux)
        if lups:
            mag = nmgy2lups(flux)
        else:
            mag = nmgy2mag(flux)
        return mag



def dered_fluxes(extinction, flux, ivar=None):
    """
    Adam says this is more accurate at the faint end.
    """
    exponent = 0.4*extinction
    flux_correct = (10.0**exponent)*flux

    if ivar is None:
        return flux_correct
    else:
        exponent = -0.8*extinction
        ivar_correct = (10.0**exponent)*ivar

        return flux_correct, ivar_correct










def sdss_wrap(ra):
    """
    Take all points with ra > 300 and put them at negative ra.
    This makes it easier to plot up ra and dec in the SDSS
    """
    ranew = numpy.array(ra, dtype='f8', copy=True, ndmin=1)
    w,=numpy.where(ra > 300)
    if w.size > 0:
        ranew[w] = ra[w]-360.0

    return ranew


