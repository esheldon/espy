"""
Module:
    imsim

Purpose:
    Tools to create fake images.

Star and galaxy models:

    mom2disk: 
        Create images from a specified model with specified covariance matrix.
    mom2dgauss:
        Double gaussians

Classes:
    ObjectSimulator:
        Create an object convolved with a PSF.  Can use models for objects
        and psf and also, if sdsspy is available, an actual random SDSS
        psf.

"""
import numpy
import images

try:
    import scipy.signal
except:
    pass

def mom2dgauss(Irr1, Irc1, Icc1, Tratio, fluxfrac1, dims, 
               counts=1.0, cen=None, all=False, dtype='f4'):
    """
    Make a double gaussian image.

    Parameters
    ----------
    Irr1,Irc1,Icc1: scalar
        The moments of the first gaussian.

    Tratio: scalar
        The moment ratio of the second to first gaussians.
        I2 = Tratio*I1

    fluxfrac1: scalar
        The fraction of the flux in the first gaussian

    dims: sequence, length 2
        The dimensions of the image [nrow,ncol]

    counts: scalar, optional
        The total counts in the image.  Default 1.0
    cen: sequence, length 2, optional
        the center of the image, default [(dims[0]-1)/2, (dims[1]-1)/2]

    all: boolean, optional
        If True, return a tuple (image, image1, image). Default False

    dtype: data-type, optional
        The data type of the image, default 'f4'

    Returns
    -------
    image: 2-d array
        The image. If all=True, it will be a tuple of images
        (image, image1, image2)
    """

    if fluxfrac1 > 1.0 or fluxfrac1 < 0.0:
        raise ValueError("fluxfrac1 must be in [0,1]")

    counts1 = counts*fluxfrac1
    counts2 = counts*(1.0-fluxfrac1)

    image1 = mom2disk('gauss', Irr1, Irc1, Icc1, dims, 
                      counts=counts1, cen=cen, dtype=dtype)
    image2 = mom2disk('gauss', Irr1*Tratio, Irc1*Tratio, Icc1*Tratio, dims, 
                      counts=counts2, cen=cen, dtype=dtype)

    image = image1 + image2

    if all:
        return image, image1, image2
    else:
        return image


def mom2disk(model, Irr, Irc, Icc, dims, counts=1., cen=None, dtype='f4'):
    """
    Name:
        mom2disk
    Purpose:
        Create an image from the input covariance matrix assuming a model, e.g.
        gaussian exponential disk, dev.

    Calling Sequence:
        image = mom2disk(model, Irr, Irc, Icc, dims, counts=1, cen=None, dtype='f4')
    Inputs:
        model: 'gauss', 'exp', or 'dev'
        Irr,Irc,Icc: The moments <I*row**2>, <I*row*col>, <I*col**2>
            r stands for row
            c stands for column
        dims: The image size [nrows, ncolumns]
    Keywords:
        counts: The total counts in the image, default is 1.0
        cen: The center of the image, default is [(nrows-1)/2., (ncolumns-1)/2.]
        dtype: The data type, default is 'f4'

    Output:
        The image returned as a numpy array.  
        
        Numpy arrays are row-major storage and are addressed image[row,column].
        This means in C you would address it with:

                array[row,col] -> (array + row*ncols + col).   
                
        Note the rightmost index, the column, varies fastest.

    """

    det = Irr*Icc - Irc**2
    if det == 0.0:
        raise RuntimeError("Determinant is zero")

    if len(dims) != 2:
        raise ValueError("dims must have two elements")
    nrows=int(dims[0])
    ncols=int(dims[1])

    if cen is None:
        rowcen = (nrows-1.)/2.
        colcen = (ncols-1.)/2.
        cen=[rowcen,colcen]
    else:
        rowcen = cen[0]
        colcen = cen[1]

    row, col, rr, index = get_flatgrid(Irr, Irc, Icc, dims, cen, dtype=dtype)

    model = model.lower()
    if model == 'gauss':
        rr = 0.5*rr/det
    elif model == 'exp':
        rr = numpy.sqrt(rr/det)
    elif model == 'dev':
        rr = numpy.sqrt(rr)
        rr = 7.67*( (rr/det)**(.25) -1 )
    else: 
        raise ValueError("model must be one of gauss, exp, or dev")


    disk=numpy.zeros(nrows*ncols, dtype=dtype)
    w, = numpy.where(rr < 10.8) 
    if w.size > 0:
        disk[index[w]] = numpy.exp(-rr[w])

    norm = disk.sum()
    if norm > 0:
        disk *= (counts/norm)

    disk = disk.reshape(nrows, ncols)
    return disk

def get_flatgrid(Irr, Irc, Icc, dims, cen, dtype='f4'):
    nrows=dims[0]
    ncols=dims[1]
    index=numpy.arange(nrows*ncols, dtype='i4')
    row = numpy.zeros(nrows*ncols, dtype=dtype)
    col = numpy.zeros(nrows*ncols, dtype=dtype)
    rr = numpy.zeros(nrows*ncols, dtype=dtype)

    col[:] = index % ncols
    row[:] = index/ncols

    rowcen=cen[0]
    colcen=cen[1]

    row -= rowcen
    col -= colcen
    rr[:] = row**2*Icc -2*row*col*Irc + col**2*Irr

    return row, col, rr, index

def test_gauss(Irr, Irc, Icc, dims, cen, dtype='f4'):
    """
    Test gaussian image.  Compare unweighted moments to input
    """
    row, col, rr, index=get_flatgrid(Irr, Irc, Icc, dims, cen, dtype=dtype)

    # note im is normalized
    im = mom2disk('gauss', Irr, Irc, Icc, dims, cen=cen, dtype=dtype)
    imflat=im.ravel()
    
    Irr_meas = (imflat*row**2).sum()
    Irc_meas = (imflat*row*col).sum()
    Icc_meas = (imflat*col**2).sum()

    Irr_pdiff = (Irr_meas-Irr)/Irr
    Irc_pdiff = (Irc_meas-Irc)/Irc
    Icc_pdiff = (Icc_meas-Icc)/Icc
    print("""
    Irr: {Irr} meas: {Irr_meas} %%diff: {Irr_pdiff}
    Irc: {Irc} meas: {Irc_meas} %%diff: {Irc_pdiff}
    Irr: {Icc} meas: {Icc_meas} %%diff: {Icc_pdiff}
          """.format(Irr=Irr,Irc=Irc,Icc=Icc,
                     Irr_meas=Irr_meas,Irc_meas=Irc_meas,Icc_meas=Icc_meas,
                     Irr_pdiff=Irr_pdiff,Irc_pdiff=Irc_pdiff,Icc_pdiff=Icc_pdiff))


def plot_dgauss(Irr1, Irc1, Icc1, Tratio, fluxfrac1, dims):
    import biggles
    im, im1, im2 = mom2dgauss(Irr1,Irc1,Icc1,Tratio,fluxfrac1,dims,
                              all=True)
    cen=[(im.shape[0]-1)/2.0, (im.shape[1]-1)/2.]

    # cross-section across rows
    imrows = im[:, cen[1]]
    imcols = im[cen[0], :]

    imrows1 = im1[:, cen[1]]
    imcols1 = im1[cen[0], :]

    imrows2 = im2[:, cen[1]]
    imcols2 = im2[cen[0], :]

    tab=biggles.Table(2,2)

    # cross-section of rows
    rowplt = biggles.FramedPlot()
    hrows = biggles.Histogram(imrows, color='blue')
    hrows.label = 'image'
    hrows1 = biggles.Histogram(imrows1, color='red')
    hrows1.label = 'image1'
    hrows2 = biggles.Histogram(imrows2, color='orange')
    hrows2.label = 'image2'

    rowkey = biggles.PlotKey(0.1, 0.9, [hrows, hrows1, hrows2])

    rowplt.add(hrows,hrows1,hrows2,rowkey)
    rowplt.xlabel = 'rows through center'

    # cross-section of rows
    colplt = biggles.FramedPlot()
    hcols = biggles.Histogram(imcols, color='blue')
    hcols1 = biggles.Histogram(imcols1, color='red')
    hcols2 = biggles.Histogram(imcols2, color='orange')

    colplt.add(hcols,hcols1,hcols2)
    colplt.xlabel = 'cols through center'

    # contours of images
    cplt =biggles.FramedPlot()
    levels = 7
    c = biggles.Contours(im, color='red')
    c.label = 'image'
    c.levels = levels
    c1 = biggles.Contours(im1, color='blue')
    c1.label = 'image1'
    c1.levels = levels
    c2 = biggles.Contours(im2, color='orange')
    c2.label = 'image2'
    c2.levels = levels

    ckey = biggles.PlotKey(0.1, 0.9, [c,c1,c2])

    #cplt.add(c,c1,c2,ckey)
    cplt.add(c,c1,c2)

    tab[0,0] = rowplt
    tab[0,1] = colplt
    tab[1,0] = cplt
    tab.show()


class ConvolvedImage(dict):
    """

    Simulate an object convolved with a psf.  

    The psfmodel can either be an image or a string describing
    a model, 'gauss' or 'degauss' for now.

    If psfmodel is an image , it's stats are measured if not entered

    If the psfmodel is an image, it must be normalized and sky
    subtracted
    
    For exponential disks, a the image is created to be 2*7*sigma
    wide, where sigma is determined from the expected size after
    convolution.  For gaussians, a this is 2*4*sigma

    """

    def __init__(self, 
                 objmodel, Irr, Irc, Icc, 
                 psfmodel, # can be an image or a string
                 Irr_psf=None,Irc_psf=None,Icc_psf=None, 
                 Tratio=None, fluxfrac1=None, # these for double gaussians
                 verbose=False):

        self.verbose=verbose
        self['objmodel'] = objmodel
        self['Irr'] = Irr
        self['Irc'] = Irc
        self['Icc'] = Icc

        # dims, cen to be determined
        #self['dims'] = dims
        #self['cen'] = cen
        self['counts'] = 1000.0

        self['psfmodel'] = psfmodel
        self['Irr_psf'] = Irr_psf
        self['Irc_psf'] = Irc_psf
        self['Icc_psf'] = Icc_psf
        self['Tratio'] = Tratio
        self['fluxfrac1'] = fluxfrac1

        self.make_psf()
        self.make_image0()
        self.make_image()


    def make_psf(self):
        if isinstance(self['psfmodel'], (str,unicode)):
            self._make_model_psf()
        else:
            self._stat_image_psf()

    def _stat_image_psf(self):
        import admom
        # measure the moments using admom if they were
        # not sent

        psf = self['psfmodel']
        self.psf = psf

        # note we assume the center is the middle.
        cen = [(psf.shape[0]-1)/2., (psf.shape[1]-1)/2.]
        self['psfcen'] = cen

        if self['Irr_psf'] is None:
            import admom

            if (psf.shape[0] % 2) == 0 or (psf.shape[1] % 2) == 0:
                raise ValueError("psf must be odd size in each dim")

            out = admom.admom(psf, cen[0], cen[1])
            if out['whyflag'] != 0:
                raise ValueError("Failed to characterize PSF")
            self['Irr_psf'] = out['Irr']
            self['Irc_psf'] = out['Irc']
            self['Icc_psf'] = out['Icc']

        if self.verbose:
            print("Input PSF model stats:")
            print("  Irr_psf:",self['Irr_psf'])
            print("  Irc_psf:",self['Irc_psf'])
            print("  Icc_psf:",self['Icc_psf'])
            print("  dims:   ",self.psf.shape)



    def _make_model_psf(self):
        psfmodel = self['psfmodel']
        Irr_psf = self['Irr_psf']
        Irc_psf = self['Irc_psf']
        Icc_psf = self['Icc_psf']

        Tpsf = Irr_psf + Icc_psf
        psf_sigma = self.mom2sigma(Tpsf)

        if self.verbose:
            print("Making PSF model '%s'" % psfmodel)
            print("  Irr_psf:",Irr_psf)
            print("  Irc_psf:",Irc_psf)
            print("  Icc_psf:",Icc_psf)

        if psfmodel == 'dgauss':
            Tratio = self['Tratio']
            fluxfrac1 = self['fluxfrac1']

            if Tratio is None or fluxfrac1 is None:
                raise ValueError("Enter Tratio and fluxfrac1 for "
                                 "double gauss psf")

            psf_sigma2 = self.mom2sigma(Tratio*Tpsf)
            psf_sigma_max = max(psf_sigma, psf_sigma2)

            psf_imsize = int( round(2*4*psf_sigma_max) )
            # MUST BE ODD!
            if (psf_imsize % 2) == 0:
                psf_imsize+=1
            psfdims = [psf_imsize,psf_imsize]

            if self.verbose:
                print("  Tratio:",Tratio)
                print("  fluxfrac1:",fluxfrac1)
                print("  dims:",psfdims)

            self.psf = mom2dgauss(Irr_psf, 
                                  Irc_psf, 
                                  Icc_psf, 
                                  Tratio, 
                                  self['fluxfrac1'], 
                                  psfdims)
            self['psfcen'] = [(psfdims[0]-1)/2., (psfdims[1]-1)/2.]
        elif psfmodel == 'gauss':
            psf_imsize = int( round(2*4*psf_sigma) )
            # MUST BE ODD!
            if (psf_imsize % 2) == 0:
                psf_imsize+=1
            psfdims = [psf_imsize,psf_imsize]

            if self.verbose:
                print("  dims:",psfdims)

            self.psf = mom2disk('gauss', Irr_psf, Irc_psf, Icc_psf, psfdims)
            self['psfcen'] = [(psfdims[0]-1)/2., (psfdims[1]-1)/2.]
        elif psfmodel == 'sdss':
            raise ValueError("You must call _make_sdss_psf for psfmodel='sdss'")
        else:
            raise ValueError("Unknown psf model: '%s'" % psfmodel)


    def make_image0(self):
        # run make_psf first!
        #
        # unconvolved model
        # even though it is unconvolved, make it as large as we 
        # would expect if it were convolved with the PSF
        # also make sure the psf image size is still smaller

        objmodel = self['objmodel']
        if objmodel == 'gauss':
            sigfac = 4.5
        elif objmodel == 'exp':
            #sigfac = 5.0
            sigfac = 7.0
        else:
            raise ValueError("Unknown obj model: '%s'" % objmodel)
        Irr=self['Irr']
        Irc=self['Irc']
        Icc=self['Icc']

        # expected size
        Irr_exp = Irr + self['Irr_psf']
        Icc_exp = Icc + self['Icc_psf']

        sigma = self.mom2sigma(Irr_exp+Icc_exp)
        dim = int( round(2.0*sigfac*sigma) )

        # might still be too small, just in case
        if dim < self.psf.shape[0]:
            dim = self.psf.shape[0]
        dims = [dim,dim]
        cen = (dim-1)/2.0
        cen = [cen,cen]

        self['dims'] = dims
        self['cen'] = cen

        if self.verbose:
            print("Making object model '%s'" % objmodel)
            print("  Irr: ",Irr)
            print("  Irc: ",Irc)
            print("  Icc: ",Icc)
            print("  dims:",dims)
            print("  cen: ",cen)

        self.image0 = mom2disk(objmodel, 
                               Irr, 
                               Irc, 
                               Icc, 
                               dims, 
                               cen=cen, 
                               counts=self['counts'])


    def make_image(self):

        # expand the image if needed to accomodate psf
        # if image0 is bigger than psf, it is returned unchanged

        if self.verbose:
            print("Convolving to final image")

        dims = self['dims']
        image0 = self.image0
        psf = self.psf

        # we don't have to do a numerical convolution if both are gaussians
        if self['psfmodel'] == 'gauss' and self['objmodel'] == 'gauss':
            Irr=self['Irr']+self['Irr_psf']
            Irc=self['Irc']+self['Irc_psf']
            Icc=self['Icc']+self['Icc_psf']
            image = mom2disk('gauss', 
                             Irr, 
                             Irc, 
                             Icc, 
                             self['dims'], 
                             cen=self['cen'],
                             counts=self['counts'])
        else:

            # this should be un-necessary
            image0_expand = images.expand(self.image0, self.psf.shape, verbose=self.verbose)

            image = scipy.signal.fftconvolve(image0_expand, self.psf, mode='same')

            if image.shape[0] > dims[0] or image.shape[1] > dims[1]:
                if self.verbose:
                    print("  Trimming back to requested size")
                image = image[ 0:self['dims'][0], 0:self['dims'][1] ]

        self.image = image

    def mom2sigma(self, T):
        return numpy.sqrt(T/2)

    def show(self):
        import biggles

        levels=7
        tab = biggles.Table(2,2)
        tab[0,0] = images.view(self.image0, show=False, levels=levels, min=0)
        tab[0,1] = images.view(self.psf, show=False, levels=levels, min=0)
        tab[1,0] = images.view(self.image, show=False, levels=levels, min=0)

        # cross-section across rows
        cen = self['cen']
        imrows = self.image[:, cen[1]]
        imcols = self.image[cen[0], :]

        crossplt = biggles.FramedPlot()
        hrows = biggles.Histogram(imrows, color='blue')
        hrows.label = 'Center rows'
        hcols = biggles.Histogram(imcols, color='red')
        hcols.label = 'Center columns'

        key = biggles.PlotKey(0.1, 0.9, [hrows, hcols])

        crossplt.add(hrows, hcols, key)
        crossplt.aspect_ratio=1
        yr = crossplt._limits1().yrange()
        yrange = (yr[0], yr[1]*1.2)
        crossplt.yrange = yrange

        tab[1,1] = crossplt

        tab.show()
