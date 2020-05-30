import os
def shift_image(image, shift):
    """
    Shift the image

    bi-linear interpolation is used

    parameters
    ----------
    shift: sequence[2]
        shift in row,column directions.  Can be non-integer

    outputs
    -------
    the shifted image

    dependencies
    ------------
    scipy.ndimage
    """
    import scipy.ndimage

    output = scipy.ndimage.interpolation.shift(image, shift, 
                                               output='f4',
                                               order=1,
                                               mode='constant',
                                               cval=0.0)

    return output

def get_jpg_dir(**keys):
    run=keys.get('run',None)
    camcol=keys.get('camcol',None)
    if run is None or camcol is None:
        raise ValueError("send run= and camcol=")
    d=os.path.expanduser('~/oh/boss-frames')
    d=os.path.join(d,str(run), str(camcol))
    return d

def get_jpg_path(**keys):
    d=get_jpg_dir(**keys)
    run=keys['run']
    camcol=keys['camcol']
    field=keys.get('field',None)
    if field is None:
        raise ValueError("send run=, camcol=, field=")
    fname='frame-irg-%06d-%s-%04d.jpg' % (run,camcol,field)
    path=os.path.join(d,fname)
    return path


def make_frame_jpg(**keys):
    from PIL import Image
    from numpy import array, flipud
    import sdsspy
    import esutil as eu
    import images
    #import biggles

    g,gh=sdsspy.read('frame', filter='g', header=True, **keys)
    r,rh=sdsspy.read('frame', filter='r', header=True, **keys)
    i,ih=sdsspy.read('frame', filter='i', header=True, **keys)

    gwcs=eu.wcsutil.WCS(gh)
    rwcs=eu.wcsutil.WCS(rh)
    iwcs=eu.wcsutil.WCS(ih)

    # convert r-band center to ra,dec
    rowcen_r=r.shape[0]/2.
    colcen_r=r.shape[1]/2.

    racen_r,deccen_r = rwcs.image2sky(colcen_r, rowcen_r)

    # now see where that falls in the g and i images
    gcol,grow=gwcs.sky2image(racen_r, deccen_r)
    icol,irow=iwcs.sky2image(racen_r, deccen_r)

    gshift=[rowcen_r-grow, colcen_r-gcol]
    ishift=[rowcen_r-irow, colcen_r-icol]

    g = shift_image(g, gshift)
    i = shift_image(i, ishift)

    rescale=0.4
    scales= rescale*array([5., 6.5, 9.])
    nonlinear=0.15

    colorim=images.get_color_image(i, r, g, 
                                   scales=scales,
                                   nonlinear=nonlinear,
                                   satval=30.)


    fname=get_jpg_path(**keys)
    print 'writing:',fname
    dname=os.path.dirname(fname)
    if not os.path.exists(dname):
        print 'making dirs:',dname
        os.makedirs(dname)

    #colorim = (colorim*255).astype('u1')
    colorim = images.bytescale(colorim)
    colorim=flipud(colorim)

    pim=Image.fromarray(colorim)
    pim.save(fname, quality=90)
