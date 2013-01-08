def view_atlas(**keys):
    """
    View atlas images for the input object

    parameters
    ----------
    cat: optional
        A full structure with run,rerun,camcol,field,id
    run: optional
        Run id
    rerun: optional
        rerun id
    camcol: optional
        camcol id
    field: opt
        field id
    id: opt
        id in field
    band: optional
        Only display this band
    """
    import images
    import biggles
    import sdsspy

    band=keys.get('band',None)
    cat=keys.get('cat',None)
    show=keys.get('show',True)

    if cat is not None:
        if len(cat['run'].shape)==0:
            tcat=cat
        else:
            tcat=cat[0]
        for k in ['run','camcol','field','id']:
            keys[k] = tcat[k]

    imdict=sdsspy.read('fpAtlas', **keys)
    
    if not imdict:
        raise RuntimeError("no images returned")

    if band is not None:
        bnum=sdsspy.util.FILTERNUM[band]
        bchar=sdsspy.util.FILTERCHAR[band]

        keys['bchar']=bchar
        title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d-%(bchar)s' % keys
        im=imdict['images'][bnum]
        multi=keys.get('multiview',False)
        if multi:
            plt=images.multiview(im,title=title,show=show)
        else:
            plt=images.view(im,title=title,show=show)
        return plt

    
    nrow=2
    ncol=3


    # force common range
    imr=imdict['images'][2].astype('f4').transpose()-1000
    immin=imr.min()
    immax=imr.max()

    s=imr.shape

    doarr=True
    if doarr:
        plt=biggles.FramedArray(2,3)
        plt.aspect_ratio = s[1]/float(s[0])*2./3.
    else:
        plt=biggles.Table(2,3)

    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % keys
    plt.title=title
    for i in xrange(5):
        bchar=sdsspy.util.FILTERCHAR[i]
        row=i/ncol
        col=i % ncol

        if doarr:
            im0=imdict['images'][i].astype('f4')-1000

            im =images.scale(im0, min=immin, max=immax, reverse=True)
            im = im.transpose()
            s=im.shape
            ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))

            d = biggles.Density(im, ranges)

            plti = biggles.Plot()
            plti.add(d)

            lab_bright=biggles.PlotLabel(0.9,0.9,bchar,halign='right',
                                         color='yellow',size=4)
            lab_dark=biggles.PlotLabel(0.9+0.005,0.9-0.005,bchar,halign='right',
                                       color='black',size=4)
            plt[row,col].add(d,lab_dark,lab_bright)

        else:
            im=imdict['images'][i].astype('f4')-1000
            plti=images.view(im, min=immin, max=immax, show=False)
            plt[row,col]=plti

    if show:
        plt.show()
    return plt 
    
