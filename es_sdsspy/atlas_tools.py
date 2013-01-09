from numpy import zeros, where

def get_std(im):
    from esutil.stat import sigma_clip
    mn,std = sigma_clip(im)
    return std

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

    alpha=0.02
    nonlinearity=8.0

    if cat is not None:
        if len(cat['run'].shape)==0:
            tcat=cat
        else:
            tcat=cat[0]
        for k in ['run','camcol','field','id']:
            keys[k] = tcat[k]

    keys['trim']=False
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
            plt=images.multiview(im,title=title,show=False,
                                 alpha=alpha, nonlinearity=nonlinearity)
        else:
            plt=images.view(im,title=title,show=False,
                            alpha=alpha, nonlinearity=nonlinearity)

        labs=get_shaded_label(0.9,0.9,bchar)
        plt.add(*labs)
        if show:
            plt.show()
        return plt
    
    #s1=imdict['images'][2].transpose().shape
    s1=imdict['images'][2].shape

    imtot=zeros((s1[0], 5*s1[1]), dtype=imdict['images'][2].dtype)
    #imtot=zeros((s1[0], 5*s1[1]), dtype='f4')
    #imtot[:,:] = 1000

    for i in xrange(5):
        #im=imdict['images'][2]
        #im=imdict['images'][i].astype('f4')-1000.0
        im=imdict['images'][i]
        i1=i*s1[1]
        i2=(i+1)*s1[1]
        imtot[:,i1:i2] = im

    #w1000=where(imtot == 1000)

    #imtot = imtot.astype('f4')-1000

    #imtot = images.scale(imtot, alpha=alpha, nonlinearity=nonlinearity, reverse=True)

    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % keys
    plt=images.view(imtot,title=title,show=False,
                    alpha=alpha, nonlinearity=nonlinearity)

    """
    imtot_trans=imtot.transpose()
    s=imtot_trans.shape
    ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
    d = biggles.Density(imtot_trans, ranges)

    plt=biggles.FramedPlot()
    plt.add(d)
    plt.aspect_ratio = s[1]/float(s[0])
    """
    """
    nrow=2
    ncol=3


    # force common range
    immin=imr.min()
    immax=imr.max()

    s=imr.shape

    plt=biggles.FramedArray(2,3)
    plt.aspect_ratio = s[1]/float(s[0])*2./3.

    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % keys
    plt.title=title
    for i in xrange(5):
        bchar=sdsspy.util.FILTERCHAR[i]
        row=i/ncol
        col=i % ncol

        im0=imdict['images'][i].astype('f4')-1000

        im =images.scale(im0, min=immin, max=immax, reverse=True)
        im = im.transpose()
        s=im.shape
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))

        d = biggles.Density(im, ranges)

        plti = biggles.Plot()
        plti.add(d)

        labs=get_shaded_label(0.9,0.9,bchar)
        plt[row,col].add(d)
        plt[row,col].add(*labs)
    """

    if show:
        plt.show()
    return plt 
    
def get_shaded_label(x, y, message, size=4, offset=0.003, color='yellow'):
    import biggles
    lab_bright=biggles.PlotLabel(x, y,
                                 message,
                                 halign='right',
                                 color=color,
                                 size=size)
    lab_dark=biggles.PlotLabel(x+offset,y-offset,
                               message,
                               halign='right',
                               color='black',
                               size=size)

    return lab_dark, lab_bright

def get_family_grid(ntot):
    if ntot==1:
        raise ValueError("expect ntot > 1")

    if ntot==2:
        return (1,2)

    if ntot in (3,4):
        return (2,2)

    if ntot in (5,6):
        return (2,3)

    if ntot in (7,8,9):
        return (3,3)

    if ntot in (10,11,12):
        return (3,4)

    if ntot in (13,14,15,16):
        return (4,4)

    if ntot in (17,18,19,20):
        return (4,5)

    if ntot in (21,22,23,24,25):
        return (5,5)

    if ntot in (26,27,28,29,30):
        return (5,6)
    elif ntot in (31,32,33,34,35,36):
        return (6,6)

    raise ValueError("don't know how to deal with ntot = %d" % ntot)


def view_family(cat=None, index=None, show=True, **keys):
    import sdsspy
    import biggles
    from sdsspy.atlas.atlas import NoAtlasImageError

    biggles.configure( 'default', 'fontsize_min', 1)
    if cat is None or index is None:
        raise ValueError("send cat= and index=")

    fam=sdsspy.util.Family(cat, index)

    children=fam.get_children()
    parent=fam.get_parent()
    bright=fam.get_bright()

    
    main_title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % cat[index]

    plots=[]

    itot=0

    if bright.size > 0:
        try:
            plti=view_atlas(cat=cat[bright[0]], show=False, **keys)
            plti.title='bright %d' % cat['id'][bright[0]]
            plots.append(plti)
        except NoAtlasImageError:
            pass

    if parent.size > 0:
        try:
            plti=view_atlas(cat=cat[parent[0]], show=False, **keys)
            plti.title='parent %d' % cat['id'][parent[0]]
            plots.append(plti)
        except NoAtlasImageError:
            pass

    for i in xrange(children.size):
        try:
            plti=view_atlas(cat=cat[children[i]], show=False, **keys)
            plti.title='child %d' % cat['id'][children[i]]
            plots.append(plti)
        except NoAtlasImageError:
            pass

    ntot = len(plots)
    if ntot == 0:
        raise RuntimeError("none of the objects had atlas images")

    if ntot==1:
        plt=plots[0]
        plt.title = main_title + ' ' + plt.title
    else:
        nrow,ncol=get_family_grid(ntot)
        plt=biggles.Table(nrow,ncol)
        plt.title=main_title
        plt.cellpadding=0
        plt.cellspacing=0
        for i in xrange(len(plots)):
            row=i/ncol
            col=i % ncol

            plt[row,col]=plots[i]

    if show:
        plt.show()

    return plt

def view_family_test(cat=None, index=None, band=None, show=True, **keys):
    import sdsspy
    import biggles
    import images
    from sdsspy.atlas.atlas import NoAtlasImageError

    biggles.configure( 'default', 'fontsize_min', 1)
    if cat is None or index is None or band is None:
        raise ValueError("send cat= and index= and band=")

    alpha=0.02
    nonlinearity=8.0
    bnum=sdsspy.util.FILTERNUM[band]
    bchar=sdsspy.util.FILTERCHAR[band]
    fam=sdsspy.util.Family(cat, index)

    children=fam.get_children()
    parent=fam.get_parent()
    bright=fam.get_bright()

    
    main_title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % cat[index]
    main_title += ' %s' % bchar

    run=cat['run'][index]
    camcol=cat['camcol'][index]
    field=cat['field'][index]
    rerun=cat['rerun'][index]

    imlist=[]
    titles=[]
    if bright.size > 0:
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][bright[0]])
            im=imd['images'][bnum]
            imlist.append(im)

            title='bright %d' % cat['id'][bright[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    if parent.size > 0:
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][parent[0]])
            im=imd['images'][bnum]
            imlist.append(im)

            title='parent %d' % cat['id'][parent[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    for i in xrange(children.size):
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][children[i]])
            im=imd['images'][bnum]
            imlist.append(im)

            title='child %d' % cat['id'][children[i]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    ntot = len(imlist)
    if ntot == 0:
        raise RuntimeError("none of the objects had atlas images")

    if ntot==1:
        im=imlist[0]
        plt=images.view(im,show=False,
                        alpha=alpha, nonlinearity=nonlinearity)
        plt.title = main_title + ' ' + titles[0]
    else:

        minval=9999.0e9
        maxval=-9999.0e9
        imslist=[]
        for im in imlist:
            ims=images.asinh_scale(im, alpha=alpha, nonlinearity=nonlinearity, dtype='f4')
            min0=ims.min()
            max0=ims.max()
            if min0 < minval:
                minval=min0
            if max0 > maxval:
                maxval=max0

            imslist.append(ims)

        nrow,ncol=get_family_grid(ntot)
        plt=biggles.FramedArray(nrow,ncol)
        s=imlist[0].transpose().shape

        plt.aspect_ratio = s[1]/float(s[0])*ncol/float(nrow)
        plt.title=main_title
        for i in xrange(ntot):
            im=imslist[i]
            row=i/ncol
            col=i % ncol

            if True:
                ims = im-minval
                ims /= (maxval-minval)
            else:
                ims = im-im.min()
                ims /= ims.max()
            ims = 1.0-ims

            ims=ims.transpose()
            s=ims.shape
            ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
            d = biggles.Density(ims, ranges)

            plt[row,col].add(d)

            lab=biggles.PlotLabel(0.05,0.95,titles[i],halign='left')
            plt[row,col].add(lab)
            #plt[row,col].title=titles[i]

    """
        nrow,ncol=get_family_grid(ntot)
        plt=biggles.FramedArray(nrow,ncol)
        plt.title=main_title
        plt.cellpadding=0
        plt.cellspacing=0
        for i in xrange(len(plots)):
            row=i/ncol
            col=i % ncol

            plt[row,col]=plots[i]
    """
    if show:
        plt.show()

    return plt
