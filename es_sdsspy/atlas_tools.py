from numpy import zeros, where, array

def get_std(im):
    from esutil.stat import sigma_clip
    mn,std = sigma_clip(im)
    return std

def view_many_atlas(cat=None, **keys):
    from sdsspy.atlas.atlas import NoAtlasImageError

    if cat is None:
        raise ValueError("send cat=")

    prompt=keys.get('prompt',True)
    for i in xrange(cat.size):
        try:
            view_atlas(cat=cat[i], **keys)
        except NoAtlasImageError as e:
            print e

        if prompt:
            key=raw_input('hit a key (q to quit): ')
            if key.lower() == 'q':
                return

def view_many_family(cat=None, band=None, **keys):
    from sdsspy.atlas.atlas import NoAtlasImageError

    if cat is None or band is None:
        raise ValueError("send cat= and band=")

    prompt=keys.get('prompt',True)
    indices=keys.get('indices',None)

    if indices is not None:
        indices=array(indices, copy=False)
        ntot=indices.size
    else:
        ntot=cat.size

    for i in xrange(ntot):
        if indices is not None:
            ii=indices[i]
        else:
            ii=i
        try:
            view_family(cat=cat, index=ii, band=band, **keys)
        except NoAtlasImageError as e:
            print e

        if prompt:
            key=raw_input('hit a key (q to quit): ')
            if key.lower() == 'q':
                return


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
    same_stretch=keys.get('same_stretch',True)
    reverse=keys.get('reverse',True)

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
        im=imdict['images'][bnum].astype('f4')-1000
        multi=keys.get('multiview',False)
        if multi:
            plt=images.multiview(im,title=title,show=False,
                                 reverse=reverse,
                                 alpha=alpha, nonlinearity=nonlinearity)
        else:
            plt=images.view(im,title=title,show=False,
                            reverse=reverse,
                            alpha=alpha, nonlinearity=nonlinearity)

        labs=get_shaded_label(0.9,0.9,bchar)
        plt.add(*labs)
        if show:
            plt.show()
        return plt
    
    imslist=scale_image_list(imdict['images'], same_stretch=same_stretch, 
                             reverse=reverse, alpha=alpha, 
                             nonlinearity=nonlinearity)

    nrow=2
    ncol=3
    plt=biggles.FramedArray(nrow,ncol)
    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % keys
    plt.title=title

    s=imslist[0].transpose().shape
    plt.aspect_ratio = s[1]/float(s[0])*2./3.

    for i in xrange(5):
        bchar=sdsspy.util.FILTERCHAR[i]
        row=i/ncol
        col=i % ncol

        im=imslist[i]
        im = im.transpose()
        s=im.shape
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))

        d = biggles.Density(im, ranges)

        labs=get_shaded_label(0.9,0.9,bchar)
        plt[row,col].add(d)
        plt[row,col].add(*labs)

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

def view_family(cat=None, index=None, band=None, show=True, **keys):
    import sdsspy
    import biggles
    import images
    from sdsspy.atlas.atlas import NoAtlasImageError

    same_stretch=keys.get('same_stretch',True)
    reverse=keys.get('reverse',True)
    verbose=keys.get('verbose',False)

    biggles.configure( 'default', 'fontsize_min', 1)
    if cat is None or index is None or band is None:
        raise ValueError("send cat= and index= and band=")

    alpha=0.02
    nonlinearity=8.0
    bnum=sdsspy.util.FILTERNUM[band]
    bchar=sdsspy.util.FILTERCHAR[band]
    fam=sdsspy.util.Family(cat, index, verbose=verbose)

    children=fam.get_children()
    parent=fam.get_parent()
    grandparent=fam.get_grandparent()

    
    main_title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % cat[index]
    main_title += ' %s' % bchar

    run=cat['run'][index]
    camcol=cat['camcol'][index]
    field=cat['field'][index]
    rerun=cat['rerun'][index]

    imlist=[]
    titles=[]
    if grandparent.size > 0:
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][grandparent[0]])
            im=imd['images'][bnum]
            imlist.append(im)

            title='grandparent %d' % cat['id'][grandparent[0]]
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
        raise NoAtlasImageError("none of the objects had atlas images")

    imslist=scale_image_list(imlist, same_stretch=same_stretch, 
                             reverse=reverse, alpha=alpha, 
                             nonlinearity=nonlinearity)

    if ntot==1:
        im=imslist[0]
        imt=im.transpose()
        s=im.shape
        st=imt.shape

        plt=biggles.FramedPlot()
        ranges = ((-0.5, -0.5), (st[0]-0.5, st[1]-0.5))
        d = biggles.Density(im, ranges)
        plt.add(d)
        plt.title = main_title + ' ' + titles[0]
        plt.aspect_ratio = st[0]/float(st[1])
    else:

        nrow,ncol=get_family_grid(ntot)
        plt=biggles.FramedArray(nrow,ncol)
        #s=imslist[0].transpose().shape
        s=imslist[0].shape
        plt.aspect_ratio = s[0]/float(s[1])*nrow/float(ncol)
        plt.title=main_title
        for i in xrange(ntot):
            im=imslist[i]
            row=i/ncol
            col=i % ncol

            im=im.transpose()
            s=im.shape
            ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
            d = biggles.Density(im, ranges)

            plt[row,col].add(d)

            lab=biggles.PlotLabel(0.05,0.95,titles[i],halign='left')
            plt[row,col].add(lab)

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

def scale_image_list(imlist, alpha=0.02, nonlinearity=8.0, 
                     same_stretch=True, reverse=True):
    import images
    minval=None
    maxval=None
    imslist=[]
    for im in imlist:
        #print im.max()-1000
        ims=im.astype('f4')-1000
        ims=images.asinh_scale(ims, alpha=alpha, nonlinearity=nonlinearity)

        min0=ims.min()
        max0=ims.max()
        if same_stretch:
            if minval is None or min0 < minval:
                minval=min0
            if maxval is None or max0 > maxval:
                maxval=max0
        else:
            ims -= min0
            ims /= (max0-min0)

        imslist.append(ims)

    if same_stretch or reverse:
        for i in xrange(len(imslist)):
            ims=imslist[i]

            if same_stretch:
                maxval=2**16-1-1000
                maxval=images.asinh_scale(maxval, alpha=alpha, nonlinearity=nonlinearity)
                ims -= minval
                ims /= (maxval-minval)

            if reverse:
                ims *=-1
                ims += 1.0
            imslist[i] = ims

    return imslist
