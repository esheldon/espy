from numpy import zeros, where, array

griscales = array([4.9,5.7,7.8])
#griscales /= griscales.min()

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


def get_id_info(**keys):
    import sdsspy
    out={}
    cat=keys.get('cat',None)
    if cat is not None:
        if len(cat['run'].shape)==0:
            tcat=cat
        else:
            tcat=cat[0]
        for k in ['run','camcol','field','id']:
            out[k] = tcat[k]
    elif 'objid' in keys:
        # for the title
        ids=sdsspy.util.objid_extract(keys['objid'])
        for k in ids:
            out[k] = ids[k]
    elif 'photoid' in keys:
        # for the title
        ids=sdsspy.util.photoid_extract(keys['photoid'])
        for k in ids:
            out[k] = ids[k]

    return out


def get_color_image(img, imr, imi, maxval=None,
                    alpha=0.02, nonlinearity=8.0):
    """
    Get the color image.  Note the images will already be
    transposed and scaled appropriately.
    """
    imlist = scale_image_list([img,imr,imi], 
                              scales=griscales,
                              maxval=maxval,
                              same_stretch=True, 
                              reverse=False, alpha=alpha, 
                              zero_clip=True,
                              nonlinearity=nonlinearity)

    img=imlist[0].transpose()
    imr=imlist[1].transpose()
    imi=imlist[2].transpose()
    st=imr.shape
    colorim=zeros( (st[0], st[1], 3) )

    colorim[:,:,0] = imi[:,:] # "r"
    colorim[:,:,1] = imr[:,:] # "g"
    colorim[:,:,2] = img[:,:] # "b"

    return colorim

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

    biggles.configure( 'default', 'fontsize_min', 1)

    band=keys.get('band',None)
    cat=keys.get('cat',None)
    show=keys.get('show',True)
    same_stretch=keys.get('same_stretch',True)
    reverse=keys.get('reverse',True)
    color=keys.get('color',False)

    alpha=0.02
    nonlinearity=8.0

    idinfo=get_id_info(**keys)

    keys['trim']=False
    #imdict=sdsspy.read('fpAtlas', **keys)
    imdict=sdsspy.read('fpAtlas', **idinfo)
    
    if not imdict:
        raise RuntimeError("no images returned")

    if band is not None:
        bnum=sdsspy.util.FILTERNUM[band]
        bchar=sdsspy.util.FILTERCHAR[band]

        idinfo['bchar']=bchar
        title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d-%(bchar)s' % idinfo
        im=imdict['images'][bnum]
        s=im.shape

        imslist=scale_image_list([im], same_stretch=same_stretch, 
                                 reverse=reverse, alpha=alpha, 
                                 nonlinearity=nonlinearity)

        ims=imslist[0].transpose()
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
        plt=biggles.FramedPlot()
        d = biggles.Density(ims, ranges)

        plt.add(d)
        plt.title=title
        plt.aspect_ratio=s[0]/float(s[1])
        labs=get_shaded_label(0.9,0.9,bchar)
        plt.add(*labs)
        if show:
            plt.show()
        return plt
    elif color:
        img=imdict['images'][1]
        imr=imdict['images'][2]
        imi=imdict['images'][3]
        s=imr.shape

        colorim=get_color_image(img, imr, imi)
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))

        d = biggles.Density(colorim, ranges)

        plt=biggles.FramedPlot()
        plt.add(d)
        plt.aspect_ratio=s[0]/float(s[1])
        title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % idinfo
        plt.title=title
        if show:
            plt.show()
        return plt


    imslist=scale_image_list(imdict['images'], same_stretch=same_stretch, 
                             reverse=reverse, alpha=alpha, 
                             nonlinearity=nonlinearity)

    nrow=2
    ncol=3
    plt=biggles.FramedArray(nrow,ncol)
    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % idinfo
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
    
def get_shaded_label(x, y, message, size=4, offset=0.003, color='red'):
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
    if ntot in (31,32,33,34,35,36):
        return (6,6)

    if ntot in (37,38,39,40,41,42):
        return (6,7)
    if ntot in (43,44,45,46,47,48,49):
        return (7,7)

    raise ValueError("don't know how to deal with ntot = %d" % ntot)

def view_family(cat=None, index=None, band=None, show=True, **keys):
    """
    If cat/index not sent the data will be read if the id info is
    sent
    """
    import sdsspy
    import biggles
    import images
    from sdsspy.atlas.atlas import NoAtlasImageError

    same_stretch=keys.get('same_stretch',True)
    reverse=keys.get('reverse',True)
    verbose=keys.get('verbose',False)
    color=keys.get('color',False)

    biggles.configure( 'default', 'fontsize_min', 0.75)
    if band is None and not color:
        raise ValueError("send band= or color=True")
    if cat is None or index is None:
        # try to read the data
        cat=sdsspy.read('photoObj', lower=True, **keys)
        if cat is None:
            raise ValueError("send cat= and index= or valid ids")
        idinfo=get_id_info(**keys)
        index,=where(cat['id'] == idinfo['id'])
        if index.size==0:
            raise ValueError("problem finding id in catalog")
        index=index[0]

    alpha=0.02
    nonlinearity=8.0
    if band is not None:
        bnum=sdsspy.util.FILTERNUM[band]
        bchar=sdsspy.util.FILTERCHAR[band]
    fam=sdsspy.util.Family(cat, index, verbose=verbose)

    children=fam.get_children()
    parent=fam.get_parent()
    grandparent=fam.get_grandparent()

    
    main_title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % cat[index]
    if band is not None:
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
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)
                #im=get_color_image(*imd['images'][1:1+3], alpha=alpha, 
                #                   nonlinearity=nonlinearity)


            title='grandparent %d' % cat['id'][grandparent[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    if parent.size > 0:
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][parent[0]])
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)
                #im=get_color_image(*imd['images'][1:1+3], alpha=alpha, 
                #                   nonlinearity=nonlinearity)


            title='parent %d' % cat['id'][parent[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    for i in xrange(children.size):
        try:
            imd=imdict=sdsspy.read('fpAtlas', run=run,camcol=camcol,field=field,
                                   rerun=rerun, id=cat['id'][children[i]])
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)
                #im=get_color_image(*imd['images'][1:1+3], alpha=alpha, 
                #                   nonlinearity=nonlinearity)


            title='child %d' % cat['id'][children[i]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    ntot = len(imlist)
    if ntot == 0:
        raise NoAtlasImageError("none of the objects had atlas images")

    if band is not None:
        imslist=scale_image_list(imlist, same_stretch=same_stretch, 
                                 reverse=reverse, alpha=alpha, 
                                 nonlinearity=nonlinearity)
        imslist=[im.transpose() for im in imslist]
    else:
        maxval=None
        minval=None
        for imd in imlist:
            iml=imd['images']
            img=iml[1]*griscales[0]
            imr=iml[2]*griscales[1]
            imi=iml[3]*griscales[2]
            tmaxval=array([img.max(), imr.max(), imi.max()]).max()
            tminval=array([img.min(), imr.min(), imi.min()]).min()
            if maxval is None or tmaxval > maxval:
                maxval=tmaxval
            if minval is None or tminval < minval:
                minval=tminval
        imslist=[]
        for imd in imlist:
            cim=get_color_image(*imd['images'][1:1+3], alpha=alpha, 
                                nonlinearity=nonlinearity,
                                maxval=maxval)
            imslist.append(cim)

    if ntot==1:
        im=imslist[0]
        s=im.shape
        st=imt.shape

        plt=biggles.FramedPlot()
        ranges = ((-0.5, -0.5), (st[0]-0.5, st[1]-0.5))
        d = biggles.Density(im, ranges)
        plt.add(d)
        plt.title = main_title + ' ' + titles[0]
        plt.aspect_ratio = st[1]/float(st[0])
    else:

        nrow,ncol=get_family_grid(ntot)
        plt=biggles.FramedArray(nrow,ncol)
        s=imslist[0].shape
        plt.aspect_ratio = s[1]/float(s[0])*nrow/float(ncol)
        plt.title=main_title
        for i in xrange(ntot):
            im=imslist[i]
            row=i/ncol
            col=i % ncol

            s=im.shape
            ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
            d = biggles.Density(im, ranges)

            plt[row,col].add(d)

            if color:
                lcolor='gray75'
            else:
                lcolor='black'
            lab=biggles.PlotLabel(0.075,0.9,titles[i],halign='left',
                                  color=lcolor)
            plt[row,col].add(lab)

    if show:
        plt.show()

    return plt

def scale_image_list(imlist, alpha=0.02, nonlinearity=8.0, 
                     same_stretch=True, 
                     zero_clip=False,
                     maxval=None,
                     reverse=True,
                     scales=None):
    """
    Use zero_clip=True for color images for sure.

    Universal maxval would be good if we were making images for the whole survey

    Can use something like
        maxval=2**16-1
    """
    import images
    minv=None
    maxv=None
    imslist=[]
    for i,im in enumerate(imlist):
        #print im.max()-1000
        ims=im.astype('f4')-1000
        if zero_clip:
            ims=ims.clip(0.0, ims.max())
        if scales is not None:
            ims *= scales[i]
        ims=images.asinh_scale(ims, alpha=alpha, nonlinearity=nonlinearity)

        min0=ims.min()
        max0=ims.max()
        if same_stretch:
            if minv is None or min0 < minv:
                minv=min0
            if maxv is None or max0 > maxv:
                maxv=max0
        else:
            ims -= min0
            ims /= (max0-min0)

        imslist.append(ims)

    if same_stretch or reverse:
        for i in xrange(len(imslist)):
            ims=imslist[i]

            if same_stretch:
                if maxval is not None:
                    maxv=maxval-1000.0
                    maxv=images.asinh_scale(maxv, alpha=alpha, nonlinearity=nonlinearity)
                ims -= minv
                ims /= (maxv-minv)

            if reverse:
                ims *=-1
                ims += 1.0
            imslist[i] = ims

    return imslist
