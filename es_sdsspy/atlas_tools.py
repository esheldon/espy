import os
from numpy import zeros, where, array, arcsinh, sqrt

# this is for scaled color
_default_nonlinear=0.15

# from sdss_frame_jpg (u and z are made up though)
_scales = 0.4*array([9.0,9.0,6.5,5.0,5.0])

# crude calibration in case calib files are not available
_nmgypercount=array([0.007897,
                     0.003552,
                     0.004577,
                     0.006181,
                     0.035935])

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

    shift:
        Apply sub-pixel astronometric shifts, default False
    """
    import images
    import biggles
    import sdsspy

    biggles.configure( 'default', 'fontsize_min', 1)

    band=keys.get('band',None)
    cat=keys.get('cat',None)
    show=keys.get('show',True)
    color=keys.get('color',False)
    shift=keys.get('shift',False)

    idinfo=get_id_info(**keys)

    imdict=read_atlas(**idinfo)

    nonlinear=keys.get('nonlinear', _default_nonlinear)

    if band is not None:
        bnum=sdsspy.util.FILTERNUM[band]
        bchar=sdsspy.util.FILTERCHAR[band]

        idinfo['bchar']=bchar
        title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d-%(bchar)s' % idinfo
        im=imdict['images'][bnum]

        plt=images.view(im, nonlinear=nonlinear,show=False)

        plt.title=title
        lab=biggles.PlotLabel(0.9,0.9,bchar,color='yellow')
        plt.add(lab)
        if show:
            plt.show()
        return plt
    elif color:
        img=imdict['images'][1].transpose()
        imr=imdict['images'][2].transpose()
        imi=imdict['images'][3].transpose()

        if shift:
            img,imi=shift_images(imdict, img, imi)

        colorim=images.get_color_image(imi,imr,img,
                                       satval=30.,
                                       nonlinear=nonlinear)

        s=imdict['images'][0].shape
        ranges = ((-0.5, -0.5), (s[1]-0.5, s[0]-0.5))

        d = biggles.Density(colorim, ranges)

        plt=biggles.FramedPlot()
        plt.add(d)
        plt.aspect_ratio=s[0]/float(s[1])
        title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % idinfo
        plt.title=title
        if show:
            plt.show()
        return plt


    imslist=scale_image_list(imdict['images'], nonlinear)

    nrow=2
    ncol=3
    plt=biggles.FramedArray(nrow,ncol)
    title='%(run)06d-%(camcol)d-%(field)04d-%(id)05d' % idinfo
    plt.title=title

    s=imslist[0].shape
    plt.aspect_ratio = s[0]/float(s[1])*2./3.

    for i in xrange(5):
        bchar=sdsspy.util.FILTERCHAR[i]
        row=i/ncol
        col=i % ncol

        im=imslist[i]
        s=im.shape

        im = im.transpose()
        ranges = ((-0.5, -0.5), (s[1]-0.5, s[0]-0.5))

        d = biggles.Density(im, ranges)

        #labs=get_shaded_label(0.9,0.9,bchar)
        lab=biggles.PlotLabel(0.9,0.9,bchar,color='yellow')
        plt[row,col].add(d)
        plt[row,col].add(lab)

    if show:
        plt.show()
    return plt 
    

def view_family(**keys):
    """
    If cat/index not sent the data will be read if the id info is
    sent
    """
    import sdsspy
    import biggles
    import images
    from sdsspy.atlas.atlas import NoAtlasImageError

    band=keys.get('band',None)
    color=keys.get('color',False)
    show=keys.get('show',True)
    nonlinear=keys.get('nonlinear', _default_nonlinear)
    minval=keys.get('min',1.0)

    biggles.configure( 'default', 'fontsize_min', 0.75)
    if band is None and not color:
        raise ValueError("send band= or color=True")

    fam,cat,index=read_family(**keys)

    if band is not None:
        bnum=sdsspy.util.FILTERNUM[band]
        bchar=sdsspy.util.FILTERCHAR[band]

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
            imd=read_atlas(run=run,camcol=camcol,field=field,
                           rerun=rerun, id=cat['id'][grandparent[0]])
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)

            title='grandparent %d' % cat['id'][grandparent[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    if parent.size > 0:
        try:
            imd=read_atlas(run=run,camcol=camcol,field=field,
                           rerun=rerun, id=cat['id'][parent[0]])
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)


            title='parent %d' % cat['id'][parent[0]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    for i in xrange(children.size):
        try:
            imd=read_atlas(run=run,camcol=camcol,field=field,
                           rerun=rerun, id=cat['id'][children[i]])
            if band is not None:
                im=imd['images'][bnum]
                imlist.append(im)
            else:
                imlist.append(imd)

            title='child %d' % cat['id'][children[i]]
            titles.append(title)
        except NoAtlasImageError:
            pass

    ntot = len(imlist)
    if ntot == 0:
        raise NoAtlasImageError("none of the objects had atlas images")

    if band is not None:
        imslist=scale_image_list(imlist, nonlinear)
        imslist=[im.transpose() for im in imslist]
    else:
        imslist=[]
        for imd in imlist:
            img=imd['images'][1].transpose()
            imr=imd['images'][2].transpose()
            imi=imd['images'][3].transpose()

            cim=images.get_color_image(imi,imr,img,
                                       satval=30.,
                                       nonlinear=nonlinear)

            imslist.append(cim)

    nrow,ncol=get_family_grid(ntot)
    plt=biggles.FramedArray(nrow,ncol)
    s=imslist[0].shape
    plt.aspect_ratio = s[1]/float(s[0])*nrow/float(ncol)
    plt.title=main_title
    lcolor='gray75'
    for i in xrange(ntot):
        im=imslist[i]
        row=i/ncol
        col=i % ncol

        # transposed shape
        s=im.shape
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
        d = biggles.Density(im, ranges)

        plt[row,col].add(d)

        lab=biggles.PlotLabel(0.075,0.9,titles[i],halign='left',
                              color=lcolor)
        plt[row,col].add(lab)

    if show:
        plt.show()

    return plt


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

def view_many_family(cat=None, band=None, color=False, **keys):
    from sdsspy.atlas.atlas import NoAtlasImageError

    if cat is None or (band is None and color is None):
        raise ValueError("send cat= and band= or color=")

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
            view_family(cat=cat, index=ii, band=band, color=color, **keys)
        except NoAtlasImageError as e:
            print e

        if prompt:
            key=raw_input('hit a key (q to quit): ')
            if key.lower() == 'q':
                return


def get_family_grid(ntot):
    sq=int(sqrt(ntot))
    if ntot==sq*sq:
        return (sq,sq)
    elif ntot <= sq*(sq+1):
        return (sq,sq+1)
    else:
        return (sq+1,sq+1)



def scale_image_list(imlist, nonlinear):
    import images
    imslist=[]
    for im in imlist:
        ims=images.asinh_scale(im, nonlinear)
        imslist.append(ims)

    return imslist

def get_test_data():
    from numpy import log10
    import sdsspy
    from . import select
    t=sdsspy.read('calibobj.gal', 756, 3, 300, rerun='301', lower=True)
    sel=select.Selector(t)
    mag=22.5-2.5*log10( t['modelflux'][:,2].clip(0.001, t['modelflux'][:,2].max()))
    fl=sel.flag_logic()

    w,=where(fl & (mag < 15) & (mag > 14))

    return t,w

def read_atlas(**keys):
    import sdsspy
    imd=sdsspy.read('fpAtlas', **keys)
    calib=sdsspy.read('calibPhotom', lower=True, **keys)
    if calib is None:
        print 'calibPhotom not found: using default rouch calib'
        nmgypercount=_nmgypercount
    else:
        w,=where(calib['field']==keys['field'])
        if w.size==0:
            raise ValueError("calib doesn't have field %d" % keys['field'])
        nmgypercount=calib['nmgypercount'][w[0],:]


    if not imd:
        raise RuntimeError("no images returned")

    sky=imd['SOFT_BIAS']

    for i in xrange(5):
        im = imd['images'][i].astype('f4')-sky

        im *= nmgypercount[i]*_scales[i]

        imd['images'][i] = im

    return imd

def read_family(**keys):
    import sdsspy
    # try to read the data
    cat=sdsspy.read('photoObj', lower=True, **keys)
    if cat is None:
        raise ValueError("Could not read photoObj for your object")
    idinfo=get_id_info(**keys)
    index,=where(cat['id'] == idinfo['id'])
    if index.size==0:
        raise ValueError("problem finding id in catalog")
    index=index[0]

    fam=sdsspy.util.Family(cat, index)
    return fam,cat,index

def get_id_info(**keys):
    import sdsspy
    import copy
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
    elif 'run' in keys:
        out['run']=keys['run']
        if 'rerun' in keys:
            out['rerun']=keys['rerun']
        if 'camcol' in keys:
            out['camcol']=keys['camcol']
        if 'field' in keys:
            out['field']=keys['field']
        if 'id' in keys:
            out['id']=keys['id']

    return out


def shift_images(imd, img, imi):
    """
    This isn't so useful yet because color terms
    are probably needed
    """
    import sdsspy
    import scipy.ndimage

    if False:
        return img, imi

    run=imd['run']
    camcol=imd['camcol']
    path=sdsspy.filename('photoField', run=run, camcol=camcol)

    if not os.path.exists(path):
        print 'astrometry not found, no shifts applied'
        return img, imi

    ast=sdsspy.astrom.Astrom(run=run,camcol=camcol)

    idrowg,idrowr,idrowi=imd['drow'][1:1+3]
    idcolg,idcolr,idcoli=imd['dcol'][1:1+3]

    # use corner of atlas image as reference
    rowr,colr= imd['row0'][2], imd['col0'][2]
    ra_r,dec_r = ast.pix2eq(imd['field'],'r',rowr, colr)

    rowg,colg= ast.eq2pix(imd['field'], 'g', ra_r, dec_r)
    rowi,coli= ast.eq2pix(imd['field'], 'i', ra_r, dec_r)


    gshift=array([rowg-rowr, colg-colr], dtype='f4')
    ishift=array([rowi-rowr, coli-colr], dtype='f4')

    gshift[0] -= idrowg
    gshift[1] -= idcolg
    ishift[0] -= idrowi
    ishift[1] -= idcoli

    gout = scipy.ndimage.interpolation.shift(img, gshift, 
                                             output='f4',
                                             order=1,
                                             mode='constant',
                                             cval=0.0)

    iout = scipy.ndimage.interpolation.shift(imi, ishift, 
                                             output='f4',
                                             order=1,
                                             mode='constant',
                                             cval=0.0)


    return gout,iout
