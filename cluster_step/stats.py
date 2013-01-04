from numpy import zeros, median
from esutil.stat import wmom, wmedian, histogram

SHAPE_NOISE=0.32/2


def bin_shear_data(data, bin_field, nperbin, use_median=False, use_wmedian=False):
    """
    median or wmedian are not an improvement
    """
    h,rev=histogram(data[bin_field], nperbin=nperbin, rev=True)

    nbin=len(h)
    dt,fields=get_binned_dtype(bin_field)
    bindata=zeros(nbin, dtype=dt)

    for i in xrange(nbin):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]

            wts=get_weights(data['gcov'], ind=w)

            for field in fields:
                if field == 'g1':
                    fdata=data['g'][w,0]
                elif field=='g2':
                    fdata=data['g'][w,1]
                elif field=='g1sens':
                    fdata=data['gsens'][w,0]
                elif field=='g2sens':
                    fdata=data['gsens'][w,1]
                else:
                    fdata=data[field][w]

                err_field=field+'_err'
                wmean,werr=wmom(fdata, wts, calcerr=True)
                if use_median:
                    bindata[field][i] = median(fdata)
                elif use_wmedian:
                    bindata[field][i] = wmedian(fdata,wts)
                else:
                    bindata[field][i] = wmean
                bindata[err_field][i] = werr

    bindata['g1'] /= bindata['g1sens']
    bindata['g1_err'] /= bindata['g1sens']
    bindata['g2'] /= bindata['g2sens']
    bindata['g2_err'] /= bindata['g2sens']

    return bindata

def get_weights(gcov, ind=None):
    if ind is not None:
        err2=gcov[ind,0,0] + gcov[ind,1,1]
    else:
        err2=gcov[:,0,0] + gcov[:,1,1]

    # get average of the errors, since the
    # shape noise is per component
    err2 /= 2
    wts = 1.0/(SHAPE_NOISE**2 + err2)
    return wts

def get_mean_shear(g, gsens, weight):
    g_mean,g_err=wmom(g, weight, calcerr=True)
    gsens_mean,gsens_err=wmom(gsens, weight, calcerr=True)

    g_mean /= gsens_mean
    g_err /= gsens_mean

    return {'g':g_mean, 
            'g_err':g_err, 
            'gsens':gsens_mean, 
            'gsens_err':gsens_err}

def get_mean_shear_shstruct(data):
    """
    Input is the full output struct from the shear code
    """
    wts=get_weights(data['gcov'])
    g1res=get_mean_shear(data['g'][:,0], data['gsens'][:,0], wts)
    g2res=get_mean_shear(data['g'][:,1], data['gsens'][:,1], wts)

    out={}
    out['g1']=g1res['g']
    out['g1_err']=g1res['g_err']
    out['g1sens']=g1res['gsens']
    out['g1sens_err']=g1res['gsens_err']

    out['g2']=g2res['g']
    out['g2_err']=g2res['g_err']
    out['g2sens']=g2res['gsens']
    out['g2sens_err']=g2res['gsens_err']

    return out


def get_binned_dtype(bin_field):
    fields=[bin_field,'g1','g2','g1sens','g2sens']
    dt=[]
    for f in fields:
        dt += [(f, 'f8'),(f+'_err','f8')]
    return dt, fields


