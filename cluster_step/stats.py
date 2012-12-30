from numpy import zeros
from esutil.stat import wmom, histogram

SHAPE_NOISE=0.32/2


def bin_shear_data(data, bin_field, nperbin):
    h,rev=histogram(data[bin_field], nperbin=nperbin, rev=True)

    nbin=len(h)
    dt,fields=get_binned_dtype(bin_field)
    bindata=zeros(nbin, dtype=dt)

    for i in xrange(nbin):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]

            err2=data['gcov'][w,0,0] + data['gcov'][w,1,1]
            wts = 1.0/(SHAPE_NOISE**2 + err2)

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
                bindata[field][i] = wmean
                bindata[err_field][i] = werr

    bindata['g1'] /= bindata['g1sens']
    bindata['g1_err'] /= bindata['g1sens']
    bindata['g2'] /= bindata['g2sens']
    bindata['g2_err'] /= bindata['g2sens']


    return bindata

def get_binned_dtype(bin_field):
    fields=[bin_field,'g1','g2','g1sens','g2sens']
    dt=[]
    for f in fields:
        dt += [(f, 'f8'),(f+'_err','f8')]
    return dt, fields


