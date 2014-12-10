"""
code do bin lens outputs
"""
from __future__ import print_function
import numpy
from numpy import where, log10
from .files import *
from . import averaging

def bin_run(run, bin_conf_name, jackreg_col=None):
    """
    Do the binning and write out a file
    """
    import fitsio

    data = read_collated(run)

    binner=Binner(bin_conf_name)
    res = binner.bin(data, jackreg_col=jackreg_col)

    if jackreg_col is None:
        fname=get_binned_file(run, bin_conf_name)
        d=get_binned_dir(run, bin_conf_name)
    else:
        fname=get_jack_file(run, bin_conf_name)
        d=get_jack_dir(run, bin_conf_name)

    if not os.path.exists(d):
        print("Making dir:",d)
        os.makedirs(d)
    
    print("writing:",fname)
    fitsio.write(fname, res, clobber=True)

class Binner(dict):
    """
    bin any
    """

    def __init__(self, bin_conf_name):
        self['name']=bin_conf_name

        conf=read_config(bin_conf_name)
        self.update(conf)

        self['nbin'] = len(self['bin_info'])

    def get_name(self):
        """
        the name for this bin scheme, just the identifier for this
        bin in the config name
        """
        return self['name']

    def get_nbin(self):
        """
        get the number of bins
        """
        return self['nbin']

    def get_label(self, binnum):
        """
        get the label for the bin
        """
        return self['bin_info'][binnum]['label']

    def get_bintags(self):
        """
        get the tags to average in each bin
        """
        bin_info=self['bin_info']
        
        bi0 = bin_info[0]['bins']
        bintags = [ri[0] for ri in bi0]

        return bintags

    def bin(self, data, jackreg_col=None):
        """
        bin the data and return the result
        """
        nrbin = data['rsum'][0].size
        if 'dsensum' in data.dtype.names:
            shear_style='lensfit'
        else:
            shear_style='reduced'

        bin_info=self['bin_info']
        bintags=self.get_bintags()

        print("bintags:",bintags)

        bs = averaging.lensbin_struct(nrbin,
                                      shear_style=shear_style,
                                      bintags=bintags,
                                      n=self['nbin'])

        print("len(data):",len(data))
        print("len(bin_info):",len(bin_info))
        for i,bi in enumerate(bin_info):
            range_info = bi['bins']

            comb,w = reduce_from_ranges_many(data,
                                             range_info,
                                             jackreg_col=jackreg_col,
                                             getind=True)
        
            #bs['nlenses'][i] = w.size
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            for bi in range_info:
                field_name=bi[0]
                mn='%s_mean' % field_name
                en='%s_err' % field_name
                sn='%s_sdev' % field_name
                rn='%s_range' % field_name

                # now the things we are averaging by lens weight
                mean,err,sdev = averaging.lens_wmom(data,
                                                    field_name,
                                                    ind=w,
                                                    sdev=True)
                bs[mn][i] = mean
                bs[en][i] = err
                bs[sn][i] = sdev
                bs[rn][i] = data[field_name][w].min(), data[field_name][w].max()

            i+=1

        return bs

    def bin_one(self, data, binnum, jackreg_col=None):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        range_info = self['bin_info'][binnum]['bins']
        comb = reduce_from_ranges_many(data,
                                       range_info,
                                       jackreg_col=jackreg_col,
                                       getind=True)
 
        return comb

    def select_bin(self, data, binnum):
        """
        Although not used by bin(), this is useful for other programs

        """

        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        range_info = self['bin_info'][binnum]['bins']

        logic = get_range_logic_many(data, range_info)

        w,=where(logic)
        return w



    def bin_label(self, binnum):
        return self.bin_info[binnum]['label']


def reduce_from_ranges_many(data, tags_and_ranges, jackreg_col=None, getind=False):
    """

    parameters
    ----------
    data: ndarray with fields
        The lensum data
    tags_and_ranges: list
        each element in the list is of the form
            (tagname, range, range_type)
    getind: bool
        If True, get the indices from data that are in range
    """

    logic = get_range_logic_many(data, tags_and_ranges)
    w,=where(logic)

    comb = averaging.average_lensums(data[w],jackreg_col=jackreg_col)

    if getind:
        return comb, w
    else:
        return comb
 
def get_range_logic_many(data, tags_and_ranges):
    """
    logic for multiple range sets
    """
    logic = numpy.ones(data.size, dtype='bool')
    for tag_range in tags_and_ranges:
        #print("    ",tag_range)
        tag, range, range_type = tag_range
        logic = logic & get_range_logic(data, tag, range, range_type)

    return logic

def get_range_logic(data, tag, brange, type):
    """
    logic for single range set
    """
    minval,maxval = brange
    if minval is None:
        minval = data[tag].min()
    if maxval is None:
        maxval = data[tag].max()

    #print(minval,maxval)
    if type == '[]':
        logic = (data[tag] >= minval) & (data[tag] <= maxval)
    elif type == '[)':
        logic = (data[tag] >= minval) & (data[tag] < maxval)
    elif type == '(]':
        logic = (data[tag] > minval) & (data[tag] <= maxval)
    elif type == '()':
        logic = (data[tag] > minval) & (data[tag] < maxval)
    else:
        raise ValueError("Bad range type: '%s'" % type)

    return logic

def define_bins(var_input, lastmin, alpha=0.6666, visible=False, prompt=False):
    """

    define bins assuming the data look similar to a schechter function.  The
    trick is dealing with the exponential cutoff: pick a "lastmin" such that
    defines the lower edge of the last bin.  For lower values we assume a power
    law, working downward keeping N*var^alpha = constant

    parameters
    ----------
    var: array
        sample data with wich to define bins
    lastmin: float
        lower edge of last bin
    alpha: float, optional
        Assumed power law at the lower end of distribution. Default 0.6666
    visible: bool, optional
        make a plot if True
    prompt: bool, optional
        prompt the user if set True
    """
    
    from biggles import FramedPlot, Curve, Point, PlotLabel
    var = var_input.copy()
    var.sort()

    # reverse sort, biggest to smallest
    var = numpy.fromiter(reversed(var), dtype='f8')

    ind = numpy.arange(var.size,dtype='i8')

    w_last, = where(var > lastmin)
    mvar_last = var[w_last].sum()/w_last.size

    print("wlast.size:",w_last.size,"mvar_last:",mvar_last)

    cref = 0.5*log10(w_last.size) + alpha*log10(mvar_last)

    # now look at mean var*N for rest by using cumulative
    # sums


    var = var[w_last[-1]:]
    var_last = lastmin
    binnum=0

    minbin=[]
    maxbin=[]
    while 1:
        nd = 1+numpy.arange(var.size)

        mvar = var.cumsum()/nd
        cval = 0.5*log10(nd) + alpha*log10(mvar)

        wthis, = where(cval < cref)

        var_min = var[wthis[-1]]
        var_max = var_last

        minbin.append(var_min)
        maxbin.append(var_max)
        if visible:
            print("numleft:",var.size,"wthis.size",wthis.size)

            plt=FramedPlot()
            curve = Curve(var, cval)
            plt.add(curve)

            oc = Curve([var.min(), var.max()],[cref,cref],color='blue')
            plt.add(oc)

            p = Point(var[wthis[-1]], cval[wthis[-1]], color='orange', type='filled circle')
            plt.add(p)

            binnum -= 1
            blab=PlotLabel(0.9,0.9,'bin %d' % binnum, halign='right')

            rlab=PlotLabel(0.9,0.85,r'%0.2f $ < var < $ %0.2f' % (var_min,var_max), halign='right')
            nlab=PlotLabel(0.9,0.80,'N: %d' % wthis.size, halign='right')

            plt.add(blab)
            plt.add(rlab)
            plt.add(nlab)
            plt.show()

            if prompt:
                key=raw_input('hit a key: ')
                if key == 'q':
                    return
        var_last = var_min
        var=var[wthis[-1]+1:]


        if len(var) == 0:
            break
    
    minbin = list(reversed(minbin))
    maxbin = list(reversed(maxbin))
    minbin.append(maxbin[-1])
    maxbin.append(None)
    
    minstr=[]
    maxstr=[]
    for i in xrange(len(minbin)):
        if minbin[i] is not None:
            minstr.append('%0.1f' % minbin[i])
        else:
            minstr.append('None')
        if maxbin[i] is not None:
            maxstr.append('%0.1f' % maxbin[i])
        else:
            maxstr.append('None')

    minstr = '[' + ', '.join(minstr) +']'
    maxstr = '[' + ', '.join(maxstr) +']'

    #for i in xrange(len(minbin)):
    #    print('%i %0.1f %0.1f' % (i+1,minbin[i],maxbin[i]))

    print("nbin:",len(minbin))
    print(minstr)
    print(maxstr)


