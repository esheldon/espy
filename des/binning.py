from .files_common import *
from . import averaging

def bin_run(run, bin_conf_name):
    """
    Do the binning and write out a file
    """
    import fitsio

    data = read_collated(run)

    binner=binner(bin_conf_name)
    res = binner.bin(data)

    fname=get_binned_file(run, bin_conf_name)
    d=get_binned_dir(run, bin_conf_name)
    if not os.path.exists(d):
        print("Making dir:",d)
        os.makedirs(d)
    
    print("writing:",fname)
    fitsio.write(fname, res, clobber=True)

class Binner(dict):
    """
    bin any
    """

    def __init__(self, bin_conf_name, **keys):
        self.update(keys)
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


    def bin(self, data):
        """
        bin the data and return the result
        """
        nrbin = data['rsum'][0].size

        bin_info=self['bin_info']

        bi0 = bin_info[0]['bins']
        bintags = [ri[0] for ri in bi0]
        print("bintags:",bintags)
        bs = averaging.lensbin_struct(nrbin,
                                      bintags=bintags,
                                      n=self['nbin'])

        print("len(data):",len(data))
        print("len(bin_info):",len(bin_info))
        for i,bi in enumerate(bin_info):
            range_info = bi['bins']

            comb,w = reduce_from_ranges_many(data,
                                             range_info,
                                             getind=True)
        
            bs['nlenses'][i] = w.size
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

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        range_info = self['bin_info'][binnum]['bins']
        comb = reduce_from_ranges_many(data,
                                       range_info,
                                       getind=True)
 
        return comb



    def bin_label(self, binnum):
        return self.bin_info[binnum]['label']


def reduce_from_ranges_many(data, tags_and_ranges, getind=False):
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

    logic = numpy.ones(data.size, dtype='bool')
    for tag_range in tags_and_ranges:
        print("    ",tag_range)
        tag, range, range_type = tag_range
        logic = logic & get_range_logic(data, tag, range, range_type)

    w=where1(logic)

    comb = averaging.average_lensums(data[w])

    if getind:
        return comb, w
    else:
        return comb
 
def get_range_logic(data, tag, brange, type):
    minval,maxval = brange
    if minval is None:
        minval = data[tag].min()
    if maxval is None:
        maxval = data[tag].max()

    print(minval,maxval)
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


