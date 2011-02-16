import os
from sys import stdout,stderr
import pprint

import numpy
from numpy import where

import esutil as eu
from esutil.ostools import path_join, expand_path

import columns
import sdssgal

import zphot

# select the zphot input data from a columns database
class ColumnSelector:
    """

    zcs = ColumnSelector(sample, nchunk=None)

    sample is any string, e.g. '01'

    zcs.select()
    zcs.write()

    """
    def __init__(self, photo_sample=None, nchunk=None):
        """
        E.g.  procrun='prim04'
        """

        self.prefix='zinput'
        self.photo_sample = photo_sample
        self.nchunk = nchunk

        self.data = None
        self.dir = None
        self.keep_indices=None

        self.conf = None

        if self.photo_sample is not None:
            self.init(photo_sample)

    def init(self, photo_sample):
        self.photo_sample = photo_sample

        pzdir = zphot.photoz_dir()
        self.dir = path_join(pzdir,'inputs',self.photo_sample)
        stdout.write("photo dir: %s\n" % self.dir)

        self.conf = zphot.read_config('zinput',photo_sample)
        if self.conf['photo_sample'] != photo_sample:
            raise ValueError("photo_sample in config is not '%s'" \
                             % photo_sample)
        pprint.pprint(self.conf)

        if self.conf['weighting']:
            self.conf['filetype'] = 'dat'

    def load_cols(self, procrun):
        self.cols = sdssgal.open_columns(procrun)

    def write(self, chunk=None):
        if self.keep_indices is None:
            raise ValueError("run select() first")

        if chunk is None and self.nchunk is not None:
            # loop and write the chunks
            for chunk in xrange(self.nchunk):
                self.write(chunk)
            return

        if not os.path.exists(self.dir):
            os.makedirs(self.dir)

        outfile=self.filename(chunk=chunk)
        stdout.write("\nWill write to file: %s\n" % outfile)

        self.load_input_columns(chunk=chunk)

        if self.conf['filetype'] == 'fits':
            stdout.write("Ensuring to big endian in place\n")
            eu.numpy_util.to_big_endian(self.data, inplace=True)
        else:
            stdout.write("Ensuring native in place\n")
            eu.numpy_util.to_native(self.data, inplace=True)

        stdout.write("Writing to file: %s\n" % outfile)
        stdout.flush()


        if self.conf['filetype'] == 'dat':
            r = eu.recfile.Open(outfile, "w", delim=' ')
            r.write(self.data)
            r.close()
        else:
            eu.io.write(outfile, self.data, clobber=True)
 
    def get_chunk_indices(self, chunk):
        ntot = len(self.keep_indices)
        nperchunk = ntot/self.nchunk
        beg = chunk*nperchunk
        if chunk == (self.nchunk-1):
            end = ntot
        else:
            end = beg + nperchunk
        ind = self.keep_indices[beg:end]
        return ind


    def select(self):
        if self.conf is None:
            raise ValueError("run init with a photo_sample name")

        # for shorthand notation
        self.load_cols(self.conf['procrun'])

        c = self.cols

        stdout.write("\nReading flags\n")
        flags = c['flags'][:]

        stdout.write('Getting binned logic\n')
        binned_logic = sdssgal.binned_logic_ri(flags)
        w,=where(binned_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del flags

        logic = binned_logic
        
        stdout.write("Reading objc_flags\n")
        objc_flags = c['objc_flags'][:]
        stdout.write('Getting objc_flags1 logic\n')
        objc_flags1_logic = sdssgal.object1_logic(objc_flags)
        w,=where(objc_flags1_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del objc_flags

        logic = logic & objc_flags1_logic

        # always do model mag logic
        mmin = self.conf['modelmag_min']
        mmax = self.conf['modelmag_max']

        stdout.write("Cutting to reasonable model mags: [%s,%s]\n" \
                     % (mmin,mmax))
        stdout.write("Reading modelmag_dered\n")
        mag = c['modelmag_dered'][:]
        stdout.write('Getting mag logic\n')
        mag_logic = \
                  (mag[:,0] > mmin) & (mag[:,0] < mmax) \
                & (mag[:,1] > mmin) & (mag[:,1] < mmax) \
                & (mag[:,2] > mmin) & (mag[:,2] < mmax) \
                & (mag[:,3] > mmin) & (mag[:,3] < mmax) \
                & (mag[:,4] > mmin) & (mag[:,4] < mmax)
        w,=where(mag_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del mag

        logic = logic & mag_logic


        rmin = self.conf['cmodel_rmin']
        rmax = self.conf['cmodel_rmax']

        stdout.write("Cutting to cmodel_r: [%s,%s]\n" % (rmin,rmax))
        stdout.write("Reading cmodelmag_dered\n")
        cmag = c['cmodelmag_dered'][:]
        stdout.write('Getting cmag logic\n')
        cmag_logic = (cmag[:,2] > rmin) & (cmag[:,2] < rmax)
        w,=where(cmag_logic)
        stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
        del cmag

        logic = logic & cmag_logic


        masktype = self.conf['masktype']
        if masktype is not None:
            mask_colname = 'in'+masktype
            stdout.write("Reading %s (and good)\n" % mask_colname)
            inmask = c[mask_colname][:]
            stdout.write("Getting %s mask logic\n" % masktype)
            mask_logic = (inmask == 1)
            w,=where(mask_logic)
            stdout.write("    %s/%s passed\n" % (w.size,c['flags'].size))
            del inmask

            logic = logic & mask_logic
        
        w,=where(logic)
        stdout.write("A total of %s/%s passed\n" % (w.size,c['flags'].size))

        self.keep_indices = w


    def filename(self, chunk=None):

        outdir = self.dir

        format='%s-%s'
        fpars=[self.prefix,
               self.photo_sample]

        if chunk is not None:
            format += '-chunk%03i'
            fpars.append(chunk)

        format += '.%s'
        fpars.append(self.conf['filetype'])

        outfile = format % tuple(fpars)

        outfile = path_join(outdir,outfile) 
        return outfile

    def load_input_columns(self, chunk=None):

        if self.keep_indices is None:
            raise ValueError("run select() first")

        if self.nchunk is not None:
            if chunk is None:
                raise ValueError("please send chunk")
            indices = self.get_chunk_indices(chunk)
        else:
            indices = self.keep_indices

        del self.data
        if self.conf['weighting']:
            colnames = ['photoid',
                        'cmodelmag_dered',
                        'modelmag_dered']

            stdout.write("Reading columns: %s\n" % colnames)
            tmp = self.cols.read_columns(colnames, 
                                         rows=indices, 
                                         verbose=True)

            # we want the r-band cmodelmag and the model colors
            dt = zphot.weighting.photo_dtype()
            data = numpy.zeros(tmp.size, dtype=dt)

            data['photoid'] = tmp['photoid']
            data['cmodelmag_dered_r'] = tmp['cmodelmag_dered'][:,2]
            data['model_umg'] = \
                tmp['modelmag_dered'][:,0] - tmp['modelmag_dered'][:,1]
            data['model_gmr'] = \
                tmp['modelmag_dered'][:,1] - tmp['modelmag_dered'][:,2]
            data['model_rmi'] = \
                tmp['modelmag_dered'][:,2] - tmp['modelmag_dered'][:,3]
            data['model_imz'] = \
                tmp['modelmag_dered'][:,3] - tmp['modelmag_dered'][:,4]
            del tmp
            self.data = data
        else:

            colnames = ['photoid',
                        'modelmag_dered',
                        'modelmag_dered_err']

            stdout.write("Reading columns: %s\n" % colnames)
            self.data = self.cols.read_columns(colnames, 
                                               rows=indices,
                                               verbose=True)



