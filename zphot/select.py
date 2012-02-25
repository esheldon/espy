from __future__ import print_function
import os
from sys import stdout,stderr
import pprint

import numpy
from numpy import where

import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join, expand_path

import columns
import es_sdsspy

import zphot

# select the zphot input data from a columns database


class ColumnSelector:
    """

    zcs = ColumnSelector(sample, nchunk=None)

    sample is any string, e.g. '01'

    zcs.select()
    zcs.write()

    """
    def __init__(self, photo_sample, nchunk=None):
        """
        E.g.  procrun='prim04'
        """

        self.prefix='zinput'
        self.photo_sample = photo_sample
        self.nchunk = nchunk

        self.data = None
        self.keep_indices=None

        pzdir = zphot.photoz_dir()
        self.dir = path_join(pzdir,'inputs',self.photo_sample)
        print("photo dir: ",self.dir)

        self.conf = zphot.read_config('zinput',photo_sample)
        if self.conf['photo_sample'] != photo_sample:
            raise ValueError("photo_sample in config is not '%s'" \
                             % photo_sample)
        pprint.pprint(self.conf)

        if self.conf['weighting']:
            self.conf['filetype'] = 'dat'

    def load_cols(self, primary=True):
        if primary:
            self.cols = es_sdsspy.sweeps_collate.open_columns('primgal')
        else:
            self.cols = es_sdsspy.sweeps_collate.open_columns('gal')

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
        print("\nWill write to file:", outfile)

        self.load_input_columns(chunk=chunk)

        if self.conf['filetype'] == 'fits':
            print("Ensuring to big endian in place")
            eu.numpy_util.to_big_endian(self.data, inplace=True)
        else:
            print("Ensuring native in place")
            eu.numpy_util.to_native(self.data, inplace=True)

        print("Writing to file:",outfile)
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


    def select(self, primary=True, for_training=False):
        """

        DONT USE for_training
        If for_training=True, a looser cut on the r magnitude is applied so we
        don't have a sharp edge in the histogram at the same spot as the photometric
        sample. might not help, we'll see.

        """
        if self.conf is None:
            raise ValueError("run init with a photo_sample name")

        # for shorthand notation
        self.load_cols(primary=primary)
        c = self.cols

        #
        # first cut everything down to mask and r < 21.8
        #

        masktype = self.conf['masktype']
        mask_colname = 'in'+masktype
        print("Reading %s (and good)" % mask_colname)
        inmask = c[mask_colname][:]
        print("Getting %s mask logic" % masktype)
        mask_logic = (inmask == 1)
        wmask,=where(mask_logic)
        ntot=inmask.size
        print("    %s/%s passed (%0.2f)" % (wmask.size,ntot,wmask.size/float(ntot)))


        rmin = self.conf['cmodel_rmin']
        if for_training:
            self.conf['cmodel_rmax'] = 22.1
        rmax = self.conf['cmodel_rmax']

        print("Cutting to cmodel_r: [%s,%s]" % (rmin,rmax))
        print("Reading cmodelmag_dered")
        cmag = c['cmodelmag_dered_r'][:]
        print('Getting cmag logic:',rmax)
        cmag_logic = (cmag > rmin) & (cmag < rmax)
        wrmag,=where(cmag_logic)
        print("    %s/%s passed (%0.2f)" % (wrmag.size,ntot,wrmag.size/float(ntot)))

        wmagmask, = where(mask_logic & cmag_logic)
        print("    %s/%s passed both (%0.2f)" % (wmagmask.size,ntot,wmagmask.size/float(ntot)))

        del inmask
        del cmag

        #
        # Now watch the rest of the cuts to see how it breaks down
        #

        print("**\n** NOTE From here on, we have cut to the mask and rmag\n**")

        print("\nReading flags,objc_flags,calib_status")
        data = c.read_columns(['objc_flags',
                               'flags_r','flags_i',
                               'calib_status_u',
                               'calib_status_g',
                               'calib_status_r',
                               'calib_status_i',
                               'calib_status_z'], rows=wmagmask)
        ntot = data.size
                               
        selector = es_sdsspy.select.Selector(data)

        print("Getting binned logic")
        binned_logic = selector.binned_logic()
        w,=where(binned_logic)
        print("    %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))

        print("Getting calib logic")
        calib_logic = selector.calib_logic()
        w,=where(calib_logic)
        print("    %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))

        print("Getting object1 logic")
        object1_logic = selector.object1_logic()
        w,=where(object1_logic)
        print("    %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))


        del data


        # always do model mag logic
        mmin = self.conf['modelmag_min']
        mmax = self.conf['modelmag_max']

        print("Cutting to reasonable model mags: [%s,%s]" \
                     % (mmin,mmax))
        print("Reading modelmag_dered")
        mag_u = c['modelmag_dered_u'][wmagmask]
        mag_g = c['modelmag_dered_g'][wmagmask]
        mag_r = c['modelmag_dered_r'][wmagmask]
        mag_i = c['modelmag_dered_i'][wmagmask]
        mag_z = c['modelmag_dered_z'][wmagmask]
        print('Getting mag logic')
        mag_logic = \
                  (mag_u > mmin) & (mag_u < mmax) \
                & (mag_g > mmin) & (mag_g < mmax) \
                & (mag_r > mmin) & (mag_r < mmax) \
                & (mag_i > mmin) & (mag_i < mmax) \
                & (mag_z > mmin) & (mag_z < mmax)
        w,=where(mag_logic)
        print("    %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))
        print("        Breakdown:")

        w=where1((mag_u > mmin) & (mag_u < mmax))
        print("            u %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))
        w=where1((mag_g > mmin) & (mag_g < mmax))
        print("            g %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))
        w=where1((mag_r > mmin) & (mag_r < mmax))
        print("            r %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))
        w=where1((mag_i > mmin) & (mag_i < mmax))
        print("            i %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))
        w=where1((mag_z > mmin) & (mag_z < mmax))
        print("            z %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))




        print("Cutting thing_id >= 0")
        thing_id = c['thing_id'][wmagmask]
        thing_id_logic = (thing_id >= 0)
        w,=where(thing_id_logic)
        print("    %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))


        logic = \
            binned_logic \
            & calib_logic \
            & object1_logic \
            & mag_logic \
            & thing_id_logic


        
        w,=where(logic)
        print("A total of %s/%s passed (%0.2f)" % (w.size,ntot,w.size/float(ntot)))

        self.keep_indices = wmagmask[w]


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
                        'cmodelmag_dered_r',
                        'modelmag_dered_u',
                        'modelmag_dered_g',
                        'modelmag_dered_r',
                        'modelmag_dered_i',
                        'modelmag_dered_z']

            extra_columns = self.conf.get('extra_columns',[])
            for cdict in extra_columns:
                colnames.append(cdict['name'])



            print("Reading columns:")
            pprint.pprint(colnames)
            tmp = self.cols.read_columns(colnames, 
                                         rows=indices, 
                                         verbose=True)

            # we want the r-band cmodelmag and the model colors
            dt = zphot.weighting.photo_dtype(self.photo_sample)
            data = numpy.zeros(tmp.size, dtype=dt)

            data['photoid'] = tmp['photoid']
            data['cmodelmag_dered_r'] = tmp['cmodelmag_dered_r']
            data['model_umg'] = \
                tmp['modelmag_dered_u'] - tmp['modelmag_dered_g']
            data['model_gmr'] = \
                tmp['modelmag_dered_g'] - tmp['modelmag_dered_r']
            data['model_rmi'] = \
                tmp['modelmag_dered_r'] - tmp['modelmag_dered_i']
            data['model_imz'] = \
                tmp['modelmag_dered_i'] - tmp['modelmag_dered_z']

            for cdict in extra_columns:
                name=cdict['name']
                data[name] = tmp[name]
            del tmp

            self.data = data
        else:

            raise ValueError("fix the non-weighting version")
            colnames = ['photoid',
                        'modelmag_dered',
                        'modelmag_dered_err']

            print("Reading columns: ", colnames)
            self.data = self.cols.read_columns(colnames, 
                                               rows=indices,
                                               verbose=True)


class ColumnSelectorOld:
    """

    zcs = ColumnSelector(sample, nchunk=None)

    sample is any string, e.g. '01'

    zcs.select()
    zcs.write()

    """
    def __init__(self, photo_sample, nchunk=None):
        """
        E.g.  procrun='prim04'
        """

        self.prefix='zinput'
        self.photo_sample = photo_sample
        self.nchunk = nchunk

        self.data = None
        self.keep_indices=None

        pzdir = zphot.photoz_dir()
        self.dir = path_join(pzdir,'inputs',self.photo_sample)
        print("photo dir: ",self.dir)

        self.conf = zphot.read_config('zinput',photo_sample)
        if self.conf['photo_sample'] != photo_sample:
            raise ValueError("photo_sample in config is not '%s'" \
                             % photo_sample)
        pprint.pprint(self.conf)

        if self.conf['weighting']:
            self.conf['filetype'] = 'dat'

    def load_cols(self, procrun):
        import sdssgal
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
        import sdssgal
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
            dt = zphot.weighting.photo_dtype(self.photo_sample)
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



