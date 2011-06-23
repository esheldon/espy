"""

Find objects with > N observations.  Start by limiting to the stripe82 runs,
and then just read through them. Save a collated file of 

    combined flux in each band
    single epoch fluxes.
    seeing

Based on the combined flux, we can do excellent s/g separation.  Then bin the
objects by their single-epoch flux and seeing and look at the distribution of
concentration.

order of things
---------------
collate_stripe82_runs


"""
from __future__ import print_function

import os
import numpy
from numpy import log10, unique, sqrt, exp
import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join
import sdsspy
from . import stomp_maps

import columns

def collate_stripe82_epochs(nrunmin=10, rlim=1.0):
    se=Stripe82Epochs(nrunmin, rlim)
    se.collate()

class Stripe82Epochs:
    def __init__(self, nrunmin=10, rlim=1.0):
        """
        nrunmin: scalar, optional
            nrunmin is a first limit on the number of runs we demand
            default 10
        rlim: scalar, optional
            limit on not extinction corrected r-band flux. Default is 1.0 (22.5)

        Run methods in this order
            collate()  collate the primary objects and epochs
            extract_psf_fwhm()  get psf_fwhm from the sweeps for the epochs
            collate_psf_fwhm()  collate the psf_fwhm into the epochs

        """
        self.nrunmin=nrunmin
        self.rlim=rlim

        # use the new sweeps with combined fluxes
        #os.environ['PHOTO_SWEEP'] = \
        #    path_join(os.environ['PHOTO_SWEEP_BASE'],'dr8_final_fcombine')

    def dir(self):
        outdir=os.environ['SGSEP_DIR']
        return outdir
    def coldir(self, type):
        d=self.dir()
        d = path_join(d, type+'.cols')
        return d
    def open_columns(self, type):
        d=self.coldir(type)
        return columns.Columns(d)

    def match_rg(self, rgrun):
        import lensing
        cols=self.open_columns('primary82')
        rgcols=lensing.regauss.open_columns(rgrun)
        print("reading photoid")
        photoid = cols['photoid'][:]
        print("  read",photoid.size)

        print("reading rg photoid")
        rg_photoid = rgcols['photoid'][:]
        print("  read",rg_photoid.size)

        print("doing match")
        #mrg,mc = eu.numpy_util.match(rg_photoid, photoid)
        mc,mrg = eu.numpy_util.match(photoid, rg_photoid)
        mc.tofile('/home/users/esheldon/tmp/mc.dat')
        mrg.tofile('/home/users/esheldon/tmp/mrg.dat')
        print("  matched",mrg.size)

        rgindex = numpy.zeros(photoid.size, dtype='i4') - 9999
        rgindex[mc] = numpy.array(mrg, dtype='i4')

        colname = 'rgmatch%s' % rgrun

        cols.write_column(colname, rgindex, create=True)

        """
        print("    cutting R > 1/3")
        R = rgcols['R_rg_r'][mrg]
        wkeep=where1(R > 0.3333)
        print("    kept",wkeep.size)
        mrg=mrg[wkeep]
        mc=mc[wkeep]

        return mc
        """

    def extract_psf_fwhm(self):
        """
        This could have been avoided had I added psf_fwhm to the epochs
        Not all will match the sweeps, we will have to collate the results
        """
        import es_sdsspy
        sgdir = self.dir()
        extract_file = path_join(sgdir, 'epochs82-psf-fwhm.rec')

        cols = self.open_columns('epochs82')

        extract_cols = ['run','rerun','camcol','field','id','psf_fwhm']

        pid = cols['photoid'][:]
        se=es_sdsspy.sweeps.SweepExtractor(pid, allow_nomatch=True)
        se.extract(extract_cols, extract_file)

    def collate_psf_fwhm(self):
        """
        Collate the results from matching epochs to sweeps to get psf_fwhm.
        """
        print("opening columns")
        cols = self.open_columns('epochs82')
        print("getting photoid")
        pid = cols['photoid'][:]

        sgdir = self.dir()
        extract_file = path_join(sgdir, 'epochs82-psf-fwhm.rec')
        print("reading psf_fwhm from:",extract_file)
        fwhm_data = eu.io.read(extract_file)



        print("matching")
        fwhm_pid = sdsspy.photoid(fwhm_data)
        m,ms = eu.numpy_util.match(pid, fwhm_pid)

        print("  matched %i/%i" % (m.size,pid.size))

        fwhm = numpy.zeros(pid.size, dtype='f4')

        print("writing psf_fwhm_g")
        fwhm[:] = -9999.
        fwhm[m] = fwhm_data['psf_fwhm'][ms,1]
        cols.write_column('psf_fwhm_g', fwhm, create=True)

        print("writing psf_fwhm_r")
        fwhm[:] = -9999.
        fwhm[m] = fwhm_data['psf_fwhm'][ms,2]
        cols.write_column('psf_fwhm_r', fwhm, create=True)

        print("writing psf_fwhm_i")
        fwhm[:] = -9999.
        fwhm[m] = fwhm_data['psf_fwhm'][ms,3]
        cols.write_column('psf_fwhm_i', fwhm, create=True)



    def collate(self):
        """

        Keep objects if have nrunmin matches in both model and psf flux, in
        *any* one of g,r,i

        """
        runs = stripe82_run_list()
        #runs=runs[where1(runs == 259)]
        
        pcoldir = self.coldir('primary82')
        ecoldir = self.coldir('epochs82')
        if os.path.exists(pcoldir):
            raise ValueError("primary coldir exists: '%s'\n    Start Fresh" % pcoldir)
        if os.path.exists(ecoldir):
            raise ValueError("epochs coldir exists: '%s'\n    Start Fresh" % ecoldir)

        pcols = columns.Columns(pcoldir)
        pcols.create()
        ecols = columns.Columns(ecoldir)
        ecols.create()

        # count in the epochs columns
        epochs_index = 0
        for run in runs:
            for camcol in [1,2,3,4,5,6]:
                fname=sdsspy.filename('photoCombineCamcol', run, camcol)
                if os.path.exists(fname):
                    data=sdsspy.read('photoCombineCamcol', run, camcol,
                                     verbose=True,lower=True)

                    # first select by number of runs matched, then do flag
                    # selection and such
                    print("  selecting r-band nrunmin >",self.nrunmin)

                    model_nuse = data['cmodel_nuse'][:,2]
                    psf_nuse   = data['psf_nuse'][:,2]
                    w=where1( (psf_nuse >= self.nrunmin) & (model_nuse >= self.nrunmin) )

                    print("    kept:",w.size)
                    if w.size > 0:
                        data = data[w]
                        print("  getting sweep matches for flags")
                        data = self.match_sweeps(data)

                        print("      found",data.size,"sweep matches")
                        if data.size > 0:
                            print("  selecting flags and rflux >",self.rlim)
                            w = where1( get_select_logic(data, self.rlim) )
                            print("    kept:",w.size)

                            if w.size > 0:
                                data = data[w]

                                print("    making primary output")
                                poutput = self.make_primary_output(data)

                                print("    gathering epochs")
                                epochs, epochs_indices, epochs_count = \
                                        self.gather_epochs(data, epochs_index)
                                poutput['epochs_index'] = epochs_indices
                                poutput['epochs_count'] = epochs_count

                                nexpected = self.calc_expected_nepoch(data)
                                if epochs.size != nexpected:
                                    #if nexpected-epochs.size == 1:
                                    #print("FOUND ERROR nexpected greater than epochs by 1. primary bug?")
                                    #else:
                                    #    raise ValueError("Expected",nexpected,"epochs, but found",epochs.size)
                                    print("ERROR: Expected",nexpected,"epochs, but found",epochs.size)

                                print("      found:",epochs.size)
                                print("    making epochs output")
                                eoutput = self.make_epochs_output(epochs)

                                pcols.write_columns(poutput)
                                ecols.write_columns(eoutput)

                                # keep track of index in epochs database
                                epochs_index += epochs.size



    def match_sweeps(self, data):
        gal=sdsspy.read('calibobj.gal',data['run'][0], data['camcol'][0],
                        lower=True)
        star=sdsspy.read('calibobj.star',data['run'][0], data['camcol'][0],
                         lower=True)
        photoid = sdsspy.photoid(data)
        gphotoid = sdsspy.photoid(gal)
        sphotoid = sdsspy.photoid(star)

        mdata_g, mg = eu.numpy_util.match(photoid, gphotoid)
        print("      gal matches:",mg.size)
        mdata_s, ms = eu.numpy_util.match(photoid, sphotoid)
        print("      star matches:",ms.size)

        nmatch = mdata_g.size + mdata_s.size
        if nmatch == 0:
            return numpy.array([])

        newdt = [('photoid','i8'),
                 ('objc_flags','i4'),
                 ('flags','i4',5),
                 ('flags2','i4',5),
                 ('calib_status','i4',5),
                 ('modelflux','f4',5),
                 ('modelflux_ivar','f4',5),
                 ('psfflux','f4',5),
                 ('psfflux_ivar','f4',5)]

        dtype = data.dtype.descr + newdt
        out = numpy.zeros(nmatch, dtype=dtype)

        #data_keep = numpy.zeros(nmatch, dtype=data.dtype)
        if mg.size > 0:
            gal = gal[mg]
            out['photoid'][0:mg.size] = gphotoid[mg]
            out['flags'][0:mg.size] = gal['flags']
            out['flags2'][0:mg.size] = gal['flags2']
            out['objc_flags'][0:mg.size] = gal['objc_flags']
            out['calib_status'][0:mg.size] = gal['calib_status']

            out['modelflux'][0:mg.size] = gal['modelflux']
            out['modelflux_ivar'][0:mg.size] = gal['modelflux_ivar']
            out['psfflux'][0:mg.size] = gal['psfflux']
            out['psfflux_ivar'][0:mg.size] = gal['psfflux_ivar']

            for n in data.dtype.names:
                out[n][0:mg.size] = data[n][mdata_g]

        if ms.size > 0:
            star=star[ms]
            out['photoid'][mg.size:] = sphotoid[ms]
            out['flags'][mg.size:] = star['flags']
            out['flags2'][mg.size:] = star['flags2']
            out['objc_flags'][mg.size:] = star['objc_flags']
            out['calib_status'][mg.size:] = star['calib_status']

            out['modelflux'][mg.size:] = star['modelflux']
            out['modelflux_ivar'][mg.size:] = star['modelflux_ivar']
            out['psfflux'][mg.size:] = star['psfflux']
            out['psfflux_ivar'][mg.size:] = star['psfflux_ivar']

            for n in data.dtype.names:
                out[n][mg.size:] = data[n][mdata_s]

        return out

    def calc_expected_nepoch(self,data):
        model_nuse = data['cmodel_nuse'][:,2]
        psf_nuse = data['psf_nuse'][:,2]
        nuse = numpy.where(model_nuse > psf_nuse, model_nuse, psf_nuse)
        return nuse.sum()

    def gather_epochs(self, data, epochs_index_current):
        photoid = sdsspy.photoid(data)

        fh, frev=eu.stat.histogram(data['field'], binsize=1, rev=True)

        epochlist = []
        epochs_indices = numpy.zeros(data.size)
        epochs_count = numpy.zeros(data.size)

        epochs_index = epochs_index_current
        for i in xrange(fh.size):
            if frev[i] != frev[i+1]:
                wfield = frev[ frev[i]:frev[i+1] ]
                field = data['field'][wfield[0]]

                epochs = sdsspy.read('photoEpochs', 
                                     data['run'][0], 
                                     data['camcol'][0], 
                                     field=field,
                                     verbose=False, 
                                     lower=True)

                # first extract those that were actually used according
                # to the criteria for selecting the primaries above
                w=where1(epochs['cmodel_used'][:,2] == 1)
                if w.size == 0:
                    raise ValueError("found none that were used!")

                epochs = epochs[w]
                    
                # extract those that match by photoid
                for j in xrange(wfield.size):
                    pid = photoid[wfield[j]]
                    w=where1(epochs['primary_photoid'] == pid)
                    if w.size == 0:
                        raise ValueError("no matches for pid",pid)

                    epochs_indices[wfield[j]] = epochs_index
                    epochs_count[wfield[j]] = w.size
                    epochs_index += w.size

                    epochlist.append(epochs[w])

        if len(epochlist) == 0:
            raise ValueError("Expected to find some epochs!")

        epochs = eu.numpy_util.combine_arrlist(epochlist)

        if (epochs_index-epochs_index_current) != epochs.size:
            raise ValueError("epoch index not incremented right amount:"
                             "%i vs %i" % (epochs_index-epochs_index_current,epochs.size))

        return epochs, epochs_indices, epochs_count
                
    def make_primary_output(self, data):
        dtype = self.prim_dtype()
        out = numpy.zeros(data.size, dtype=dtype)

        out['photoid'] = sdsspy.photoid(data)
        out['objc_type'] = data['objc_type']

        for ftype in ['psf','cmodel']:
            for filter in ['g','r','i']:
                fnum=sdsspy.FILTERNUM[filter]

                if ftype == 'cmodel':
                    # we couldn't get cmodel because stars don't have dev,exp
                    out['modelflux_'+filter] = data['modelflux'][:,fnum]
                    out['modelflux_ivar_'+filter] = data['modelflux_ivar'][:,fnum]
                else:
                    out[ftype+'flux_'+filter] = data[ftype+'flux'][:,fnum]
                    out[ftype+'flux_ivar_'+filter] = data[ftype+'flux_ivar'][:,fnum]

                out[ftype+'_nuse_'+filter] = data[ftype+'_nuse'][:,fnum]

                out[ftype+'flux_mean_'+filter] = data[ftype+'flux_mean'][:,fnum]
                out[ftype+'flux_mean_ivar_'+filter] = data[ftype+'flux_mean_ivar'][:,fnum]
        return out


    def make_epochs_output(self, data):
        dtype = self.epochs_dtype()
        out = numpy.zeros(data.size, dtype=dtype)

        out['primary_photoid'] = data['primary_photoid']

        out['photoid'] = sdsspy.photoid(data)
        out['objc_type'] = data['objc_type']

        for ftype in ['psf','cmodel']:
            for filter in ['g','r','i']:
                fnum=sdsspy.FILTERNUM[filter]
                out[ftype+'flux_'+filter] = data[ftype+'flux'][:,fnum]
                out[ftype+'flux_ivar_'+filter] = data[ftype+'flux_ivar'][:,fnum]

                out[ftype+'_used_'+filter] = data[ftype+'_used'][:,fnum]

        return out

    def epochs_dtype(self):
        dtype = [('photoid','i8'),
                 ('primary_photoid','i8'),

                 ('objc_type','i1'),

                 ('cmodel_used_g','i1'),
                 ('cmodelflux_g','f4'),
                 ('cmodelflux_ivar_g','f4'),

                 ('cmodel_used_r','i1'),
                 ('cmodelflux_r','f4'),
                 ('cmodelflux_ivar_r','f4'),

                 ('cmodel_used_i','i1'),
                 ('cmodelflux_i','f4'),
                 ('cmodelflux_ivar_i','f4'),

                 ('psf_used_g','i1'),
                 ('psfflux_g','f4'),
                 ('psfflux_ivar_g','f4'),

                 ('psf_used_r','i1'),
                 ('psfflux_r','f4'),
                 ('psfflux_ivar_r','f4'),

                 ('psf_used_i','i1'),
                 ('psfflux_i','f4'),
                 ('psfflux_ivar_i','f4')]
        return dtype

    def prim_dtype(self):
        dtype = [('photoid','i8'),
                 ('epochs_index','i4'), # index into the epochs
                 ('epochs_count','i4'), # number of objects in the epochs

                 ('objc_type','i1'),

                 # can't use cmodel because stars don't have exp/dev!
                 ('modelflux_g','f4'),
                 ('modelflux_ivar_g','f4'),
                 ('modelflux_r','f4'),
                 ('modelflux_ivar_r','f4'),
                 ('modelflux_i','f4'),
                 ('modelflux_ivar_i','f4'),

                 ('psfflux_g','f4'),
                 ('psfflux_ivar_g','f4'),
                 ('psfflux_r','f4'),
                 ('psfflux_ivar_r','f4'),
                 ('psfflux_i','f4'),
                 ('psfflux_ivar_i','f4'),

                 ('cmodel_nuse_g','i2'),
                 ('cmodel_nuse_r','i2'),
                 ('cmodel_nuse_i','i2'),

                 ('cmodelflux_mean_g','f4'),
                 ('cmodelflux_mean_ivar_g','f4'),
                 ('cmodelflux_mean_r','f4'),
                 ('cmodelflux_mean_ivar_r','f4'),
                 ('cmodelflux_mean_i','f4'),
                 ('cmodelflux_mean_ivar_i','f4'),

                 ('psf_nuse_g','i2'),
                 ('psf_nuse_r','i2'),
                 ('psf_nuse_i','i2'),

                 ('psfflux_mean_g','f4'),
                 ('psfflux_mean_ivar_g','f4'),
                 ('psfflux_mean_r','f4'),
                 ('psfflux_mean_ivar_r','f4'),
                 ('psfflux_mean_i','f4'),
                 ('psfflux_mean_ivar_i','f4')]
        return dtype

    def load_primary_avg_gri(self, indices=None):
        cols = self.open_columns('primary82')

        cmodel_mean_g = cols['cmodelflux_mean_g'][:]
        cmodel_mean_ivar_g = cols['cmodelflux_mean_ivar_g'][:]

        psf_mean_g = cols['psfflux_mean_g'][:]
        psf_mean_ivar_g = cols['psfflux_mean_ivar_g'][:]

        cmodel_mean_r = cols['cmodelflux_mean_r'][:]
        cmodel_mean_ivar_r = cols['cmodelflux_mean_ivar_r'][:]

        psf_mean_r = cols['psfflux_mean_r'][:]
        psf_mean_ivar_r = cols['psfflux_mean_ivar_r'][:]

        cmodel_mean_i = cols['cmodelflux_mean_i'][:]
        cmodel_mean_ivar_i = cols['cmodelflux_mean_ivar_i'][:]

        psf_mean_i = cols['psfflux_mean_i'][:]
        psf_mean_ivar_i = cols['psfflux_mean_ivar_i'][:]

        cmodelflux, cmodelflux_ivar = avg_gri(cmodel_mean_g, cmodel_mean_ivar_g,
                                            cmodel_mean_r, cmodel_mean_ivar_r,
                                            cmodel_mean_i, cmodel_mean_ivar_i)
        psfflux, psfflux_ivar = avg_gri(psf_mean_g, psf_mean_ivar_g,
                                        psf_mean_r, psf_mean_ivar_r,
                                        psf_mean_i, psf_mean_ivar_i)

        cbad=where1(numpy.isnan(cmodelflux))
        if cbad.size != 0:
            raise ValueError("nan found in cmodelflux")
        pbad=where1(numpy.isnan(psfflux))
        if pbad.size != 0:
            raise ValueError("nan found in psfflux")



        if indices is not None:
            cmodelflux=cmodelflux[indices]
            cmodelflux_ivar=cmodelflux_ivar[indices]
            psfflux=psfflux[indices]
            psfflux_ivar=psfflux_ivar[indices]
        return cmodelflux, cmodelflux_ivar, psfflux, psfflux_ivar


    def plot_nuse(self, cmodel_nuse, psf_nuse):
        import biggles

        hcmodel=eu.stat.histogram(cmodel_nuse,more=True)
        hpsf=eu.stat.histogram(psf_nuse,more=True)
        plt=biggles.FramedPlot()

        p_hcmodel = biggles.Histogram(hcmodel['hist'], x0=hcmodel['low'][0], binsize=1, color='blue')
        p_hpsf = biggles.Histogram(hpsf['hist'], x0=hpsf['low'][0], binsize=1, color='red')

        p_hcmodel.label = 'cmodel'
        p_hpsf.label = 'psf'

        key=biggles.PlotKey(0.9,0.9,[p_hcmodel,p_hpsf], halign='right')

        plt.add(p_hcmodel, p_hpsf, key)
        plt.xlabel = 'number of exposures'
        plt.show()

    def primary_gaussmix(self, prompt=True, show=True, ngauss=4, 
                         nrunmin=10, combine_gri=True, nperbin=50000, rg=None):
        """
        Use em algorithm to get a gaussian mixture model.  Default 4 gaussians

        flux to mags of 1.6 is about 22
        flux log10(flux) mags
        1    0.0         22.5
        1.6  0.2         22.0
        3    0.48        21.3
        4    0.6         21.0
        """
        import biggles
        import scipy.optimize
        import lensing

        from scikits.learn import mixture

        cols = self.open_columns('primary82')

        cmodel_nuse = cols['cmodel_nuse_r'][:]
        psf_nuse = cols['psf_nuse_r'][:]

        # make plot *before* cut
        #self.plot_nuse(cmodel_nuse, psf_nuse)

        w=where1((cmodel_nuse >= nrunmin) & (psf_nuse >= nrunmin))
        if w.size == 0:
            raise ValueError("none passed nrunmim of",nrunmin)

        if rg is not None:
            print("matching rg")
            rgmcol = 'rgmatch%s' % rg
            rgm = cols[rgmcol][w]
            wrg = where1(rgm >= 0)
            w=w[wrg]
            rgm=rgm[wrg]

            print("    with matches:",w.size)

            rgcols=lensing.regauss.open_columns(rg)
            print("    cutting R > 1/3")
            R = rgcols['R_rg_r'][rgm]
            wkeep=where1(R > 0.3333)
            print("      kept",wkeep.size)

            w=w[wkeep]


        cmodel_nuse=cmodel_nuse[w]
        psf_nuse=psf_nuse[w]


        if combine_gri:
            print("getting combined gri for mean flux")
            fname='mean flux gri'
            cmodel_mean_r, cmodel_mean_ivar_r, psf_mean_r, psf_mean_ivar_r = \
                self.load_primary_avg_gri(indices=w)
        else:
            fname='mean flux r'
            print("getting r-band for mean flux")
            cmodel_mean_r = cols['cmodelflux_mean_r'][w]
            psf_mean_r = cols['psfflux_mean_r'][w]

        model_r = cols['modelflux_r'][w]
        psf_r = cols['psfflux_r'][w]


        logcmodel_mean_r = log10( cmodel_mean_r.clip(0.001,cmodel_mean_r.max()) )

        # 1.9 is about 21.8
        minc=-.2
        maxc=1.
        cbinsize=0.005
        cmean = 1.0-psf_mean_r/cmodel_mean_r
        c = 1.0-psf_r/model_r
        xtitle='1-psf/cmodel'

        w=where1((cmean > minc) & (cmean < maxc) )

        hdict = eu.stat.histogram(logcmodel_mean_r[w], 
                                  min=0.2, 
                                  nperbin=nperbin, 
                                  more=True)
        h=hdict['hist']
        rev=hdict['rev']

        nruntext = biggles.PlotLabel(0.9,0.9,'nexp >= %d' % nrunmin, halign='right')
        nbin = hdict['hist'].size
        output=None
        for i in xrange(hdict['hist'].size):
            if rev[i] != rev[i+1]:


                plt=biggles.FramedPlot()

                wbin=rev[ rev[i]:rev[i+1] ]
                wbin = w[wbin]


                min_logr = hdict['low'][i]
                max_logr = hdict['high'][i]
                max_magr = 22.5-2.5*min_logr
                min_magr = 22.5-2.5*max_logr

                cmean_h=eu.stat.histogram(cmean[wbin], min=minc, max=maxc, binsize=cbinsize, more=True)

                p_cmean_h=biggles.Histogram(cmean_h['hist'], x0=minc, binsize=cbinsize)
                p_cmean_h.label = fname

                c_h=eu.stat.histogram(c[wbin], min=minc, max=maxc, binsize=cbinsize, more=True)
                p_c_h=biggles.Histogram(c_h['hist'], x0=minc, binsize=cbinsize, color='grey')
                p_c_h.label = 'SE flux'




                key = biggles.PlotKey(0.9,0.5,[p_cmean_h,p_c_h], halign='right')
                plt.add(p_c_h, p_cmean_h)

                plt.add(key)
                

                plt.xtitle=xtitle

                ftext=biggles.PlotLabel(0.9,0.85,r'$%0.2f < log_{10}(f) < %0.2f$' % (min_logr,max_logr),
                                        halign='right')
                mtext=biggles.PlotLabel(0.9,0.8,r'$%0.2f > mag > %0.2f$' % (22.5-2.5*min_logr,22.5-2.5*max_logr),
                                        halign='right')
                plt.add(ftext)
                plt.add(mtext)
                plt.add(nruntext)


                gmm = mixture.GMM(n_states=ngauss)
                data = cmean[wbin]
                data.shape = (data.size,1)
                gmm.fit(data, n_iter=10, min_covar=1.e-6)

                mg = MultiGauss(gmm)


                yfit = mg.eval(cmean_h['center'])
                renorm = wbin.size/float(yfit.sum())
                yfit *= renorm
                plt.add(biggles.Curve(cmean_h['center'], yfit))

                if ngauss <= 5:
                    colors = ['red','magenta','blue','purple','orange']
                else:
                    import pcolors
                    colors = pcolors.rainbow(ngauss,'hex')
                for i in xrange(ngauss):
                    yfit = mg.evalone(cmean_h['center'],i)
                    yfit *= renorm
                    plt.add(biggles.Curve(cmean_h['center'], yfit, color=colors[i]))

                if show:
                    plt.show()

                pfile=gaussmix_plotfile(ngauss,min_logr, max_logr, ext='png')
                print("writing plot file:",pfile)
                plt.write_img(800,800,pfile)

                if prompt and show is not False:
                    k=raw_input('hit a key (q to quit): ')
                    if k == 'q':
                        return

                weights = gmm.weights
                means = gmm.means[:,0]
                sigmas = [sqrt(gmm.covars[i][0][0]) for i in xrange(ngauss)]

                if output is None:
                    output = numpy.zeros(nbin, dtype=[('min_logr','f8'),
                                                      ('max_logr','f8'),
                                                      ('min_magr','f8'),
                                                      ('max_magr','f8'),
                                                      ('clow','f8',cmean_h['hist'].size),
                                                      ('chigh','f8',cmean_h['hist'].size),
                                                      ('chist_mean','i8',cmean_h['hist'].size),
                                                      ('chist_se','i8',cmean_h['hist'].size),
                                                      ('weights','f8',ngauss),
                                                      ('means','f8',ngauss),
                                                      ('sigmas','f8',ngauss)])


                output['min_logr'][i] = min_logr
                output['max_logr'][i] = max_logr
                output['min_magr'][i] = min_magr
                output['max_magr'][i] = max_magr

                output['clow'][i] = c_h['low']
                output['chigh'][i] = c_h['high']
                output['chist_se'][i] = c_h['hist']
                output['chist_mean'][i] = cmean_h['hist']

                output['weights'][i] = weights
                output['means'][i] = means
                output['sigmas'][i] = sigmas

        outfile=gaussmix_file(ngauss)
        eu.io.write(outfile, output, verbose=True)



    def primary_twogauss(self, prompt=True, show=True, 
                         nrunmin=10, combine_gri=True, nperbin=50000, linear=True, rg=None):
        """
        fit a double gaussian to the stacked fluxes for primary objects

        flux to mags of 1.6 is about 22
        flux log10(flux) mags
        1    0.0         22.5
        1.6  0.2         22.0
        3    0.48        21.3
        4    0.6         21.0
        """
        import biggles
        import scipy.optimize
        import lensing

        cols = self.open_columns('primary82')

        cmodel_nuse = cols['cmodel_nuse_r'][:]
        psf_nuse = cols['psf_nuse_r'][:]

        # make plot *before* cut
        #self.plot_nuse(cmodel_nuse, psf_nuse)

        w=where1((cmodel_nuse >= nrunmin) & (psf_nuse >= nrunmin))
        if w.size == 0:
            raise ValueError("none passed nrunmim of",nrunmin)

        if rg is not None:
            print("matching rg")
            rgmcol = 'rgmatch%s' % rg
            rgm = cols[rgmcol][w]
            wrg = where1(rgm >= 0)
            w=w[wrg]
            rgm=rgm[wrg]

            print("    with matches:",w.size)

            rgcols=lensing.regauss.open_columns(rg)
            print("    cutting R > 1/3")
            R = rgcols['R_rg_r'][rgm]
            wkeep=where1(R > 0.3333)
            print("      kept",wkeep.size)

            w=w[wkeep]


        cmodel_nuse=cmodel_nuse[w]
        psf_nuse=psf_nuse[w]


        objc_type = cols['objc_type'][w]
        if combine_gri:
            print("getting combined gri for mean flux")
            fname='mean flux gri'
            cmodel_mean_r, cmodel_mean_ivar_r, psf_mean_r, psf_mean_ivar_r = \
                self.load_primary_avg_gri(indices=w)
        else:
            fname='mean flux r'
            print("getting r-band for mean flux")
            cmodel_mean_r = cols['cmodelflux_mean_r'][w]
            psf_mean_r = cols['psfflux_mean_r'][w]

        model_r = cols['modelflux_r'][w]
        psf_r = cols['psfflux_r'][w]


        logcmodel_mean_r = log10( cmodel_mean_r.clip(0.001,cmodel_mean_r.max()) )

        # 1.9 is about 21.8
        if linear:
            minc=-.1
            maxc=1.
            cbinsize=0.01
            cmean = 1.0-psf_mean_r/cmodel_mean_r
            c = 1.0-psf_r/model_r
            xtitle='1-psf/cmodel'
        else:
            minc=-0.05
            maxc=0.4
            cbinsize=0.005
            cmean = -log10(psf_mean_r/cmodel_mean_r)
            c     = -log10(psf_r/model_r)

            #cmean -= 0.03
            #c -= 0.03
            off=0.1
            power=0.4
            cmean = (cmean+off)**power - off**power
            c = (c+off)**power - off**power
            cbinsize = 0.002
            minc = -0.1
            maxc=0.6
            xtitle=r'$log_{10}(cmodel/psf)$'

        w=where1((cmean > minc) & (cmean < maxc) )

        hdict = eu.stat.histogram(logcmodel_mean_r[w], 
                                  min=0.2, 
                                  nperbin=nperbin, 
                                  more=True)
        h=hdict['hist']
        rev=hdict['rev']

        nruntext = biggles.PlotLabel(0.9,0.9,'nexp >= %d' % nrunmin, halign='right')
        nbin = hdict['hist'].size
        output=None
        for i in xrange(hdict['hist'].size):
            if rev[i] != rev[i+1]:


                plt=biggles.FramedPlot()

                wbin=rev[ rev[i]:rev[i+1] ]
                wbin = w[wbin]


                min_logr = hdict['low'][i]
                max_logr = hdict['high'][i]
                max_magr = 22.5-2.5*min_logr
                min_magr = 22.5-2.5*max_logr

                cmean_h=eu.stat.histogram(cmean[wbin], min=minc, max=maxc, binsize=cbinsize, more=True)

                p_cmean_h=biggles.Histogram(cmean_h['hist'], x0=minc, binsize=cbinsize)
                p_cmean_h.label = fname

                c_h=eu.stat.histogram(c[wbin], min=minc, max=maxc, binsize=cbinsize, more=True)
                p_c_h=biggles.Histogram(c_h['hist'], x0=minc, binsize=cbinsize, color='grey')
                p_c_h.label = 'SE flux'


                if output is None:
                    output = numpy.zeros(nbin, dtype=[('min_logr','f8'),
                                                      ('max_logr','f8'),
                                                      ('min_magr','f8'),
                                                      ('max_magr','f8'),
                                                      ('clow','f8',cmean_h['hist'].size),
                                                      ('chigh','f8',cmean_h['hist'].size),
                                                      ('chist_mean','i8',cmean_h['hist'].size),
                                                      ('chist_se','i8',cmean_h['hist'].size),
                                                      ('pars','f8',6)])



                key = biggles.PlotKey(0.9,0.5,[p_cmean_h,p_c_h], halign='right')
                plt.add(p_c_h, p_cmean_h)

                wstar=where1(objc_type[wbin] == 6); 
                if wstar.size > 0:
                    wstar=wbin[wstar]
                    cstar_h=eu.stat.histogram(cmean[wstar], min=minc, max=maxc, binsize=cbinsize, more=True)
                    p_cstar_h=biggles.Histogram(cstar_h['hist'], x0=minc, binsize=cbinsize, color='blue')
                    plt.add(p_cstar_h)


                wgal=where1(objc_type[wbin] == 3); 
                if wgal.size > 0:
                    wgal=wbin[wgal]
                    cgal_h=eu.stat.histogram(cmean[wgal], min=minc, max=maxc, binsize=cbinsize, more=True)
                    p_cgal_h=biggles.Histogram(cgal_h['hist'], x0=minc, binsize=cbinsize, color='red')
                    plt.add(p_cgal_h)


                plt.add(key)
                

                plt.xtitle=xtitle

                ftext=biggles.PlotLabel(0.9,0.85,r'$%0.2f < log_{10}(f) < %0.2f$' % (min_logr,max_logr),
                                        halign='right')
                mtext=biggles.PlotLabel(0.9,0.8,r'$%0.2f > mag > %0.2f$' % (22.5-2.5*min_logr,22.5-2.5*max_logr),
                                        halign='right')
                plt.add(ftext)
                plt.add(mtext)
                plt.add(nruntext)



                yerr=numpy.sqrt(cmean_h['hist'])
                yerr=yerr.clip(1.0, yerr.max())
                gf = TwoGaussClass(cmean_h['center'], cmean_h['hist'], yerr)
                parguess=(0.5*wbin.size*cbinsize, 0.025, 0.01, 
                          0.5*wbin.size*cbinsize, 0.125, 0.05)
                #res = scipy.optimize.fmin(gf.chi2, parguess)
                res = scipy.optimize.fmin_powell(gf.chi2, parguess)
                print(res)

                yfit = gf.twogauss(res)
                yfit1 = gf.gauss(res[0:3])
                yfit2 = gf.gauss(res[3:])

                plt.add(biggles.Curve(cmean_h['center'], yfit))
                if res[1] < res[4]:
                    plt.add(biggles.Curve(cmean_h['center'], yfit1, color='blue'))
                    plt.add(biggles.Curve(cmean_h['center'], yfit2, color='red'))
                else:
                    plt.add(biggles.Curve(cmean_h['center'], yfit2, color='blue'))
                    plt.add(biggles.Curve(cmean_h['center'], yfit1, color='red'))
 
               
                

                if show:
                    plt.show()

                pfile=twogauss_plotfile(min_logr, max_logr, ext='png')
                print("writing plot file:",pfile)
                plt.write_img(800,800,pfile)

                if prompt:
                    k=raw_input('hit a key (q to quit): ')
                    if k == 'q':
                        return

                output['min_logr'][i] = min_logr
                output['max_logr'][i] = max_logr
                output['min_magr'][i] = min_magr
                output['max_magr'][i] = max_magr

                output['clow'][i] = c_h['low']
                output['chigh'][i] = c_h['high']
                output['chist_se'][i] = c_h['hist']
                output['chist_mean'][i] = cmean_h['hist']

                output['pars'][i] = res

        outfile='~/oh/star-galaxy-separation/twogauss-fits/twogauss-fits.rec'
        eu.io.write(outfile, output, verbose=True)


def gaussmix_dir():
    dir=os.path.expanduser('~/oh/star-galaxy-separation/gaussmix')
    return dir
def gaussmix_file(ngauss):
    dir=gaussmix_dir()
    if not os.path.exists(dir):
        os.makedirs(dir)
    f=os.path.join(dir, 'gaussmix-%d.rec' % (ngauss,))
    return f
def gaussmix_plotfile(ngauss, fmin, fmax, ext='png'):
    dir=gaussmix_dir()
    dir=os.path.join(dir, "plots")
    if not os.path.exists(dir):
        os.makedirs(dir)
    f=os.path.join(dir, 'gaussmix-%d-%0.2f-%0.2f.%s' % (ngauss,fmin,fmax,ext))
    return f



def twogauss_dir():
    dir=os.path.expanduser('~/oh/star-galaxy-separation/twogauss-fits')
    return dir
def twogauss_plotfile(fmin, fmax, ext='png'):
    dir=twogauss_dir()
    f=os.path.join(dir, "plots", 'twogauss-%0.2f-%0.2f.%s' % (fmin,fmax,ext))
    return f

class MultiGauss:
    def __init__(self, gmm):
        """
        Takes a gmm object

        eval(x) returns normalized evaluation of gaussians
        """
        self.gmm = gmm

    def eval(self, x):
        model = numpy.zeros(x.size, dtype='f8')

        gmm = self.gmm
        for i in xrange(gmm.n_states):
            mean = gmm.means[i,0]
            var = gmm.covars[i][0][0]
            weight = gmm.weights[i]
            g = self.gauss(x, mean, var)
            model[:] += weight*g


        return model

    def evalone(self, x, i):
        """
        Just evaluate one of the gaussians
        """

        gmm = self.gmm
        mean = gmm.means[i,0]
        var = gmm.covars[i][0][0]
        weight = gmm.weights[i]
        return weight*self.gauss(x, mean, var)


    def gauss(self, x, mean, var):
        siginv2 = 1./var

        g = exp(-0.5*(x-mean)**2*siginv2)
        norm = sqrt(siginv2/2./numpy.pi)

        g *= norm
        return g



class TwoGaussClass:
    def __init__(self, x, y, yerr=None):
        self.x = numpy.array(x, dtype='f8', copy=False)
        self.y = numpy.array(y, dtype='f8', copy=False)

        if yerr is None:
            self.yerr = numpy.ones(x.size, dtype='f8')
        else:
            self.yerr=numpy.array(yerr, dtype='f8', copy=False)
        self.ivar = 1.0/self.yerr**2

    def chi2(self, pars):
        if pars[0] <= 1.e-6 or pars[2] <= 1.e-6 or pars[3] <= 1.e-6 or pars[5] <= 1.e-6:
            return 1.e15

        yfunc = self.twogauss(pars)
        chi2 = self.ivar*(self.y-yfunc)**2
        return chi2.sum()


    def twogauss(self, pars):
        g1 = self.gauss(pars[0:3] )
        g2 = self.gauss(pars[3:] )

        return g1+g2

    def gauss(self, pars):
        siginv2 = 1./pars[2]**2

        g = exp(-0.5*(self.x-pars[1])**2*siginv2)

        norm = pars[0]*sqrt(siginv2/2./numpy.pi)
        #norm = A

        g *= norm
        return g


def stripe82_run_list():
    """
    We need to get these fcombined runs
    """
    win = sdsspy.window.Window()
    fl = win.read('flist')

    map=stomp_maps.load('boss','basic')
    inmap = map.Contains(fl['ra'],fl['dec'],'eq')

    w=where1(  (fl['dec'] >= -1.25) 
             & (fl['dec'] <= 1.25) 
             & ((fl['ra'] > 300) | (fl['ra'] < 60)) 
             & (fl['score'] > 0.1) 
             & (inmap==1))
    uruns = unique(fl['run'][w])
    print("number of unique runs:",uruns.size)
    return uruns


def avg_gri(flux_g, ivar_g, 
            flux_r, ivar_r, 
            flux_i, ivar_i):
    """

    NOTE: you should use the same weights (ivar) for the psf and modelfluxes
    for things like the concentration.

    if you don't want to use an object in a particular band,  you should
    set the ivar to zero *before* calling this function

    fluxes are clipped between [0.001,1.e4], which is magnitude [12.5,30]
    """

    # clip the fluxes on the high and low end
    # this is mag between 12.5 and 30
    flux_g = flux_g.clip(0.001, 1.e4)
    flux_r = flux_r.clip(0.001, 1.e4)
    flux_i = flux_i.clip(0.001, 1.e4)

    # clip ivar as well, although this should not really be a problem as ivar
    # seems to always be well behaved

    ivar_g = ivar_g.clip(0.0, 70)
    ivar_r = ivar_r.clip(0.0, 70)
    ivar_i = ivar_i.clip(0.0, 70)

    ivarsum = ivar_g + ivar_r + ivar_i

    fsum = flux_g*ivar_g + flux_r*ivar_r + flux_i*ivar_i

    flux = fsum/ivarsum

    return flux, ivarsum


def calc_c(cmodelflux, psfflux, log=False):
    if log:
        return -log10(psfflux/cmodelflux)
    else:
        return 1.0-psfflux/cmodelflux


def get_select_logic(objs, rflux_lim):
    """
    flux_lim:
        1 is 22.5
        3 is 21.3
        4 is 21.0
    """
    from . import select

    sel = select.Selector(objs)
    fl = sel.flag_logic()


    #flux = sdsspy.make_cmodelflux(objs, doivar=False)
    rflux = objs['modelflux'][:,2]
    flux_logic = rflux > rflux_lim

    return flux_logic & fl


def read_test_data():
    import esutil as eu
    from . import select
    gal=eu.io.read('~/data/boss/calibObj-000756-3-gal.fits',lower=True)
    star=eu.io.read('~/data/boss/calibObj-000756-3-star.fits',lower=True)

    gflags = get_select_logic(gal,4.0)
    sflags = get_select_logic(star,4.0)

    gw = where1(gflags)
    sw = where1(sflags)

    ntot = gw.size + sw.size
    dt=[('origtype','i4'),
        ('modelflux','f4',5),
        ('modelflux_ivar','f4',5),
        ('psfflux','f4',5),
        ('psfflux_ivar','f4',5)]

    data = numpy.zeros(ntot, dtype=dt)
    data['origtype'][0:gw.size] = 3
    data['origtype'][gw.size:] = 6
    for n in ['modelflux','modelflux_ivar','psfflux','psfflux_ivar']:
        data[n][0:gw.size] = gal[n][gw]
        data[n][gw.size:] = star[n][sw]

    return data

def test(data=None, logc=False):
    from biggles import Histogram, FramedPlot, PlotKey, Table
    if data is None:
        data = read_test_data()

    
    modelflux, modelflux_ivar = avg_gri(data['modelflux'][:,1],data['modelflux_ivar'][:,1],
                                        data['modelflux'][:,2],data['modelflux_ivar'][:,2],
                                        data['modelflux'][:,3],data['modelflux_ivar'][:,3])
    #psfflux, psfflux_ivar = avg_gri(data['psfflux'][:,1],data['psfflux_ivar'][:,1],
    #                                data['psfflux'][:,2],data['psfflux_ivar'][:,2],
    #                                data['psfflux'][:,3],data['psfflux_ivar'][:,3])
    psfflux, psfflux_ivar = avg_gri(data['psfflux'][:,1],data['modelflux_ivar'][:,1],
                                    data['psfflux'][:,2],data['modelflux_ivar'][:,2],
                                    data['psfflux'][:,3],data['modelflux_ivar'][:,3])

    fmin_log=1.e-3
    fmax_log=4.


    tab = Table(2,2)
    binsize=0.05
    col=0
    for type in ['modelflux','psfflux']:

        flux_plt = FramedPlot()

        h = eu.stat.histogram(log10(data[type][:,1]), min=fmin_log, max=fmax_log, binsize=binsize)
        gmod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='green')
        gmod_h.label = 'g '+type

        h = eu.stat.histogram(log10(data[type][:,2]), min=fmin_log, max=fmax_log, binsize=binsize)
        rmod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='red')
        rmod_h.label = 'r '+type

        h = eu.stat.histogram(log10(data[type][:,3]), min=fmin_log, max=fmax_log, binsize=binsize)
        imod_h = Histogram(h, x0=fmin_log, binsize=binsize, color='magenta')
        imod_h.label = 'i '+type


        if type == 'modelflux':
            h = eu.stat.histogram(log10(modelflux), min=fmin_log, max=fmax_log, binsize=binsize)
        else:
            h = eu.stat.histogram(log10(psfflux), min=fmin_log, max=fmax_log, binsize=binsize)

        mod_h = Histogram(h, x0=fmin_log, binsize=binsize, width=2)
        mod_h.label = 'combined '+type

        key = PlotKey(0.5,0.9,[gmod_h, rmod_h, imod_h, mod_h])
        
        flux_plt.add(gmod_h, rmod_h, imod_h, mod_h, key)
        flux_plt.xlabel = 'flux'

        tab[0,col] = flux_plt

        col += 1



    col=0
    for logc in [False,True]:
        if logc:
            xmin=-0.1
            #xmax=1
            xmax=0.6
            binsize = 0.01
        else:
            xmin=-0.1
            xmax=1.0
            binsize = 0.01

        gc = calc_c(data['modelflux'][:,1], data['psfflux'][:,1], log=logc)
        rc = calc_c(data['modelflux'][:,2], data['psfflux'][:,2], log=logc)
        ic = calc_c(data['modelflux'][:,3], data['psfflux'][:,3], log=logc)
        allc = calc_c(modelflux, psfflux, log=logc)
        
        c_plt = FramedPlot()

        h = eu.stat.histogram(gc, min=xmin, max=xmax, binsize=binsize)
        gch = Histogram(h, x0=xmin, binsize=binsize, color='green')
        gch.label = 'g'

        h = eu.stat.histogram(rc, min=xmin, max=xmax, binsize=binsize)
        rch = Histogram(h, x0=xmin, binsize=binsize, color='red')
        rch.label = 'r'

        h = eu.stat.histogram(ic, min=xmin, max=xmax, binsize=binsize)
        ich = Histogram(h, x0=xmin, binsize=binsize, color='magenta')
        ich.label = 'i'


        h = eu.stat.histogram(allc, min=xmin, max=xmax, binsize=binsize)
        ch = Histogram(h, x0=xmin, binsize=binsize, width=2)
        ch.label = 'combined '

        key = PlotKey(0.7,0.9,[gch, rch, ich, ch])
        
        c_plt.add(gch, rch, ich, ch, key)
        if logc:
            c_plt.xlabel = r'$-log_{10}(psf/model)$'
        else:
            c_plt.xlabel = '1-psf/model'

        tab[1,col] = c_plt
        col += 1

    tab.show()



def test_fit2gauss():
    import biggles
    import mcmc
    import scipy.optimize

    N1 = 10000
    N2 = 10000

    m1 = 0.0
    sig1 = 0.05
    m2 = 0.7
    sig2 = 0.6

    data1 = numpy.random.standard_normal(N1)*sig1 + m1
    data2 = numpy.random.standard_normal(N2)*sig2 + m2
    data = numpy.zeros(N1+N2, dtype='f8')
    data[0:N1] = data1
    data[N1:] = data2

    binsize=sig1*0.2
    hdict = eu.stat.histogram(data, binsize=binsize, more=True)
    hdict1 = eu.stat.histogram(data1, binsize=binsize, more=True)
    hdict2 = eu.stat.histogram(data2, binsize=binsize, more=True)

    hp = biggles.Histogram(hdict['hist'], x0=hdict['low'][0], binsize=binsize)
    hp1 = biggles.Histogram(hdict1['hist'], x0=hdict1['low'][0], binsize=binsize,color='blue')
    hp2 = biggles.Histogram(hdict2['hist'], x0=hdict2['low'][0], binsize=binsize,color='red')

    data_hplt=biggles.FramedPlot()
    data_hplt.add(hp, hp1, hp2)

    x = numpy.array(hdict['center'], dtype='f8')
    y = numpy.array(hdict['hist'], dtype='f8')
    yerr = sqrt(hdict['hist'])
    #yerr = yerr.clip(1., yerr.max())
    yerr=numpy.ones(x.size,dtype='f8')

    # normalization needs the delta x in it!
    dx = x[1]-x[0]
    parguess = [N1*dx, m1, sig1, N2*dx, m2, sig2]

    gf = TwoGaussClass(x, y, yerr)
    parguess=tuple(parguess)

    res = scipy.optimize.fmin(gf.chi2, parguess)

    best_N1 = res[0]
    best_m1 = res[1]
    best_sig1 = res[2]
    best_N2 = res[3]
    best_m2 = res[4]
    best_sig2 = res[5]

    print("parguess:",parguess)
    print("best N1:",best_N1)
    print("best m1:",best_m1)
    print("best sig1:",best_sig1)
    print("best N2:",best_N2)
    print("best m2:",best_m2)
    print("best sig2:",best_sig2)

    yfit = gf.twogauss(res)
    if res[1] > res[4]:
        yfit1 = gf.gauss(res[3:])
        yfit2 = gf.gauss(res[0:3])
    else:
        yfit1 = gf.gauss(res[0:3])
        yfit2 = gf.gauss(res[3:])
    data_hplt.add(biggles.Curve(x, yfit))
    data_hplt.add(biggles.Curve(x, yfit1,color='blue'))
    data_hplt.add(biggles.Curve(x, yfit2, color='red'))
    data_hplt.show()


    return res




