"""
code to match the redshift distributions of randoms to data
"""
from __future__ import print_function
import numpy
from .files import *
from . import binning
from . import averaging

import weighting

class RandomMatcher(dict):
    """
    match randoms by z or selection

    example
    -------
    rm = RandomMatcher(lens_run, rand_run, bin_scheme)
    rm.match()

    outputs
    -------
    weights files and binned files for the randoms are written
    """
    def __init__(self, lens_run, rand_run, bin_scheme, binsize=0.0, show=False):
        self['lens_run']=lens_run
        self['rand_run']=rand_run
        self['bin_scheme']=bin_scheme

        self['binsize']=binsize

        self['lrun_conf'] = cascade_config(lens_run)
        self['rrun_conf'] = cascade_config(rand_run)

        self._set_match_type()

        self['bin_conf'] = read_config(bin_scheme)

        self['show']=show

        self.binner=binning.Binner(self['bin_scheme'])
        
        self._make_dirs()
        self._set_data()

    def match(self):
        """
        match according to the requested match type
        """
        import fitsio

        nbin=self.binner.get_nbin()
        bs=self._get_binned_struct()

        fname=get_match_weights_file(self['lens_run'],self['rand_run'],self['bin_scheme'])
        print("opening:",fname)
        with fitsio.FITS(fname,'rw',clobber=True) as fits:
            for binnum in xrange(nbin):
                print("-"*70)
                print("matching bin: %d/%d" % (binnum+1,nbin))

                w,weights,comb=self._run_matcher(binnum)
                fits.write(weights)

                self._print_effnum(weights)
                self._plot_zhist(self.rz, self.z[w], weights, binnum)

                for n in comb.dtype.names:
                    bs[n][binnum] = comb[n][0]

        self._write_binned(bs)
        self._write_html()

    def _run_matcher(self, binnum):
        """
        run the appropriate matcher
        """
        if self['match_type']=='z':
            w,weights=self._match_one_by_z(binnum)
        elif self['match_type']=='selection':
            w,weights=self._match_one_by_selection(binnum)
        else:
            w,weights=self._match_one_by_selection_and_z(binnum)

        print("    combining")
        jackreg_col=self['rrun_conf']['lens_conf'].get('jackreg_col',None)
        print("        using jackknifes:",jackreg_col)
        comb = averaging.average_lensums(self.rdata,
                                         weights=weights,
                                         jackreg_col=jackreg_col)
        return w,weights,comb

    def _match_one_by_z(self, binnum):
        """
        match redshifts in this bin
        """
        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)

        weights = weighting.hist_match(rz, z[w], self['binsize'],
                                        extra_weights1=self.extra_weights)

        return w,weights

    def _match_one_by_selection(self, binnum):
        """
        match by applying the same selection to real and random data

        This is opposed to the match_by_z method which forces the redshift
        distributions to match
        """


        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)
        rw=self.binner.select_bin(self.rdata, binnum)

        weights=numpy.zeros(nrand)
        weights[rw]=1.0

        if self.extra_weights is not None:
            weights[:] *= self.extra_weights[:]

        return w,weights

    def _match_one_by_selection_and_z(self, binnum):
        """
        match by applying the same selection to real and random data, as
        well as redshift matching
        """


        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)
        rw=self.binner.select_bin(self.rdata, binnum)

        weights0=numpy.zeros(nrand)
        weights0[rw]=1.0

        if self.extra_weights is not None:
            weights0[:] *= self.extra_weights[:]


        weights = weighting.hist_match(rz, z[w], self['binsize'],
                                       extra_weights1=weights0)
        return w,weights


    def _print_effnum(self, weights):
        num=weights.size
        effnum = weights.sum()
        effperc = effnum/num
        print("    effective number: %d/%d = %0.2f" % (effnum, num, effperc))

    def _plot_zhist(self, rz, z, weights, binnum):
        """
        plot the z hist for data, randoms, and weighted randoms
        """
        pngfile=get_match_weights_file(self['lens_run'],
                                       self['rand_run'],
                                       self['bin_scheme'],
                                       binnum=binnum,
                                       ext='png')

        tit=self.binner.get_label(binnum)
        tit+=' rand: '+self['rand_run']

        print("    writing:",pngfile)

        if weights is None:
            weights=numpy.ones(rz.size)

        weighting.plot_results1d(rz, z, weights, self['binsize'], 
                                 pngfile=pngfile, title=tit,
                                 xlabel='z',
                                 label1='rand', label2='lenses',
                                 show=self['show'])

    def _write_html(self):
        """
        make an html file with all the z hist plots
        """
        from glob import glob
        d=get_match_weights_dir(self['lens_run'],
                                self['rand_run'],
                                self['bin_scheme'])
        command='cd %s; im2html -p *.png > zhist.html' % d

        print("making html file")
        os.system(command)

    def _write_binned(self, bs):
        """
        write the binned lensum struct
        """
        import fitsio
        fname=get_match_binned_file(self['lens_run'],
                                    self['rand_run'],
                                    self['bin_scheme'])

        d=os.path.dirname(fname)
        try:
            os.makedirs(d)
        except:
            pass

        print("writing binned:",fname)
        fitsio.write(fname,bs,clobber=True)



    def _make_dirs(self):
        """
        make all necessary dirs
        """
        # this will make the binned dir as well
        d=get_match_weights_dir(self['lens_run'],
                                self['rand_run'],
                                self['bin_scheme'])
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)

    def _set_data(self):
        """
        read the catalogs
        """
        self.data=read_collated(self['lens_run'])
        self.rdata=read_collated(self['rand_run'])

        self.z=self.data[self['lrun_conf']['lens_conf']['z_col']]
        self.rz=self.rdata[self['rrun_conf']['lens_conf']['z_col']]

        weight_col=self['rrun_conf']['lens_conf'].get('extra_weight_col',None)
        if weight_col is not None:
            self.extra_weights=self.rdata[weight_col]
        else:
            self.extra_weights=None

    def _set_match_type(self):
        """
        set the match type and make sure it is valud
        """
        self['match_type'] = self['rrun_conf']['lens_conf']['match_type']
        if self['match_type'] not in ['z','selection','selection-and-z']:
            raise ValueError("bad match type: '%s'" % self['match_type'])

    def _get_binned_struct(self):
        """
        get the random summed and averaged struct
        """

        nrad = self['lrun_conf']['lens_conf']['nbin']
        nbin=self.binner.get_nbin()
        bin_info=self['bin_conf']['bin_info']
        bintags=self.binner.get_bintags()

        print("bintags:",bintags)

        shear_style=self['lrun_conf']['source_conf']['shear_style']
        bs = averaging.lensbin_struct(nrad,
                                      shear_style=shear_style,
                                      bintags=bintags,
                                      n=nbin)

        return bs

class RandomMatcherRemove(RandomMatcher):
    """
    match randoms by z, using remove method

    example
    -------
    rm = RandomMatcherRemove(lens_run, rand_run, bin_scheme)
    rm.match()

    outputs
    -------
    index files and binned files for the randoms are written
    """

    def match(self):
        """
        match according to the requested match type
        """
        import fitsio

        nbin=self.binner.get_nbin()
        bs=self._get_binned_struct()

        fname=get_match_weights_file(self['lens_run'],self['rand_run'],self['bin_scheme'])
        print("opening:",fname)
        with fitsio.FITS(fname,'rw',clobber=True) as fits:
            for binnum in xrange(nbin):
                print("-"*70)
                print("matching bin: %d/%d" % (binnum+1,nbin))

                w,wr,weights,comb=self._run_matcher(binnum)
                if weights is not None:
                    print("    writing kept ones with weights")
                    st=numpy.zeros(wr.size, dtype=[('index','i8'),('weight','f8')])
                    st['index']=wr
                    st['weight']=weights
                    fits.write(st)
                else:
                    print("    writing kept ones ")
                    fits.write(wr)

                self._print_perc(wr)
                self._plot_zhist(self.rz[wr], self.z[w], weights, binnum)

                for n in comb.dtype.names:
                    bs[n][binnum] = comb[n][0]

        self._write_binned(bs)
        self._write_html()

    def _run_matcher(self, binnum):
        """
        run the appropriate matcher
        """
        if self['match_type']=='z':
            w,wr,weights=self._match_one_by_z(binnum)
        else:
            w,wr,weights=self._match_one_by_selection_and_z(binnum)

        print("    combining")
        print("        weights","None" if weights is None else "not None")
        comb = averaging.average_lensums(self.rdata[wr], weights=weights)
        return w,wr,weights,comb

    def _match_one_by_z(self, binnum):
        """
        match redshifts in this bin
        """
        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)

        wr,weights = self._do_z_match(rz, z[w], weights=self.extra_weights)
        return w,wr,weights

    def _match_one_by_selection_and_z(self, binnum):
        """
        match by applying the same selection to real and random data, as
        well as redshift matching
        """

        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)
        wr0=self.binner.select_bin(self.rdata, binnum)

        wr,weights=self._do_z_match(rz[wr0], z[w], weights=self.extra_weights)

        wr=wr0[wr]
        return w,wr,weights

    def _do_z_match(self, rz, z, weights=None):
        """
        all the other match methods call this one
        """

        wr = weighting.hist_match_remove(rz, z, self['binsize'],
                                         extra_weights1=weights)

        if weights is None:
            weights_out=None
        else:
            weights_out=weights[wr]

        return wr,weights_out


    def _print_perc(self, wr):
        nkeep=wr.size
        ntot=self.rz.size
        perc=nkeep/float(ntot)
        print("    percent kept: %d/%d = %0.2f" % (nkeep,ntot,perc))


    def _write_html(self):
        """
        make an html file with all the z hist plots
        """
        from glob import glob
        d=get_match_weights_dir(self['lens_run'],
                                self['rand_run'],
                                self['bin_scheme'])
        command='cd %s; im2html -p *.png > zhist.html' % d

        print("making html file")
        os.system(command)

    def _set_match_type(self):
        """
        set the match type and make sure it is valud
        """
        self['match_type'] = self['rrun_conf']['lens_conf']['match_type']
        if self['match_type'] not in ['z','selection-and-z']:
            raise ValueError("bad match type for remove: '%s'" % self['match_type'])
