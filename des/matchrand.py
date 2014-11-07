"""
code to match the redshift distributions of randoms to data
"""
from __future__ import print_function
import numpy
from .files_common import *
from . import binning

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
        nbin=self.binner.get_nbin()

        for binnum in xrange(nbin):
            print("-"*70)
            print("matching bin: %d/%d" % (binnum+1,nbin))
            if self['match_type']=='z':
                self._match_one_by_z(binnum)
            else:
                self._match_one_by_selection(binnum)
        self._make_html()

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
        weights *= ( 1.0/weights.max() )

        self._print_effnum(weights)
        self._plot_zhist(rz, z[w], weights, binnum)

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

        self._print_effnum(weights)
        self._plot_zhist(rz, z[w], weights, binnum)

    def _print_effnum(self, weights):
        num=weights.size
        effnum = weights.sum()
        effperc = effnum/num
        print("effective number: %d/%d = %0.2f" % (effnum, num, effperc))

    def _plot_zhist(self, rz, z, weights, binnum):
        """
        plot the z hist for data, randoms, and weighted randoms
        """
        pngfile=get_match_weights_file(self['lens_run'],
                                       self['rand_run'],
                                       self['bin_scheme'],
                                       binnum,
                                       ext='png')

        tit=self.binner.get_label(binnum)
        tit+=' rand: '+self['rand_run']

        weighting.plot_results1d(rz, z, weights, self['binsize'], 
                                 pngfile=pngfile, title=tit,
                                 xlabel='z',
                                 label1='rand', label2='lenses',
                                 show=self['show'])

    def _make_html(self):
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

        weight_col=self['rrun_conf']['lens_conf'].get('extra_weight_col')
        if weight_col is not None:
            self.extra_weights=self.rdata[weight_col]
        else:
            self.extra_weights=None

    def _set_match_type(self):
        """
        set the match type and make sure it is valud
        """
        self['match_type'] = self['rrun_conf']['lens_conf']['match_type']
        if self['match_type'] not in ['z','selection']:
            raise ValueError("bad match type: '%s'" % self['match_type'])
