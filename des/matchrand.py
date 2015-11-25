"""
code to match the redshift distributions of randoms to data
"""
from __future__ import print_function
import numpy
from .files import *
from . import binning
from . import averaging


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
        if self['match_type']=='selection':
            self.match_select()
        else:
            self.match_weighted()

    def match_select(self):
        """
        match according to the requested match type
        """
        import fitsio

        nbin=self.binner.get_nbin()
        bs=self._get_binned_struct()

        for binnum in xrange(nbin):
            print("-"*70)
            print("matching bin: %d/%d" % (binnum+1,nbin))

            w,wr,rand_weights,comb=self._run_matcher(binnum)

            # rand_weights should be None
            assert rand_weights is None,"rand_weights should be non for select method"
            self._plot_zhist(self.rz[wr], self.z[w], rand_weights, binnum)
            self._plot_radec(self.data[w], self.rdata[wr], binnum)

            for n in comb.dtype.names:
                bs[n][binnum] = comb[n][0]

        self._write_binned(bs)
        self._write_html()


    def match_weighted(self):
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

                w,wr,rand_weights,comb=self._run_matcher(binnum)
                fits.write(rand_weights)

                self._print_effnum(rand_weights)
                self._plot_zhist(self.rz[wr], self.z[w], rand_weights, binnum)
                self._plot_radec(self.data[w], self.rdata[wr], binnum)

                for n in comb.dtype.names:
                    bs[n][binnum] = comb[n][0]

        self._write_binned(bs)
        self._write_html()

    def _run_matcher(self, binnum):
        """
        run the appropriate matcher
        """
        col=self.get('extra_weight_col',None)
        if col is not Noen and self['match_type'] != 'selection':
            raise NotImplementedError("implement extra_weight_col for match "
                                      "type '%s'" % self['match_type'])

        if self['match_type']=='z':
            # just force the redshift distributions to match
            # weights is matched to w
            w,rand_weights=self._match_one_by_z(binnum)
            wr=numpy.arange(self.rdata.size)

        elif self['match_type']=='selection':
            print("    using selection only")
            w,wr=self._match_one_by_selection(binnum)

            if col is not None:
                print("    adding weights: '%s'" % col)
                rand_weights=self.data[col][wr]
            else:
                rand_weights=None

        elif self['match_type'] == 'selection-and-z':
            # apply selection *and* force redshift distributions to match
            # weights is matched to w
            w,wr,rand_weights=self._match_one_by_selection_and_z(binnum)

        else:
            raise ValueError("bad match_type: '%s'" % self['match_type'])

        print("    combining")

        jackreg_col=self['rrun_conf']['lens_conf'].get('jackreg_col',None)
        print("        using jackknifes:",jackreg_col)

        comb = averaging.average_lensums(self.rdata[wr],
                                         weights=rand_weights, # matches w or is None
                                         jackreg_col=jackreg_col)
        return w,wr,rand_weights,comb

    def _match_one_by_z(self, binnum):
        """
        match redshifts in this bin, but apply no selection to
        the randoms
        """
        import weighting

        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)

        rand_weights = weighting.hist_match(rz, z[w], self['binsize'])

        return w,rand_weights

    def _match_one_by_selection(self, binnum):
        """
        match by applying the same selection to real and random data

        no weights are derived to force a match on the redshift distributions
        """

        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)
        wr=self.binner.select_bin(self.rdata, binnum)

        return w,wr

    def _match_one_by_selection_and_z(self, binnum):
        """
        match by applying the same selection to real and random data, as
        well as redshift matching
        """

        import weighting

        z=self.z
        rz=self.rz
        nrand=rz.size

        w=self.binner.select_bin(self.data, binnum)
        wr=self.binner.select_bin(self.rdata, binnum)

        rand_weights = weighting.hist_match(rz[wr], z[w], self['binsize'])
        return w,wr,rand_weights


    def _print_effnum(self, weights):
        num=weights.size
        effnum = weights.sum()
        effperc = effnum/num
        print("    effective number: %d/%d = %0.2f" % (effnum, num, effperc))

    def _plot_zhist(self, rz, z, rand_weights, binnum):
        """
        plot the z hist for data, randoms, and weighted randoms
        """
        import weighting

        pngfile=get_match_weights_file(self['lens_run'],
                                       self['rand_run'],
                                       self['bin_scheme'],
                                       binnum=binnum,
                                       ext='png')

        tit=self.binner.get_label(binnum)
        tit+=' rand: '+self['rand_run']

        print("    writing:",pngfile)

        if rand_weights is None:
            rand_weights=numpy.ones(rz.size)

        weighting.plot_results1d(rz, z, rand_weights, self['binsize'], 
                                 pngfile=pngfile, title=tit,
                                 xlabel='z',
                                 label1='rand', label2='lenses',
                                 show=self['show'])

    def _plot_radec(self, data, rdata, binnum):
        """
        plot the z hist for data, randoms, and weighted randoms
        """
        import esutil as eu
        import biggles
        import converter

        pngfile=get_match_weights_file(self['lens_run'],
                                       self['rand_run'],
                                       self['bin_scheme'],
                                       binnum=binnum,
                                       ext='png')

        pngfile=pngfile.replace('weights.png','radec.png')
        epsfile=pngfile.replace('.png','.eps')

        title=self.binner.get_label(binnum)

        print("    writing:",pngfile)

        nrand=215000
        if nrand > rdata.size:
            nrand=rdata.size
        #frac=0.1
        #nrand = int( rdata.size * frac )
        print("    plotting",nrand,"randoms")
        rind = eu.random.random_indices(rdata.size, nrand)

        plt=biggles.FramedPlot()
        plt.title=title
        plt.xlabel='RA'
        plt.ylabel='DEC'

        rpts=biggles.Points(rdata['ra'][rind], rdata['dec'][rind],
                            type='dot')
        pts=biggles.Points(data['ra'], data['dec'],
                           type='filled circle',
                           color='red',
                           size=0.5)
        # randoms go on bottom
        plt.add(rpts, pts)

        fpts, frpts=eu.plotting.fake_points(['filled circle']*2,
                                            ['random','lenses'],
                                            colors=['black','red'])
        key=biggles.PlotKey(0.9,0.9, [fpts, frpts], halign='right')
        plt.add(key)

        print("writing:",epsfile)
        plt.write_eps(epsfile)
        converter.convert(epsfile,
                          dpi=150,
                          verbose=True)


    def _write_html(self):
        """
        make an html file with all the z hist plots
        """
        from glob import glob
        d=get_match_weights_dir(self['lens_run'],
                                self['rand_run'],
                                self['bin_scheme'])

        print("making html files")

        command='cd %s; im2html -p *weights*.png > zhist.html' % d
        os.system(command)

        command='cd %s; im2html -p *radec*.png > radec.html' % d
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

        if 'weight_col' in self['rrun_conf']['lens_conf']:
            raise RuntimeError("weight_col is no longer supported")

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

        raise RuntimeError("adapt to no weight col")
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

        wr,weights = self._do_z_match(rz, z[w])
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

        wr,weights=self._do_z_match(rz[wr0], z[w])

        wr=wr0[wr]
        return w,wr,weights

    def _do_z_match(self, rz, z, weights=None):
        """
        all the other match methods call this one
        """

        import weighting

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
