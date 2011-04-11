from __future__ import print_function
import numpy
import os,sys
from sys import stdout

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

import lensing

import cosmology
import copy

import zphot

def create_input(sample):
    """
        create_input('01')
    """

    conf = lensing.files.json_read('scat', sample)
    if 'dr8regauss' in conf['catalog']:
        c = DR8Catalog(sample)
    else:
        c = DESMockSrcCatalog(sample)
    c.create_objshear_input()

class DR8Catalog(dict):
    """

    Before using, make sure you have matched the regauss cols with your chosen
    photoz sample using lensing.regauss.zphot_match()

    """
    def __init__(self, sample):
        conf = lensing.files.json_read('scat', sample)
        for k in conf:
            self[k] = conf[k]

        if 'dr8regauss' not in self['catalog']:
            raise ValueError("Expected dr8regauss as catalog")
        if self['sigmacrit_style'] != 2:
            raise ValueError("Expected sigmacrit_style 2")

        self.open_all_columns()

    def open_all_columns(self):
        print("opening regauss columns for procrun:",self['procrun'])
        self.rgcols = lensing.regauss.open_columns(self['procrun'], self['sweeptype'])
        print("  #rows:",self.rgcols['photoid'].size)
        print("opening zphot columns for pzrun:",self['pzrun'])
        self.pzcols = zphot.weighting.open_pofz_columns(self['pzrun'])
        print("  #rows:",self.pzcols['photoid'].size)

    def create_objshear_input(self):
        filter=self['filter']
        zlvals=lensing.sigmacrit.make_zlvals(self['dzl'], self['zlmin'], self['zlmax'])
        nzl = zlvals.size

        keep,zphot_matches = self.select()

        print("creating output struct")
        dt = lensing.files.scat_dtype(self['sigmacrit_style'], nzl=nzl)
        output = numpy.zeros(keep.size, dtype=dt)

        print("extracting basic columns and copying into output")
        output['ra'][:] = self.rgcols['ra'][keep]
        output['dec'][:] = self.rgcols['dec'][keep]
        #output['g1'][:] = self.rgcols['e1_rg_'+filter+'_eq'][keep]
        #output['g2'][:] = self.rgcols['e2_rg_'+filter+'_eq'][keep]
        output['g1'][:] = self.rgcols['e1_rg_'+filter][keep]
        output['g2'][:] = self.rgcols['e2_rg_'+filter][keep]
        output['err'][:] = self.rgcols['uncer_rg_'+filter][keep]
        
        return output
        raise ValueError("need to add equatorial rotated shapes")

    def select(self):
        """

        Apply selections to the sources, such as regauss processing flags,
        resolution cuts and selecting those with good photoz matches.

        """

        filter = self['filter']

        match_column = 'match_zphot%s' % self['pzrun']
        if match_column not in self.rgcols:
            raise ValueError("First use regauss.zphot_match() to match zphot and regauss")

        print("reading:",match_column)
        m = self.rgcols[match_column][:]
        print("Reading R")
        R = self.rgcols['R_rg_'+filter][:]
        print("Reading corrflags")
        flags = self.rgcols['corrflags_rg_'+filter][:]


        match_logic = (m >= 0)
        wmatch=where1(match_logic)
        print("Found %s/%s zphot matches" % (wmatch.size, m.size))
        

        print("Getting R logic")
        R_logic = R > self['Rcut']
        wR = where1(match_logic & R_logic)
        print("R > %f: %i/%i" % (self['Rcut'], wR.size, wmatch.size))


        print("Getting flag logic")
        flag_logic = flags == 0
        wflag = where1(match_logic & flag_logic)
        print("corrflags_rg == 0: %i/%i" % (wflag.size, wmatch.size))


        # combined logic
        wgood=where1(match_logic & R_logic & flag_logic)

        print("Found %s/%s good objects" % (wgood.size, m.size))

        wgood.sort()
        # wgood is index into overall set
        matches = m[wgood]
        return wgood, matches



    def calc_scinv(self, pzrun, show=False):
        """

        THIS NEEDS TO BE FINISHED

        - get the correction factor N(z)/sumpofz(z)
            calls zphot.weighting.pofz_correction(pzrun)

        - for each p(z) multiply by that factor
        - generate mean inverse critical density.


        """

        pzconf = zphot.cascade_config(pzrun)
        wo = zphot.weighting.WeightedOutputs()


class DESMockSrcCatalog(dict):
    """
    This reads the mock catalog and creates an input
    catalog for objshear


    Need to have config files that tell us what kind
    of sigma crit inv we are using.
    """

    def __init__(self, sample):
        conf = lensing.files.json_read('scat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['desmocks-2.13o-f8']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

        self['cosmo'] = lensing.files.json_read('cosmo',self['cosmo_sample'])

    def file(self):
        fname = lensing.files.sample_file('scat',self['sample'])
        return fname
    def read(self):
        return lensing.files.scat_read(self['sample'])

    def add_scinv(self):
        from . import sigmacrit
        import zphot
        if 'zlmin' not in self or 'zlmax' not in self or 'dzl' not in self:
            raise ValueError("You must have zlmin,zlmax,dzl in config")


        """
        note we don't really need this because we have perfect zs, but this
        is what you would do with a real sample

        wo = zphot.WeightedOutputs()
        z_file = wo.z_file(self['pzrun'], chunk=0)
        zscols = zphot.weighting.read_z(z_file)
        zsvals = (zscols['zmax']+zscols['zmin'])/2.0
        
        sc = sigmacrit.ScinvCalculator(zsvals, 
                                       self['dzl'],
                                       self['zlmin'],
                                       self['zlmax'])
        """

        # ok, now read the original file and add to it
        # the man scinv as a function of zlens
        # for now just using exact z

        zlvals = sigmacrit.make_zlvals(self['dzl'],
                                       self['zlmin'],
                                       self['zlmax'])
        n_zlens=zlvals.size
        
        data = self.read_original()
        dt=data.dtype.descr
        new_dtype = dt + [('scinv','f8',n_zlens)]
        out = numpy.zeros(data.size, dtype=new_dtype)
        eu.numpy_util.copy_fields(data, out)


        # actually not using sc since redshifts are exact
        cosmo=self['cosmo']
        c = cosmology.Cosmo(omega_m=cosmo['omega_m'], H0=cosmo['H0'])

        print('adding scinv to each')

        for i in xrange(data.size):
            if ((i+1) % 1000) == 0:
                print("%s/%s  %s%%" % (i+1,data.size,float(i+1)/data.size*100.))

            zs = data['z'][i]
            out['scinv'][i,:] = c.sigmacritinv(zlvals, zs)

        h=copy.deepcopy(cosmo)
        h['zlvals'] = zlvals

        outf=self.original_file(scinv=True)
        print("writing scinv file:",outf)

        eu.io.write(outf, out, header=h)

    def create_objshear_input(self):
        outfile = self.file()

        print("sigmacrit_style:",self['sigmacrit_style'])
        if self['sigmacrit_style'] == 2:
            data = self.read_original(scinv=True)
            print('creating output array')
            output = self.output_array(data.size, nzl=data['scinv'][0,:].size)
        else:
            data = self.read_original()
            print('creating output array')
            output = self.output_array(data.size)


        print('copying data')
        output['ra'] = data['ora']
        output['dec'] = data['odec']

        print('not changing sign of gamma1')
        output['g1'] = data['gamma1']

        output['g2'] = data['gamma2']
        output['err'] = 0.0
        output['hpixid'] = -9999

        if self['sigmacrit_style'] == 1:
            output['z'] = data['z']
            output['dc'] = -9999
        else:
            output['scinv'] = data['scinv']

        lensing.files.scat_write(self['sample'], output)


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def original_file(self, scinv=False):
        d = self.original_dir()
        f='%s-sources.rec' % self['catalog']
        if scinv:
            f=f.replace('sources','sources-scinv')
        infile = path_join(d, f)
        return infile

    def read_original(self, scinv=False):
        infile = self.original_file(scinv=scinv)
        if not os.path.exists(infile):
            raise ValueError("File not found: %s\n" % infile)

        data = eu.io.read(infile, lower=True, verbose=True)
        return data

    def output_array(self, num, nzl=None):
        dt = lensing.files.scat_dtype(self['sigmacrit_style'], nzl=nzl)
        output = numpy.zeros(num, dtype=dt)
        return output


