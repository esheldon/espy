"""

Classes for each catalog type, e.g. desmocks

Also create_input() function to make the input files for objshear

"""
from __future__ import print_function
import numpy
from  numpy import pi as PI
import os,sys
from sys import stdout

import lensing

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

import es_sdsspy
from es_sdsspy import stomp_maps

from . import binning

import cosmology

def instantiate_sample(sample, **keys):
    conf = lensing.files.read_config('lcat',sample)

    if conf['catalog'] in ['redmapper-random','maxbcg-random']:
        return SDSSRandom(sample, **keys)

    elif conf['catalog'][0:9] == 'redmapper':
        return RedMapper(sample, **keys)

    elif conf['catalog'] == 'maxbcg-full':
        return MaxBCG(sample)

    elif conf['catalog'] == 'desmocks-2.13':
        return DESMockLensCatalog(sample)

    elif conf['catalog'] == 'sdss-voids-01':
        return SDSSVoids(sample)
    elif conf['catalog'] in ['sdss-voids-rand-01']:
        return SDSSVoidsRandom(sample)
    else:
        raise ValueError("don't know about catalog %s" % conf['catalog'])

def create_input(sample, **keys):
    """
    e.g.  create_input('01')
    """

    c = instantiate_sample(sample, **keys)
    c.create_objshear_input(**keys)

def plot_coverage(sample):
    c = instantiate_sample(sample)
    c.plot_coverage()
    c.plot_coverage(region='ngc')
    c.plot_coverage(region='sgc')

def original_file(sample):
    c = instantiate_sample(sample)
    return c.original_file()

def read_original(sample):
    c = instantiate_sample(sample)
    return c.read_original()


def output_array(num):
    dt = lensing.files.lcat_dtype()
    output = numpy.zeros(num, dtype=dt)
    return output


class LcatBase(dict):
    def __init__(self, sample, **keys):

        for k in keys:
            self[k] = keys[k]

        conf = lensing.files.read_config('lcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

    def read(self):
        return lensing.files.lcat_read(sample=self['sample'])

    def file(self):
        fname = lensing.files.sample_file(type='lcat',sample=self['sample'], ext='dat')
        return fname

    def read_original(self):
        infile = self.original_file()
        stdout.write("Reading original catalog: %s\n" % infile)
        data = eu.io.read(infile, lower=True, ensure_native=True)
        return data

    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)
    def original_file(self, ext='fits'):
        d = self.original_dir()
        f='%s.%s' % (self['catalog'],ext)
        infile = path_join(d, f)
        return infile


class SDSSVoids(LcatBase):
    def __init__(self, sample, **keys):

        # this copies each key to self[key]
        LcatBase.__init__(self, sample, **keys)

    def convert2fits(self):
        import recfile
        # ID    ra        dec         redshift      radius
        dt=[('id','i8'),('ra','f8'),('dec','f8'),('z','f8'),('radius','f8')]
        
        f=self.original_file(ext='txt')
        fout=self.original_file()
        print("reading:",f)

        with recfile.Recfile(f, skiplines=1, delim=' ',dtype=dt) as robj:
            data = robj[:]

        print("writing:",fout)
        eu.io.write(fout, data)


    def create_objshear_input(self, **keys):

        data = self.read_original()
        orig_size = data.size
        zindex = numpy.arange(orig_size,dtype='i8')

        zmin = self['zmin']
        zmax = self['zmax']

        good=where1(  (data['z'] > zmin) & (data['z'] < zmax) )
        print("  z cut: %s/%s: %s" % (data.size-good.size,
                                      orig_size,
                                      (data.size-good.size)/float(orig_size)) )
        if good.size == 0:
            stop

        print("Actually trimming the bad z for speed")
        data = data[good]
        zindex = zindex[good]

        #maskflags = self.get_maskflags(data['ra'][0:5], data['dec'][0:5], data['z'][0:5])
        maskflags = self.get_maskflags(data['ra'], data['dec'], data['z'])
        quad_logic = es_sdsspy.stomp_maps.quad_logic(maskflags)

        good = where1(quad_logic)
        print("  quad mask cut: %s/%s: %s" % (data.size-good.size,
                                              orig_size,
                                              (data.size-good.size)/float(orig_size)) )

        if good.size == 0:
            stop

        print('creating output array')
        output = output_array(good.size)

        print('copying data')
        output['zindex']    = zindex[good]
        output['ra']        = data['ra'][good]
        output['dec']       = data['dec'][good]
        output['z']         = data['z'][good]
        output['maskflags'] = maskflags[good]
        lensing.files.lcat_write(sample=self['sample'], data=output)

    def get_maskflags(self, ra, dec, z):
        """

        Run the stomp edge checking code.

        There are two checks: quadrant checking at full rmax, and a hard edge
        cut at rmax_hard.

        """

        # get radius for edge check
        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        print(cconf)
        
        c = cosmology.Cosmo(H0=cconf['H0'], omega_m=cconf['omega_m'])

        print("    Getting radii. z range of inputs is [%.2f, %.2f]" % (z.min(), z.max()))
        # Da is in Mpc
        Da = c.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI

        self.basic_map = \
            es_sdsspy.stomp_maps.load('boss','basic', maxres=2048)

        print("    Getting basic quadrant maskflags at full rmax: %0.2f" % self['rmax'])
        print("    radii are in range [%f,%f]" % (radius.min(), radius.max()))
        maskflags = self.basic_map.Contains(ra, dec, "eq", radius)
        w=es_sdsspy.stomp_maps.quad_check(maskflags)
        print("    Keeping %d/%d quad" % (w.size,ra.size))

        maskflags = numpy.array(maskflags, dtype='i8')
        return maskflags

    def plot_z_radius(self):
        import biggles
        import converter
        data=self.read_original()
        plt=eu.plotting.bscatter(data['z'],data['radius'],show=False)

        zb = binning.ZBinner(4)
        ll, hl = zb.bin_ranges()

        c1=biggles.Curve([ll[0]]*2, [6,60])
        c2=biggles.Curve([ll[1]]*2, [6,60])
        c3=biggles.Curve([ll[2]]*2, [6,60])
        c4=biggles.Curve([ll[3]]*2, [6,60])
        c5=biggles.Curve([hl[3]]*2, [6,60])

        plt.add(c1,c2,c3,c4,c5)

        w1=where1( (data['z'] > ll[0]) & (data['z'] < hl[0]) )
        w2=where1( (data['z'] > ll[1]) & (data['z'] < hl[1]) )
        w3=where1( (data['z'] > ll[2]) & (data['z'] < hl[2]) )
        w4=where1( (data['z'] > ll[3]) & (data['z'] < hl[3]) )

        type='filled circle'
        size=1
        plt.add(biggles.Points(data['z'][w1],data['radius'][w1],color='blue',type=type,size=size))
        plt.add(biggles.Points(data['z'][w2],data['radius'][w2],color='red',type=type,size=size))
        plt.add(biggles.Points(data['z'][w3],data['radius'][w3],color='magenta',type=type,size=size))
        plt.add(biggles.Points(data['z'][w4],data['radius'][w4],color='brown',type=type,size=size))
        plt.xlabel='z'
        plt.ylabel='radius'
        #plt.show()

        f=self.original_file().replace('.fits','-z-rad.eps')
        print(f)
        plt.write_eps(f)
        converter.convert(f, dpi=100, verbose=True)
        

class SDSSRandom(LcatBase):
    """
    This is used for the new dr8 cluster catalog "Red Mapper".  

    Random points are generated from the tycho stomp, edges are 
    checked against the basic map.
    
    I've also used it with the old MaxBCG for testing but it does not use the
    right mask for that catalog.

    """
    def __init__(self, sample, **keys):

        # this copies each key to self[key]
        LcatBase.__init__(self, sample, **keys)

        if self['catalog'] not in ['maxbcg-random','redmapper-random']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        self.cosmo = cosmology.Cosmo(omega_m=cconf['omega_m'], H0=cconf['H0'])

        self.zgen = eu.random.Generator(self.cosmo.dV, 
                                        xrange=[self['zmin'],self['zmax']], 
                                        nx=1000, 
                                        method='cut')
        self['mapname'] = 'boss'
        self['maptype'] = 'basic'
        self['tycho_maptype'] = 'tycho'
        # we might want to increase this!
        self['maxres'] = 2048

    def read_original(self):
        """
        For randoms, the lens input is the catalog
        """
        return lensing.files.lcat_read(sample=self['sample'])


    def load_stomp_maps(self):
        if not hasattr(self, 'basic_map'):
            self.basic_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'], 
                                          maxres=self['maxres'])
            self.tycho_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['tycho_maptype'],
                                          maxres=self['maxres'])

    def create_objshear_input(self, nrand=None, extra=None):
        """
        To work in chunks, send nrand= and extra=chunknum
        """
        if nrand is not None:
            if extra is None:
                raise ValueError("If sending nrand, also send extra=")
        else:
            nrand = self['nrand']

        self.load_stomp_maps()

        strict_edgecut = self.get('strict_edgecut',False)

        fname = self.file()

        print("Generating",nrand,"random points "
              "with z in [%0.2f,%0.2f]" % (self['zmin'],self['zmax']))
        n=0

        dt = lensing.files.lcat_dtype()
        output = numpy.zeros(nrand, dtype=dt)
        while n < nrand:
            if n == 0:
                print("  generating",nrand-n," ",end='')
            else:
                print("  re-generating",nrand-n," ",end='')
            print("z",end='')
            z = self.zgen.genrand(nrand-n)
            print(" -> ra,dec ",end='')
            ra,dec = self.tycho_map.GenerateRandomEq(nrand-n)

            if 'rmax_hard' in self:
                print(" -> maskflags_hard (%0.1f)" % self['rmax_hard'], end='')
                maskflags_hard = self.get_maskflags(ra,dec,z,hard=True)
                hard_edge_logic = es_sdsspy.stomp_maps.quad_logic(maskflags_hard, strict=True)

            print(" -> maskflags (%0.1f)" % self['rmax'], end='')
            maskflags = self.get_maskflags(ra,dec,z)
            quad_logic      = es_sdsspy.stomp_maps.quad_logic(maskflags,      strict=strict_edgecut)

            #wgood = es_sdsspy.stomp_maps.quad_check(maskflags, strict=strict_edgecut)
            if 'rmax_hard' in self:
                wgood = where1(quad_logic & hard_edge_logic)
            else:
                wgood = where1(quad_logic)

            print(" -> good ones:",wgood.size)
            if wgood.size > 0:
                output['zindex'][n:n+wgood.size] = numpy.arange(n,n+wgood.size,dtype='i8')
                output['ra'][n:n+wgood.size] = ra[wgood]
                output['dec'][n:n+wgood.size] = dec[wgood]
                output['z'][n:n+wgood.size] = z[wgood]
                output['maskflags'][n:n+wgood.size] = maskflags[wgood]
                n += wgood.size

        lensing.files.lcat_write(sample=self['sample'], data=output, extra=extra)

    def get_maskflags(self, ra, dec, z, hard=False):
        """

        Run the stomp edge checking code. This uses the basic map, while the
        points themselves are generated from the tycho map

        """

        self.load_stomp_maps()

        # Da is in Mpc
        Da = self.cosmo.Da(0.0, z)

        # radius in *degrees*
        if hard:
            rmax = self['rmax_hard']
        else:
            rmax = self['rmax']
        radius = rmax/Da*180./PI
        
        maskflags = self.basic_map.Contains(ra, dec, "eq", radius)

        return numpy.array(maskflags, dtype='i8')

    def plot_coverage(self, region='both', show=True, dops=True):
        """

        Plot a random subset of the randoms along with the
        boss survey geometry area as bounding boxes

        """
        import biggles
        from biggles import FramedPlot, Points, PlotKey

        symsize=0.5

        l = self.read()
        llam,leta = eu.coords.eq2sdss(l['ra'],l['dec'])

        lammin,lammax = (-70.,70.)

        if region=='ngc':
            # a bit high to make room for the legend
            emin,emax=(-40,60)

            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
            width=1.5

        elif region=='sgc':
            emin,emax=(125,165)
            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
            width=2
        else:
            emin,emax=(-40,165)
            biggles.configure('screen','width', 1140)
            biggles.configure('screen','height', 1140)
            width=2

        wl=where1((leta > emin) & (leta < emax))
        llam=llam[wl]
        leta=leta[wl]


        plt=FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'
        xrng = (lammin, lammax)
        yrng = (emin, emax)
        plt.xrange = xrng
        plt.yrange = yrng


        print("adding random subset of randoms")

        ii = eu.numpy_util.random_subset(llam.size, 500000)
        allp = Points(llam[ii],leta[ii],type='dot',size=symsize)
        plt.add(allp)

        plt.aspect_ratio = (yrng[1]-yrng[0])/float(xrng[1]-xrng[0])

        #es_sdsspy.stomp_maps.plot_boss_geometry(color='blue',plt=plt,show=False)
        es_sdsspy.stomp_maps.plot_boss_geometry(plt=plt,show=False,width=width)

        if show:
            plt.show()

        if dops:
            d = lensing.files.sample_dir(type='lcat',sample=self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, '%s-%s-coverage.eps' % ('lcat',self['sample']))
            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
        return plt


class SDSSVoidsRandom(LcatBase):
    """
    This is used for the randoms generated by Peter Sutter for the
    voids.  It is ~DR7
    """
    def __init__(self, sample, **keys):

        # this copies each key to self[key]
        LcatBase.__init__(self, sample, **keys)

        if self['catalog'] not in ['sdss-voids-rand-01']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        self.cosmo = cosmology.Cosmo(omega_m=cconf['omega_m'], H0=cconf['H0'])

        self.zgen = eu.random.Generator(self.cosmo.dV, 
                                        xrange=[self['zmin'],self['zmax']], 
                                        nx=1000, 
                                        method='cut')
        self['mapname'] = 'boss'
        self['maptype'] = 'basic'
        self['tycho_maptype'] = 'tycho'
        # we might want to increase this!
        self['maxres'] = 2048

    def load_stomp_maps(self):
        if not hasattr(self, 'basic_map'):
            self.basic_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'], 
                                          maxres=self['maxres'])
            self.tycho_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['tycho_maptype'],
                                          maxres=self['maxres'])
    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)

    def original_file(self):
        d = self.original_dir()
        f='%s.fits' % self['catalog']
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        stdout.write("Reading original catalog: %s\n" % infile)
        data = eu.io.read(infile, lower=True, ensure_native=True)
        return data


    def create_objshear_input(self, nsplit=None, split=None):
        """
        To work in chunks, send nrand= and extra=chunknum
        """

        if nsplit is None or split is None:
            raise ValueError("send nsplit= and split=")
        print("doing split %s: %s/%s" % (split,split+1,nsplit))

        self.load_stomp_maps()

        strict_edgecut = self['strict_edgecut']

        fname = self.file()

        n=0

        data=self.read_original()
        zindex = numpy.arange(data.size,dtype='i8')

        # do in chunks so we can see the progress
        npersplit = data.size/nsplit
        nleft = data.size % nsplit

        data = data[split*npersplit:(split+1)*npersplit]
        zindex = zindex[split*npersplit:(split+1)*npersplit]

        print("Generating z in [%0.2f,%0.2f]" % (self['zmin'],self['zmax']))
        z = self.zgen.genrand(data.size)


        print(" -> maskflags, max radius: %0.1f" % self['rmax'])
        maskflags = self.get_maskflags(data['ra'],data['dec'],z)

        quad_logic = es_sdsspy.stomp_maps.quad_logic(maskflags, strict=strict_edgecut)

        wgood = where1(quad_logic)

        print(" -> good ones:",wgood.size)

        data      = data[wgood]
        zindex    = zindex[wgood]
        z         = z[wgood]
        maskflags = maskflags[wgood]

        dt = lensing.files.lcat_dtype()
        output = numpy.zeros(wgood.size, dtype=dt)
        output['zindex'][:]    = zindex
        output['ra'][:]        = data['ra']
        output['dec'][:]       = data['dec']
        output['z'][:]         = z
        output['maskflags'][:] = maskflags

        lensing.files.lcat_write(sample=self['sample'], data=output, split=split)


    def get_maskflags(self, ra, dec, z, hard=False):
        """

        Run the stomp edge checking code. This uses the basic map, while the
        points themselves are generated from the tycho map

        """

        self.load_stomp_maps()

        # Da is in Mpc
        Da = self.cosmo.Da(0.0, z)

        # radius in *degrees*
        if hard:
            rmax = self['rmax_hard']
        else:
            rmax = self['rmax']
        radius = rmax/Da*180./PI
        
        maskflags = self.basic_map.Contains(ra, dec, "eq", radius)

        return numpy.array(maskflags, dtype='i8')



class RedMapper(LcatBase):
    def __init__(self, sample, **keys):

        # this copies each key to self[key]
        LcatBase.__init__(self, sample, **keys)

        if self['catalog'] not in ['redmapper-dr8-3.4-like',
                                   'redmapper-dr8-3.4-nord',
                                   'redmapper-dr8-3.14',
                                   'redmapper-dr8-3.14-cen2']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        self['mapname'] = 'boss'
        self['maptype'] = 'basic'
        self['tycho_maptype'] = 'tycho'
        # we might want to increase this!
        self['maxres'] = 2048


    def create_objshear_input(self, **keys):
        
        strict_edgecut = self.get('strict_edgecut',False)

        z_field = 'z_lambda'

        fname = self.file()

        data = self.read_original()

        # keep index into original data
        orig_size = data.size
        zindex = numpy.arange(orig_size,dtype='i8')

        # trim z for speed
        z_logic = self.z_logic(data[z_field])

        print("Actually trimming the bad z for speed")
        w=where1(z_logic)
        data = data[w]
        zindex = zindex[w]


        # trim poorly understood low lambda stuff
        lambda_logic = self.lambda_logic(data['lambda_chisq'])

        # make sure in the tycho window and two adjacent quadrants
        # not hitting edge (or no edge if strict=True)

        md = self.get_maskflags(data['ra'], data['dec'], data[z_field])

        tycho_logic = (md['in_tycho'] == 1)
        quad_logic = es_sdsspy.stomp_maps.quad_logic(md['maskflags'], 
                                                     strict=strict_edgecut)
        hard_edge_logic = es_sdsspy.stomp_maps.quad_logic(md['maskflags_hard'], 
                                                          strict=True)

        good = where1(lambda_logic & tycho_logic & quad_logic & hard_edge_logic)
        print("Finally kept: %d/%d" % (good.size,data.size))

        print('creating output array')
        output = output_array(good.size)

        print('copying data')
        output['zindex']    = zindex[good]
        output['ra']        = data['ra'][good]
        output['dec']       = data['dec'][good]
        output['z']         = data[z_field][good]
        output['maskflags'] = md['maskflags'][good]
        lensing.files.lcat_write(sample=self['sample'], data=output)

    def lambda_logic(self, lam):
        print("Cutting lambda > %0.2f" % self['lambda_min'])
        logic = lam > self['lambda_min']

        w=where1(logic)
        print("    Keeping %d/%d" % (w.size,lam.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic


    def z_logic(self, z):
        print("Cutting z to [%f, %f]" % (self['zmin'],self['zmax']))
        logic = (z > self['zmin']) & (z < self['zmax']) 

        w=where1(logic)
        print("    Keeping %d/%d" % (w.size,z.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic

    def load_stomp_maps(self):
        if not hasattr(self, 'basic_map'):
            self.basic_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'], 
                                          maxres=self['maxres'])
            self.tycho_map = \
                es_sdsspy.stomp_maps.load(self['mapname'],self['tycho_maptype'],
                                          maxres=self['maxres'])

    def get_maskflags(self, ra, dec, z, types=['tycho','rmax','rmax_hard']):
        """

        Run the stomp edge checking code.

        There are two checks: quadrant checking at full rmax, and a hard edge
        cut at rmax_hard.

        """

        self.load_stomp_maps()

        # get radius for edge check
        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        print(cconf)
        
        c = cosmology.Cosmo(H0=cconf['H0'], omega_m=cconf['omega_m'])

        print("    Getting radii. z range of inputs is [%.2f, %.2f]" % (z.min(), z.max()))
        # Da is in Mpc
        Da = c.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI
        radius_hard = self['rmax_hard']/Da*180./PI

        in_tycho=None
        maskflags_hard=None
        maskflags=None
        
        if 'tycho' in types:
            print("    Ensuring in tycho window")
            in_tycho = self.tycho_map.Contains(ra, dec, "eq")

            w=where1(in_tycho == 1)
            print("    Keeping %d/%d" % (w.size,ra.size))
            in_tycho = numpy.array(in_tycho, dtype='i8')

        if 'rmax_hard' in types:
            # two types of edge checks
            print("    Getting basic maskflags at rmax_hard: %0.2f" % self['rmax_hard']) 
            print("    radii are in range [%f,%f]" % (radius_hard.min(), radius_hard.max()))
            maskflags_hard = self.basic_map.Contains(ra, dec, "eq", radius_hard)
            whard=es_sdsspy.stomp_maps.quad_check(maskflags_hard, strict=True)
            print("    Keeping %d/%d from hard cut" % (whard.size,ra.size))
            maskflags_hard = numpy.array(maskflags_hard, dtype='i8')

        if 'rmax' in types:
            print("    Getting basic quadrant maskflags at full rmax: %0.2f" % self['rmax'])
            print("    radii are in range [%f,%f]" % (radius.min(), radius.max()))
            maskflags = self.basic_map.Contains(ra, dec, "eq", radius)
            w=es_sdsspy.stomp_maps.quad_check(maskflags)
            print("    Keeping %d/%d quad" % (w.size,ra.size))

            maskflags = numpy.array(maskflags, dtype='i8')

        out={}
        if in_tycho is not None:
            out['in_tycho'] = in_tycho
        if maskflags_hard is not None:
            out['maskflags_hard'] = maskflags_hard
        if maskflags is not None:
            out['maskflags'] = maskflags

        return out


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)

    def original_file(self):
        d = self.original_dir()
        f='%s.fits' % self['catalog']
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        stdout.write("Reading original catalog: %s\n" % infile)
        data = eu.io.read(infile, lower=True, ensure_native=True)
        return data


    def plot_coverage_bybin(self, binner, region='both', show=True, dops=True, rand=None):
        import pcolors
        import biggles
        import converter
        from biggles import FramedPlot, Points, PlotKey


        orig = self.read_original()
        lcat = self.read()

        all_clam,all_ceta = eu.coords.eq2sdss(orig['ra'],orig['dec'])

        l = orig[lcat['zindex']]
        clam,ceta = eu.coords.eq2sdss(lcat['ra'],lcat['dec'])

        clammin,clammax = (-70.,120.)

        if region=='ngc':
            # a bit high to make room for the legend
            emin,emax=(-40,60)

            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)

            clammin,clammax = (-70.,120.)

        elif region=='sgc':
            emin,emax=(105,165)
            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
            clammin,clammax = (-50.,90.)
        else:
            emin,emax=(-40,165)
            biggles.configure('screen','width', 1140)
            biggles.configure('screen','height', 1140)
            clammin,clammax = (-70.,120.)

        wl=where1((all_ceta > emin) & (all_ceta < emax))
        all_clam=all_clam[wl]
        all_ceta=all_ceta[wl]

        wl=where1((ceta > emin) & (ceta < emax))
        clam=clam[wl]
        ceta=ceta[wl]
        l=l[wl]


        plt=FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'
        xrng = (clammin, clammax)
        yrng = (emin, emax)
        plt.xrange = xrng
        plt.yrange = yrng


        print("adding all lenses")

        type = 'filled circle'
        symsize=0.2
        colors = pcolors.rainbow(binner['nbin'],'hex')


        if rand is not None:
            clam_r,ceta_r = eu.coords.eq2sdss(rand['ra'],rand['dec'])
            wl=where1((ceta_r > emin) & (ceta_r < emax))
            clam_r=clam_r[wl]
            ceta_r=ceta_r[wl]
            rp = Points(clam_r, ceta_r, type='dot', size=0.2)
            plt.add(rp)
        
        size_min=0.2
        size_max=4

        sizes=[]
        minlambda = l['lambda_zred'].min()
        maxlambda = l['lambda_zred'].max()
        for i in xrange(binner['nbin']):
            w=binner.select_bin(l, i)
            mlam=l['lambda_zred'][w].mean()
            # scale 0 to 1
            sz=(mlam-minlambda)/maxlambda
            # now scale size
            sz = size_min + sz*(size_max-size_min)
            sizes.append(sz)

        all_plots=[]
        labels=[]
        #for i in xrange(binner['nbin']):
        for i in reversed(xrange(binner['nbin'])):
            w=binner.select_bin(l, i)

            #points = Points(clam[w], ceta[w],type=type,size=symsize, color=colors[i])
            points = Points(clam[w], ceta[w],type=type,size=sizes[i], color=colors[i])
            labels.append(binner.bin_label(i))

            plt.add(points)


        labels.reverse()
        fakepoints = eu.plotting.fake_filled_circles(labels, colors)
        key=PlotKey(0.95,0.95,fakepoints,halign='right',size=1.5)
        plt.add(key)

        plt.aspect_ratio = (yrng[1]-yrng[0])/float(xrng[1]-xrng[0])

        es_sdsspy.stomp_maps.plot_boss_geometry(color='blue',plt=plt,show=False)

        if show:
            plt.show()

        if dops:
            d = lensing.files.sample_dir(type='lcat',sample=self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, 'lcat-%s-coverage-bybin.eps' % self['sample'])
            if rand is not None:
                epsfile=epsfile.replace('.eps','-withrand.eps')
            if region in ['sgc','ngc']:
                epsfile=epsfile.replace('.eps','-%s.eps' % region)

            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
            print("converting to png")
            converter.convert(epsfile, dpi=300)
        return plt




    def plot_coverage(self, region='both', show=True, dops=True):
        import biggles
        from biggles import FramedPlot, Points, PlotKey


        l = self.read()
        llam,leta = eu.coords.eq2sdss(l['ra'],l['dec'])
        maskflags = l['maskflags']

        lammin,lammax = (-70.,70.)

        if region=='ngc':
            # a bit high to make room for the legend
            emin,emax=(-40,60)

            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)

        elif region=='sgc':
            emin,emax=(100,165)
            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
        else:
            emin,emax=(-40,165)
            biggles.configure('screen','width', 1140)
            biggles.configure('screen','height', 1140)

        wl=where1((leta > emin) & (leta < emax))
        llam=llam[wl]
        leta=leta[wl]
        maskflags=maskflags[wl]


        plt=FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'
        xrng = (lammin, lammax)
        yrng = (emin, emax)
        plt.xrange = xrng
        plt.yrange = yrng


        print("adding all lenses")

        type = 'filled circle'
        symsize=0.2

        allp = Points(llam,leta,type=type,size=symsize)
        allp.label='all'
        plt.add(allp)

        wquad = es_sdsspy.stomp_maps.quad_check(maskflags)
        print("adding quad pass")
        quadp = Points(llam[wquad],leta[wquad],type=type,color='red',size=symsize)
        quadp.label = 'quad good'
        plt.add(quadp)

        fakepoints = eu.plotting.fake_filled_circles(['all','quad good'],['black','red'])
        key=PlotKey(0.95,0.95,fakepoints,halign='right')
        plt.add(key)

        plt.aspect_ratio = (yrng[1]-yrng[0])/float(xrng[1]-xrng[0])

        es_sdsspy.stomp_maps.plot_boss_geometry(color='blue',plt=plt,show=False)

        if show:
            plt.show()

        if dops:
            d = lensing.files.sample_dir(type='lcat',sample=self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, 'lcat-%s-coverage.eps' % self['sample'])
            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
        return plt





class MaxBCG(LcatBase):
    def __init__(self, sample, **keys):

        LcatBase.__init__(self, sample, **keys)

        if self['catalog'] not in ['maxbcg-full']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])


    def create_objshear_input(self, **keys):
        fname = self.file()

        data = self.read_original()

        # keep index into original data
        orig_size = data.size
        zindex = numpy.arange(orig_size,dtype='i8')

        # trim not well understood low ngals stuff
        ngals_logic = self.ngals_logic(data['ngals_r200'])

        # trim z for speed
        z_logic = self.z_logic(data['photoz_cts'])

        good = where1(ngals_logic & z_logic)
        print("Finally kept: %d/%d" % (good.size,data.size))

        print('creating output array')
        output = output_array(good.size)

        print('copying data')
        output['zindex']    = zindex[good]
        output['ra']        = data['ra'][good]
        output['dec']       = data['dec'][good]
        output['z']         = data['photoz_cts'][good]
        output['maskflags'] = self.get_maskflags(output['ra'],output['dec'],output['z'])
        lensing.files.lcat_write(sample=self['sample'], data=output)

    def ngals_logic(self, ngals):
        print("Cutting ngals >= %d" % self['ngals_r200_min'])
        logic = ngals >= self['ngals_r200_min']

        w=where1(logic)
        print("    Keeping %d/%d" % (w.size,ngals.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic


    def z_logic(self, z):
        print("Cutting z to [%f, %f]" % (self['zmin'],self['zmax']))
        logic = (z > self['zmin']) & (z < self['zmax']) 

        w=where1(logic)
        print("    Keeping %d/%d" % (w.size,z.size))
        if w.size == 0:
            raise ValueError("No objects passed z cut")

        return logic

    def get_maskflags(self, ra, dec, z):
        """

        Run the stomp edge checking code.

        """

        print("Getting BOSS basic maskflags")
        # get radius for edge check
        cconf = lensing.files.read_config('cosmo',self['cosmo_sample'])
        print(cconf)
        
        c = cosmology.Cosmo(H0=cconf['H0'], omega_m=cconf['omega_m'])

        print("    Getting radii. z range of inputs is [%.2f, %.2f]" % (z.min(), z.max()))
        # Da is in Mpc
        Da = c.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI

        print("    radii are in range [%f,%f]" % (radius.min(), radius.max()))

        
        map = stomp_maps.load('boss','basic')
        print("    Getting maskflags...")
        maskflags = map.Contains(ra, dec, "eq", radius)

        return numpy.array(maskflags, dtype='i8')


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)

    def original_file(self):
        d = self.original_dir()
        f='catalog_full_bcgs_orig.fit'
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        stdout.write("Reading original catalog: %s\n" % infile)
        data = eu.io.read(infile, lower=True, ensure_native=True)
        return data

    def plot_coverage(self, region='both', show=True, dops=True):
        """

        Plot the lenses, both masked and unmasked plus goemetry
        bounds

        """
        import biggles
        from biggles import FramedPlot, Points, PlotKey

        symsize=0.5

        l = self.read()
        llam,leta = eu.coords.eq2sdss(l['ra'],l['dec'])
        maskflags = l['maskflags']

        lammin,lammax = (-70.,70.)

        if region=='ngc':
            # a bit high to make room for the legend
            emin,emax=(-40,60)

            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)

        elif region=='sgc':
            emin,emax=(125,165)
            biggles.configure('screen','width', 1800)
            biggles.configure('screen','height', 1140)
        else:
            emin,emax=(-40,165)
            biggles.configure('screen','width', 1140)
            biggles.configure('screen','height', 1140)

        wl=where1((leta > emin) & (leta < emax))
        llam=llam[wl]
        leta=leta[wl]
        maskflags=maskflags[wl]


        plt=FramedPlot()
        plt.xlabel=r'$\lambda$'
        plt.ylabel=r'$\eta$'
        xrng = (lammin, lammax)
        yrng = (emin, emax)
        plt.xrange = xrng
        plt.yrange = yrng


        print("adding all lenses")

        allp = Points(llam,leta,type='dot',size=symsize)
        allp.label='all MaxBCG'
        plt.add(allp)

        wquad = es_sdsspy.stomp_maps.quad_check(maskflags)
        print("adding quad pass")
        quadp = Points(llam[wquad],leta[wquad],type='dot',color='red',size=symsize)
        quadp.label = 'quad good'
        plt.add(quadp)

        fakepoints = eu.plotting.fake_filled_circles(['all MaxBCG','quad good'],['black','red'])
        key=PlotKey(0.95,0.95,fakepoints,halign='right')
        plt.add(key)

        plt.aspect_ratio = (yrng[1]-yrng[0])/float(xrng[1]-xrng[0])

        es_sdsspy.stomp_maps.plot_boss_geometry(color='blue',plt=plt,show=False)

        if show:
            plt.show()

        if dops:
            d = lensing.files.sample_dir(type='lcat',sample=self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, 'lcat-%s-coverage.eps' % self['sample'])
            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
        return plt





class DESMockLensCatalog(dict):
    """
    Provides the interface needed by CatalogConverter
    """

    def __init__(self, sample, **keys):
        conf = lensing.files.read_config('lcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['catalog'] not in ['desmocks-2.13']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

    def file(self):
        fname = lensing.files.sample_file(type='lcat',sample=self['sample'])
        return fname

    def read(self, split=None):
        return lensing.files.lcat_read(sample=self['sample'], split=split)

    def create_objshear_input(self, **keys):
        fname = self.file()

        data = self.read_original()

        print('creating output array')
        output = output_array(data.size)

        print('copying data')
        output['zindex'] = numpy.arange(data.size,dtype='i8')
        output['ra'] = data['ra']
        output['dec'] = data['dec']
        output['z'] = data['z']
        #output['dc'] = -9999.0

        lensing.files.lcat_write(sample=self['sample'], data=output)


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return d

    def original_file(self):
        d = self.original_dir()
        f='%s-halos.fit' % self['catalog']
        infile = path_join(d, f)
        return infile

    def read_original(self):
        infile = self.original_file()
        if not os.path.exists(infile):
            raise ValueError("File not found: %s\n" % infile)

        stdout.write("Reading lens catalog: %s\n" % infile)
        data = eu.io.read(infile, lower=True, verbose=True, 
                              ensure_native=True)
        return data


