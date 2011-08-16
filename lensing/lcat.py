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

import cosmology

def instantiate_sample(sample):
    conf = lensing.files.read_config('lcat',sample)
    if conf['catalog'] in ['redmapper','ProPer']:
        # proper is the old version
        return RedMapper(sample)
    elif conf['catalog'] == 'maxbcg-full':
        return MaxBCG(sample)
    elif conf['catalog'] in ['redmapper-random','maxbcg-random']:
        return SDSSRandom(sample)
    elif conf['catalog'] == 'desmocks-2.13':
        return DESMockLensCatalog(sample)
    else:
        raise ValueError("don't know about catalog %s" % conf['catalog'])

def create_input(sample):
    """
    e.g.  create_input('01')
    """

    c = instantiate_sample(sample)
    c.create_objshear_input()

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

        conf = lensing.files.read_config('lcat',sample)
        for k in conf:
            self[k] = conf[k]

        if self['sample'] != sample:
            raise ValueError("The config sample '%s' doesn't match input '%s'" % (self['sample'],sample))

    def read(self):
        return lensing.files.lcat_read(sample=self['sample'])

    def file(self):
        fname = lensing.files.sample_file('lcat',self['sample'])
        return fname



class SDSSRandom(LcatBase):
    """
    This is used for the new dr8 cluster catalog "Red Mapper".  

    Random points are generated from the tycho stomp, edges are 
    checked against the basic map.
    
    I've also used it with the old MaxBCG for testing but it does not use the
    right mask for that catalog.

    """
    def __init__(self, sample, **keys):

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

    def read_original(self):
        """
        For randoms, the lens input is the catalog
        """
        return lensing.files.lcat_read(sample=self['sample'])


    def load_stomp_maps(self):
        if not hasattr(self, 'basic_map'):
            self.basic_map = es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'])
            self.tycho_map = es_sdsspy.stomp_maps.load(self['mapname'],self['tycho_maptype'])

    def create_objshear_input(self):
        self.load_stomp_maps()


        fname = self.file()
        nrand = self['nrand']
        #nrand = 10000

        print("Generating",self['nrand'],"random points "
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

            print(" -> maskflags ", end='')
            maskflags = self.get_maskflags(ra,dec,z)

            wgood = es_sdsspy.stomp_maps.quad_check(maskflags)
            print(" -> good ones:",wgood.size)
            if wgood.size > 0:
                output['zindex'][n:n+wgood.size] = numpy.arange(n,n+wgood.size,dtype='i8')
                output['ra'][n:n+wgood.size] = ra[wgood]
                output['dec'][n:n+wgood.size] = dec[wgood]
                output['z'][n:n+wgood.size] = z[wgood]
                output['maskflags'][n:n+wgood.size] = maskflags[wgood]
                n += wgood.size

        lensing.files.lcat_write(self['sample'], output)

    def get_maskflags(self, ra, dec, z):
        """

        Run the stomp edge checking code. This uses the basic map, while the
        points themselves are generated from the tycho map

        """

        self.load_stomp_maps()

        # Da is in Mpc
        Da = self.cosmo.Da(0.0, z)

        # radius in *degrees*
        radius = self['rmax']/Da*180./PI
        
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
            d = lensing.files.sample_dir('lcat',self['sample'])
            d = os.path.join(d,'plots')
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = os.path.join(d, '%s-%s-coverage.eps' % ('lcat',self['sample']))
            print("Writing to eps file:",epsfile)
            plt.write_eps(epsfile)
        return plt

class RedMapper(LcatBase):
    def __init__(self, sample, **keys):

        LcatBase.__init__(self, sample, **keys)

        if self['catalog'] not in ['redmapper','ProPer']:
            raise ValueError("Don't know about catalog: '%s'" % self['catalog'])

        self['mapname'] = 'boss'
        self['maptype'] = 'basic'
        self['tycho_maptype'] = 'tycho'

    def create_objshear_input(self):
        
        # this will change when we get the actual red mapper catalog
        lambda_field = 'lambda_zred'
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
        lambda_logic = self.lambda_logic(data[lambda_field])

        # make sure in the tycho window and two adjacent quadrants
        # not hitting edge

        in_tycho, maskflags = self.get_maskflags(data['ra'], data['dec'], data[z_field])

        tycho_logic = (in_tycho == 1)
        quad_logic = es_sdsspy.stomp_maps.quad_logic(maskflags)

        good = where1(lambda_logic & tycho_logic & quad_logic)
        print("Finally kept: %d/%d" % (good.size,data.size))

        print('creating output array')
        output = output_array(good.size)

        print('copying data')
        output['zindex']    = zindex[good]
        output['ra']        = data['ra'][good]
        output['dec']       = data['dec'][good]
        output['z']         = data[z_field][good]
        output['maskflags'] = maskflags[good]
        lensing.files.lcat_write(self['sample'], output)

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
            self.basic_map = es_sdsspy.stomp_maps.load(self['mapname'],self['maptype'])
            self.tycho_map = es_sdsspy.stomp_maps.load(self['mapname'],self['tycho_maptype'])

    def get_maskflags(self, ra, dec, z):
        """

        Run the stomp edge checking code.

        """

        self.load_stomp_maps()

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

        
        print("    Ensuring in tycho window")
        in_tycho = self.tycho_map.Contains(ra, dec, "eq")

        w=where1(in_tycho == 1)
        print("    Keeping %d/%d" % (w.size,ra.size))

        print("    Getting basic maskflags...")
        maskflags = self.basic_map.Contains(ra, dec, "eq", radius)

        w=es_sdsspy.stomp_maps.quad_check(maskflags)
        print("    Keeping %d/%d" % (w.size,ra.size))

        in_tycho = numpy.array(in_tycho, dtype='i8')
        maskflags = numpy.array(maskflags, dtype='i8')

        return in_tycho, maskflags


    def original_dir(self):
        catdir = lensing.files.catalog_dir()
        d = path_join(catdir, self['catalog'])
        return str(d)

    def original_file(self):
        d = self.original_dir()
        f = 'dr8_proper_v3.2_lamgt10.fit'
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
            d = lensing.files.sample_dir('lcat',self['sample'])
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
            d = lensing.files.sample_dir('lcat',self['sample'])
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


    def create_objshear_input(self):
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
        lensing.files.lcat_write(self['sample'], output)

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
            d = lensing.files.sample_dir('lcat',self['sample'])
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
        fname = lensing.files.sample_file('lcat',self['sample'])
        return fname

    def read(self, split=None):
        return lensing.files.lcat_read(sample=self['sample'], split=split)

    def create_objshear_input(self):
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

        lensing.files.lcat_write(self['sample'], output)


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


