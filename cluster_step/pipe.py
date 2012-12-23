import os
from math import ceil
from numpy import where, sqrt, random, zeros, arange
from . import files
import gmix_image

SEXTRACTOR_GAL=1
SEXTRACTOR_STAR=2

class Pipe(dict):
    def __init__(self, **keys):
        """
        parameters
        ----------
        run:
            run id
        psfnum:
            psf number
        shnum:
            shear number
        ccd:
            ccd number

        version: optional
            version id
        """

        self._check_keys(**keys)

        for k in keys:
            self[k] = keys[k]

        conf=files.read_config(self['run'])

        for k in conf:
            self[k]=conf[k]
        
        self._load_data()

    def process(self):
        pass

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'ccd' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, ccd=")

    def _load_data(self):
        self.cat=files.read_cat(**self)
        image_orig, self.hdr=files.read_image(**self)
        self.seg, self.seg_hdr=files.read_image(ftype='seg', **self)

        skypersec=self.hdr['sky']
        exptime=self.hdr['exptime']

        # assume gain same for both amplifiers
        gain=self.hdr['gaina']
        read_noise=self.hdr['rdnoisea']

        # note unit mismatch? Both check out
        self['sky'] = skypersec*exptime/gain
        self['skyvar'] = self.hdr['sky'] + read_noise
        self['skysig'] = sqrt(self['skyvar'])
        self['ivar'] = 1.0/self['skyvar']

        self.image = image_orig.astype('f8') - self['sky']
        del image_orig

    def get_cutout(self, index, size=None, with_seg=False):
        if size is None:
            size=self['cutout_size']

        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        if with_seg:
            segcut=Cutout(self.seg, cen, size)
            return cutout, segcut

        return cutout

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)


    def get_stars(self, good=True):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """
        logic=self.cat['class']==SEXTRACTOR_STAR

        if good:
            logic = logic & (self.ares['whyflag']==0)

        w,=where(logic)
        return w

    def get_gals(self, good=True):
        """
        Get things labeled as stars.  if good, also trim
        to those with good adaptive moments
        """

        logic=self.cat['class']==SEXTRACTOR_GAL

        if good:
            logic = logic & (self.ares['whyflag']==0)
        w,=where(logic)
        return w

    def get_psf_stars(self, good=True):
        """
        Get things labeled as stars in the psf star mag range.  if good, also
        trim to those with good adaptive moments
        """
        star_logic = self.cat['class']==SEXTRACTOR_STAR
        mag_logic  = ( (self.cat['mag_auto_r'] > self['star_minmag'])
                       & (self.cat['mag_auto_r'] < self['star_maxmag']))

        logic = star_logic & mag_logic

        if good:
            logic = logic & (self.ares['whyflag']==0)

        w,=where(logic)
        return w




    def run_admom(self):
        import admom

        # seed for random numbers added for pixels not
        # in seg image
        random.seed(35)
        
        ares=self.get_admom_struct()

        for i in xrange(self.cat.size):
            c,cseg=self.get_cutout(i, with_seg=True)

            im=c.subimage
            seg=cseg.subimage

            self.zero_seg(im, seg, self.cat['id'][i])
            cen=c.subcen

            res = admom.admom(im, cen[0], cen[1],
                              sigsky=self['skysig'],
                              guess=4.0,
                              nsub=1)

            if res['whyflag'] != 0:
                #print res['whyflag'],'%.16g %.16g' % (cen[0]-res['wrow'], cen[1]-res['wcol'])
                print 'index: %4d flag: %s' % (i,res['whystr'])

            ares['row_range'][i] = c.row_range
            ares['col_range'][i] = c.col_range

            ares['wrow'][i] = c.row_range[0] + res['wrow']
            ares['wcol'][i] = c.col_range[0] + res['wcol']
            
            #for n in ares.dtype.names:
            for n in res:
                if n in ares.dtype.names and n not in ['row','col','wrow','wcol']:
                    ares[n][i] = res[n]

        w, = where(ares['whyflag'] != 0)
        print 'found %d/%d bad   %s' % (w.size, ares.size, w.size/float(ares.size))
        self.ares=ares

        files.write_fits_output(data=ares,
                                ftype='admom',
                                **self)


    def plot_admom_sizemag(self, show=False):
        """
        Plot the size-mag.  Show psf stars in a different color.
        """
        import biggles

        if not hasattr(self, 'ares'):
            raise ValueError("run_admom first")

        wstar=self.get_stars()
        wgal=self.get_gals()
        wpsf=self.get_psf_stars()


        plt=biggles.FramedPlot()
        sigma=sqrt( (self.ares['Irr']+self.ares['Icc'])/2 )
        mag=self.cat['mag_auto_r']
        psize=0.7
        pstar=biggles.Points(mag[wstar], sigma[wstar], 
                             type='filled circle', color='red',
                             size=psize)
        pgal=biggles.Points(mag[wgal], sigma[wgal], 
                            type='filled diamond', color='dark green',
                            size=psize)

        ppsf=biggles.Points(mag[wpsf], sigma[wpsf], 
                            type='filled circle', color='blue',
                            size=psize)

        pstar.label='star'
        pgal.label='gal'
        ppsf.label='psf star'

        key=biggles.PlotKey(0.95,0.92, [pstar,ppsf,pgal], halign='right')

        plt.add(pstar, ppsf, pgal, key)
        plt.xlabel='mag_auto_r'
        plt.ylabel=r'$\sigma_{AM}$ [pixels]'
        label='%s-p%s-s%s-%s' % (self['run'],self['psfnum'],self['shnum'],
                                 self['ccd'])
        plt.add(biggles.PlotLabel(0.05,0.92,label,halign='left'))

        plt.xrange=[16.5,25.5]
        plt.yrange=[0.75,8]
        if show:
            plt.show()

        epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)
        self._write_plot(plt, epsfile)



    def _write_plot(self, plt, epsfile):
        import converter
        d=os.path.dirname(epsfile)
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass
        print 'writing eps file:',epsfile
        plt.write_eps(epsfile)

        converter.convert(epsfile, dpi=100)



    def show_many_cutouts(self, indices, **keys):
        for i in indices:
            self.show_cutout(i, **keys)
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                return


    def show_cutout(self, index, size=None, with_seg=False,
                    zero_seg=False):
        import biggles
        import images


        minv=-2.5*self['skysig']
        if not with_seg:
            c=self.get_cutout(index, size=size)
            subimage=c.get_subimage()


            images.multiview(subimage, min=minv)
        else:
            c,segc=self.get_cutout(index, size=size, with_seg=True)

            subimage=c.get_subimage()
            segsub=segc.get_subimage()


            implt=images.view(subimage,show=False, min=minv)
            segplt=images.view(segsub,show=False)

            implt.title='image cutout'
            segplt.title='seg cutout'

            if zero_seg:
                self.zero_seg(subimage, segsub, self.cat['id'][index])
                zplt=images.view(subimage,show=False, min=minv)
                zplt.title='zerod'

                tab=biggles.Table(2,2)
                tab[0,0]=implt
                tab[0,1]=segplt
                tab[1,0]=zplt
            else:
                tab=biggles.Table(1,2)
                tab[0,0]=implt
                tab[0,1]=segplt

            tab.show()

    def get_admom_struct(self):
        sxdt=self.cat.dtype.descr
        dt= [('wrow','f8'),
             ('wcol','f8'),
             ('row_range','f8',2),
             ('col_range','f8',2),
             ('Irr','f8'),
             ('Irc','f8'),
             ('Icc','f8'),
             ('e1','f8'),
             ('e2','f8'),
             ('rho4','f8'),
             ('a4','f8'),
             ('s2','f8'),
             ('uncer','f8'),
             ('s2n','f8'),
             ('numiter','i4'),
             ('nsub','i4'),
             ('whyflag','i4'),
             ('whystr','S10'),
             ('shiftmax','f8')]

        dt=sxdt + dt
        data=zeros(self.cat.size, dtype=dt)
        for n in self.cat.dtype.names:
            data[n][:] = self.cat[n][:]

        return data

class Cutout:
    def __init__(self, image, cen, size):
        self.image=image
        self.cen=cen
        self.size=size
        self._make_cutout()

    def get_subimage(self):
        return self.subimage

    def get_subcen(self):
        return self.subcen

    def _make_cutout(self):
        sh=self.image.shape
        cen=self.cen
        size=self.size

        if (cen[0] < 0 or cen[1] < 0
                or cen[0] > (sh[0]-1)
                or cen[1] > (sh[1]-1) ):
            mess=("center [%s,%s] is out of "
                  "bounds of image: [%s,%s] ")
            mess=mess % (cen[0],cen[1],sh[0],sh[1])
            raise ValueError(mess)

        minrow=int(     cen[0]-size/2.-0.5)
        maxrow=int(ceil(cen[0]+size/2.+0.5))
        mincol=int(     cen[1]-size/2.-0.5)
        maxcol=int(ceil(cen[1]+size/2.+0.5))

        if minrow < 0:
            minrow=0
        if maxrow > (sh[0]-1):
            maxrow=sh[0]-1

        if mincol < 0:
            mincol=0
        if maxcol > (sh[1]-1):
            maxcol=sh[1]-1
        
        # note +1 for python slices
        self.subimage=self.image[minrow:maxrow+1, mincol:maxcol+1].copy()
        self.subcen=[cen[0]-minrow, cen[1]-mincol]

        self.row_range=[minrow,maxrow]
        self.col_range=[mincol,maxcol]

