from math import ceil
from numpy import where, sqrt, random
from . import files
import gmix_image

SEXTRACTOR_GAL=1
SEXTRACTOR_STAR=2
DEFAULT_SIZE=25

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
        repeat:
            repeat number

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

    def _check_keys(self, **keys):
        if ('run' not in keys
                or 'psfnum' not in keys
                or 'shnum' not in keys
                or 'repeat' not in keys):
            raise ValueError("send run=, psfnum=, "
                             "shnum=, repeat=")

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

    def get_cutout(self, index, size=DEFAULT_SIZE, with_seg=False):
        cen=[self.cat['row'][index], self.cat['col'][index]]

        cutout=Cutout(self.image, cen, size)

        if with_seg:
            segcut=Cutout(self.seg, cen, size)
            return cutout, segcut

        return cutout

    def show_many_cutouts(self, indices, **keys):
        for i in indices:
            self.show_cutout(i, **keys)
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                return

    def zero_seg(self, im, seg, id):
        w=where( (seg != id) & (seg != 0) )
        if w[0].size != 0:
            im[w] = self['skysig']*random.randn(w[0].size)

    def show_cutout(self, index, size=DEFAULT_SIZE, with_seg=False,
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

    def get_stars(self):
        w,=where(self.cat['class']==SEXTRACTOR_STAR)
        return w

    def get_gals(self):
        w,=where(self.cat['class']==SEXTRACTOR_GAL)
        return w


    def run_admom(self):
        import admom
        for i in xrange(self.cat.size):
            c,cseg=self.get_cutout(i, with_seg=True)

            im=c.subimage
            seg=cseg.subimage

            self.zero_seg(im, seg, self.cat['id'][i])
            cen=c.subcen
            ares = admom.admom(im,
                                    cen[0],
                                    cen[1],
                                    sigsky=self['skysig'],
                                    guess=2.0,
                                    nsub=1)

            print ares['whyflag'],'%.16g %.16g' % (cen[0]-ares['wrow'], cen[1]-ares['wcol'])

    def get_admom_struct(self):
        dt= ['row':row,
             'col':col,
             'Irr':Irr,
             'Irc':Irc,
             'Icc':Icc,
             'e1':e1,
             'e2':e2,
             'rho4':rho4,
             'a4':a4, 
             's2':s2,
             'uncer':uncer,
             's2n':s2n,
             'numiter':numiter,
             'wrow':wrow,
             'wcol':wcol,
             'nsub':nsub,
             'whyflag':whyflag,
             'whystr':whystr,
             'shiftmax':shiftmax,
             'sky':sky,
             'sigsky':sigsky}


    """
    def plot_sizemag(self):
        import biggles

        wstar=self.get_stars()
        wgal=self.get_gals()

        plt=biggles.FramedPlot()
        pstar=biggles.Points(self.cat['mag_auto_r'][wstar],
                             self.cat['
    """

    def process(self):
        pass


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
        self.subimage=self.image[minrow:maxrow+1, mincol:maxcol+1]
        self.subcen=[cen[0]-minrow, cen[1]-mincol]

        self.row_range=[minrow,maxrow]
        self.col_range=[mincol,maxcol]

