import os
import deswl
import numpy
from numpy import tanh, arctanh, sqrt
import fimage
import admom
import biggles

def run_trials(nsub, rg_admom_nsub):
    n=100
    rg_image_nsub = rg_admom_nsub 
    psf_order=10
    gal_shear_order = 8

    gal_theta = numpy.linspace(0.01843243,179.13423,n)
    spacing = gal_theta[1]-gal_theta[0]
    gal_theta += numpy.random.random(n)*spacing
    gal_T = 3.1+2.7
    gal_T = gal_T*10
    gal_e = 0.3

    # reduced shear required to make round
    shear = tanh( 0.5*arctanh(gal_e) )

    dt=[('gal_theta','f8'),
        ('shear1','f8'),
        ('shear2','f8'),
        ('shear','f8'),
        ('shear1meas','f8'),
        ('shear2meas','f8'),
        ('shearmeas','f8')]

    wl_shear=numpy.zeros(n, dtype=dt)
    am_shear=numpy.zeros(n, dtype=dt)

    wl_shear['gal_theta'] = gal_theta
    wl_shear['shear'] = shear
    am_shear['gal_theta'] = gal_theta
    am_shear['shear'] = shear

    for i in xrange(n):
        galcov=fimage.conversions.ellip2mom(gal_T,e=gal_e,
                                            theta=gal_theta[i])

        psfcov=[2.0,0.0,2.0]
        #psfcov=[2.0*5,0.0,2.0*5]
        imcov=[galcov[0]+psfcov[0], 
               galcov[1]+psfcov[1], 
               galcov[2]+psfcov[2]]

        im_sigma = fimage.mom2sigma(imcov[0]+imcov[2])
        dims = int(2*6*im_sigma)
        if (dims % 2) == 0:
            dims+=1
        if dims < 31: dims=31
        dims = [dims]*2

        shear1=0.5*(galcov[2]-galcov[0])/(galcov[2]+galcov[0])
        shear2=0.5*2*galcov[1]/(galcov[2]+galcov[0])
        wl_shear['shear1'][i] = shear1
        wl_shear['shear2'][i] = shear2
        am_shear['shear1'][i] = shear1
        am_shear['shear2'][i] = shear2
    

        wlq = deswl.cwl.WLQuick(psf_order,gal_shear_order)

        sky=0.0

        psfcen=[(dims[0]-1.)/2., (dims[0]-1.)/2.]
        psf = fimage.model_image('gauss',dims,psfcen,psfcov,
                                 nsub=nsub,counts=1)

        imcen=[(dims[0]-1.)/2., (dims[0]-1.)/2.]
        im = fimage.model_image('gauss',dims,imcen,imcov,
                                 nsub=nsub,counts=10000)

        aperture = 4.0*fimage.mom2sigma(imcov[0]+imcov[2])

        skysig=1.0
        wlq.set_psf(psf,psfcen[0],psfcen[1],sky)
        wlq.set_image(im,imcen[0],imcen[1],sky,skysig,aperture)

        guess=2.0
        flags = wlq.calculate_psf_sigma(guess)
        if flags != 0:
            print 'sigma flags:',flags
            return
        #print 'input psf sigma:',fimage.mom2sigma(psfcov[0]+psfcov[2])
        #print 'measured sigma: ',wlq.get_psf_sigma()

        flags += wlq.calculate_psf_shapelets()
        if flags != 0:
            print 'psf shapelets flags:',flags
            return

        #print
        flags += wlq.calculate_shear()
        if flags != 0:
            print 'shear flags:',flags
            return

        wl_shear['shear1meas'][i] = - wlq.get_shear1()
        wl_shear['shear2meas'][i] =   wlq.get_shear2()
        wl_shear['shearmeas'][i] = sqrt(wlq.get_shear1()**2 + 
                                        wlq.get_shear2()**2)

        #wl_shear['shear1meas'][i] = - 0.5*wlq.get_e1()
        #wl_shear['shear2meas'][i] = - 0.5*wlq.get_e2()

        rgkeys={}
        rgkeys['guess'] = (galcov[0] + galcov[2])/2
        rgkeys['guess_psf'] = (psfcov[0] + psfcov[2])/2
        rgkeys['admom_nsub'] = rg_admom_nsub
        rgkeys['image_nsub'] = rg_image_nsub
        #rgkeys['verbose'] = True
        rg = admom.ReGauss(im, imcen[0], imcen[1], psf, **rgkeys)
        rg.do_all()

        #am_shear['shear1meas'][i]=0.5*rg['rgcorrstats']['e1']
        #am_shear['shear2meas'][i]=0.5*rg['rgcorrstats']['e2']
        am_shear['shear1meas'][i]=0.5*rg['corrstats']['e1']
        am_shear['shear2meas'][i]=0.5*rg['corrstats']['e2']

        ame = sqrt(rg['corrstats']['e1']**2 + rg['corrstats']['e2']**2)

        am_shear['shearmeas'][i] = tanh(0.5*arctanh(ame))

        s='theta: %g in: %g %g wl: %g %g %g %g am: %g %g %g %g'

        s = s % (gal_theta[i],
                 shear1,shear2,
                 wl_shear['shear1meas'][i], 
                 wl_shear['shear1meas'][i]/shear1-1,
                 wl_shear['shear2meas'][i], 
                 wl_shear['shear2meas'][i]/shear2-1,
                 am_shear['shear1meas'][i], 
                 am_shear['shear1meas'][i]/shear1-1,
                 am_shear['shear2meas'][i], 
                 am_shear['shear2meas'][i]/shear2-1)
        print s
    return wl_shear, am_shear

def main():
    nsub = 4
    admom_nsub=4

    wl_shear,am_shear=run_trials(nsub, admom_nsub)
    wplt=biggles.FramedPlot()
    aplt=biggles.FramedPlot()

    xlabel=r'$\theta$'
    ylabel = r'$\gamma^{Meas} - \gamma^{True}$'
    wplt.title = 'DESWL Shapelets'
    wplt.xlabel = xlabel
    wplt.ylabel = ylabel

    aplt.title = 'Adaptive Moments'
    aplt.xlabel = xlabel
    aplt.ylabel = ylabel

    wp = biggles.Points(wl_shear['gal_theta'],
                        wl_shear['shearmeas']-wl_shear['shear'],
                        type='filled circle')

    ap = biggles.Points(am_shear['gal_theta'],
                        am_shear['shearmeas']-am_shear['shear'],
                        type='filled circle')

    wplt.add(wp)
    aplt.add(ap)
    wplt.show()
    aplt.show()

    d=os.path.expanduser('~/tmp')
    
    deswl_name = 'deswl-pixel-effect-nsub%02i.eps' % nsub
    am_name = 'deswl-pixel-effect-nsub%02i-am%02i.eps' % (nsub,admom_nsub)
    print deswl_name
    wplt.write_eps(os.path.join(d,deswl_name))
    print am_name
    aplt.write_eps(os.path.join(d,am_name))
main()
