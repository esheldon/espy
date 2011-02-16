from __future__ import print_function
import numpy
from numpy import pi,cos,sin,exp,sqrt
import scipy.integrate
import scipy.signal
import scipy.ndimage

import fimage
import admom
import biggles

import pprint

def test_convolve():
    n=14
    facvals = numpy.linspace(1.0,6.0,n)
    
    sqrt_det_fft=numpy.zeros(n)
    sqrt_det_fconv=numpy.zeros(n)
    rat = numpy.zeros(n)

    for i in xrange(n):

        fac=facvals[i]
        fac2=fac**2

        dims=[41*fac,41*fac]
        cen=[20*fac,20*fac]

        covar = [1.5*fac2,0.0,1.5*fac2]
        img=fimage.model_image('gauss',dims,cen,covar)

        g1covar = [1.4*fac2,0.0,1.4*fac2]
        g2covar = [3.0*fac2,0.0,3.0*fac2]

        b = 0.1
        dg = fimage.double_gauss(dims,cen,b, g1covar, g2covar)
        imc_fft = scipy.signal.fftconvolve(img,dg,mode='same')

        # now do it by convolving the image with each gaussian and
        # adding appropriately

        im1 = fimage.fconv.gaussconv(img, g1covar)
        im2 = fimage.fconv.gaussconv(img, g2covar)

        det1 = g1covar[0]*g1covar[2] - g1covar[1]**2
        det2 = g2covar[0]*g2covar[2] - g2covar[1]**2
        s2 = sqrt(det2/det1)
        imc_fconv = (im1 + b*s2*im2)/(1+b*s2)

        mom_fft=fimage.stat.moments(imc_fft)
        mom_fconv=fimage.stat.moments(imc_fconv)

        print("moments from fft")
        pprint.pprint(mom_fft)
        print("moments from fconv")
        pprint.pprint(mom_fconv)

        cov_fft   = mom_fft['cov']
        cov_fconv = mom_fconv['cov']
        sqrt_det_fft[i] = sqrt(cov_fft[0]*cov_fft[2] - cov_fft[1]**2)
        sqrt_det_fconv[i] = sqrt(cov_fconv[0]*cov_fconv[2] - cov_fconv[1]**2)

    rd = {'sqrt_det_fft':sqrt_det_fft,'sqrt_det_fconv':sqrt_det_fconv,'facvals':facvals}

    plotconv(rd)
    return rd
        
def plotconv(resdict):
    sqrt_det_fft=resdict['sqrt_det_fft']
    sqrt_det_fconv=resdict['sqrt_det_fconv']
    facvals = resdict['facvals']

    rat = sqrt_det_fft/sqrt_det_fconv

    tab=biggles.Table(2,1)


    plt_both = biggles.FramedPlot()
    #pfft = biggles.Curve(facvals, sqrt_det_fft, color='red') 
    pfft = biggles.Curve(facvals, sqrt_det_fft/facvals**2, color='red') 
    pfft.label = 'sqrt(det fft)$'
    #pfconv = biggles.Curve(facvals, sqrt_det_fconv, color='blue') 
    pfconv = biggles.Curve(facvals, sqrt_det_fconv/facvals**2, color='blue') 
    pfconv.label = 'sqrt(det fconv)'

    bkey = biggles.PlotKey(0.1,0.9,[pfft,pfconv],halign='left')
    plt_both.add(pfft,pfconv,bkey)

    plt_both.xlabel = 'fac'
    plt_both.ylabel = r'$sqrt(det)/fac^2$'



    pltrat = biggles.FramedPlot()
    pltrat.add(biggles.Curve(facvals,rat))
    pltrat.xlabel = 'fac'
    pltrat.ylabel = 'sqrt(det_fft/det_fconv)'
    pltrat.yrange=[1,1.03]

    tab[0,0] = plt_both
    tab[1,0] = pltrat

    tab.show()



