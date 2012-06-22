from sys import stderr
import numpy
from numpy import array, logspace, linspace, zeros, ones, \
        sqrt, \
        polyfit, poly1d, log10, arange, where, sin, cos
import fimage
import esutil as eu
import images
import os

sigfac=7.
n=20
smin=0.5
smax=20.
sigma2vals=linspace(smin,smax,n)

def test_ellip():
    re=3.0
    T = 2*re**2
    re2 = re**2

    #theta = random()*360.
    theta=37.
    nsig=15.

    num_e=20
    evals=linspace(0.0,0.6,num_e)

    dt=[('theta','f8'),
        ('e1','f8'),
        ('e2','f8'),
        ('e','f8'),
        ('e1meas','f8'),
        ('e2meas','f8'),
        ('emeas','f8')]
    data=zeros(num_e,dtype=dt)
    for i,e in enumerate(evals):
        print e
        e1 = e*cos(2*theta*numpy.pi/180.0)
        e2 = e*sin(2*theta*numpy.pi/180.0)

        Irr = (1-e1)*T/2.0
        Irc = e2*T/2.0
        Icc = (1+e1)*T/2.0

        cov = [Irr,Irc,Icc]

        dim = int(2*nsig*re)
        if (dim % 2) == 0:
            dim += 1
        dims=array([dim,dim])
        cen=(dims-1)/2.


        im=fimage.model_image('dev',
                              dims,
                              cen,
                              cov,nsub=16)
        if False:
            #images.multiview(im)
            w1,=where(im[cen[0],:] > 0)
            w2,=where(im[:,cen[1]] > 0)
            plt=eu.plotting.bscatter(w1,im[cen[0],w1],type='solid',
                                     color='blue',
                                     ylog=True,show=False)
            plt=eu.plotting.bscatter(w1,im[cen[0],w1],type='filled circle',
                                     color='blue',
                                     plt=plt,show=False)
            plt=eu.plotting.bscatter(w2,im[w2,cen[1]],type='solid',
                                 color='red',plt=plt,show=False)
            eu.plotting.bscatter(w2,im[w2,cen[1]],type='filled circle',
                                 color='red',plt=plt)
            #stop
        stat=fimage.statistics.fmom(im)

        data['e'][i] = e
        data['e1'][i] = e1
        data['e2'][i] = e2
        data['e1meas'][i] = stat['e1']
        data['e2meas'][i] = stat['e2']
        data['emeas'][i] = sqrt(stat['e1']**2 + stat['e2']**2)

    plt=eu.plotting.bscatter(data['e'],
                             data['emeas']-data['e'],
                             xlabel='e requested',
                             ylabel=r'$e_{meas}-e_{req}$')

def test_flux_ap():
    import biggles
    import pcolors 
    nsub=16
    # remember, half light radii are sqrt(10.83) smaller than
    # unweighted sigma, approximately 3.3
    re = 3.0
    sigma=re*sqrt(10.83)
    e=0.3
    theta=27.

    sigma2 = sigma**2
    cov=fimage.conversions.ellip2mom(2*sigma2, e=e, theta=theta)

    dims = array( [400,400])
    print >>stderr,'dims:',dims
    cen = (dims-1)/2.
    im=fimage.model_image('dev',
                          dims,
                          cen,
                          cov,
                          nsub=nsub)

    row,col=numpy.ogrid[0:im.shape[0], 0:im.shape[1]]
    rm = array(row - cen[0], dtype='f8')
    cm = array(col - cen[1], dtype='f8')
    radm = sqrt(rm**2 + cm**2)

    radmax=im.shape[0]/2
    #colors=pcolors.rainbow(radmax-1,'hex')
    radii = numpy.arange(1,radmax+1)
    cnts=numpy.zeros(radii.size)
    for ir,radius in enumerate(radii):
        w=where(radm <= radius)
        if w[0].size > 0:
            print w[0].size,w[1].size,im[w].size
            #color=colors[ir]
            cnts[ir] = im[w].sum()

    r50 = eu.stat.interplin(radii/re, cnts/cnts.max(), 0.5)
    r85 = eu.stat.interplin(radii/re, cnts/cnts.max(), 0.85)
    r90 = eu.stat.interplin(radii/re, cnts/cnts.max(), 0.9)
    print >>stderr,"r50:",r50
    print >>stderr,"r85:",r85
    print >>stderr,"r90:",r90
    plt=eu.plotting.bscatter(radii/re, cnts/cnts.max(),
                             xlabel=r'$r/r_e$',
                             ylabel='counts/max',
                             title='dim: %d' % dims[0],
                             xlog=True,ylog=True,show=False)
    #plt.add(biggles.Point(r50,0.5,color='red'))
    #plt.add(biggles.Point(r90,0.9,color='blue'))

    fname='~/tmp/dev-fluxfrac.rec'
    dt=[('r_over_re','f8'),
        ('fluxfrac','f8')]
    data=zeros(radii.size, dtype=dt)
    data['r_over_re'] = radii/re
    data['fluxfrac'] = cnts/cnts.max()
    eu.io.write(fname,data,clobber=True,verbose=True)
    plt.show()


def test_ap(redo=False):
    import biggles
    import pcolors 
    nsub=16
    # remember, half light radii are sqrt(10.83) smaller than
    # unweighted sigma, approximately 3.3
    re = 1.0
    sigma=re*sqrt(10.83)
    e=0.3
    theta=27.

    sigma2 = sigma**2
    cov=fimage.conversions.ellip2mom(2*sigma2, e=e, theta=theta)

    #cov = array([sigma2,.0,sigma2])

    recf='~/tmp/test-devprof-ap-nsub%02d-e%.2f-theta%.2f.rec' % (nsub,e,theta)
    epsfile=recf.replace('.rec','.eps')
    if not os.path.exists(recf) or redo:
        #dimvals = arange(20,100)
        #dimvals = arange(100,200)
        #dimvals = arange(200,400,8)
        #dimvals = arange(200,400,8)
        dimvals = linspace(50.0,400.0,20).astype('i8')

        #dimvals=[5000]
        #rat = zeros(len(dimvals))
        dt=[('image_radius','f8'),
            ('cov','f8',3),
            ('cov_meas','f8',3)]
        data=zeros(len(dimvals), dtype=dt)
        for i,dim in enumerate(dimvals):
            dims = array( [dim,dim])
            print >>stderr,'dims:',dims
            cen = (dims-1)/2.
            im=fimage.model_image('dev',
                                  dims,
                                  cen,
                                  cov,
                                  nsub=nsub)

            stat=fimage.statistics.fmom(im)

            #print stat
            print dim/sigma/2,dim,stat['cov'][0], stat['cov'][0]/sigma**2

            #rat[i] = stat['cov'][0]/sigma**2
            data['image_radius'][i] = dim/2
            data['cov'][i,:] = cov # constant
            data['cov_meas'][i,:] = stat['cov']

            eu.io.write(recf, data, clobber=True)
    else:
        data=eu.io.read(recf)

    sigma2 = (data['cov'][:,0] + data['cov'][:,2])/2
    sigma2meas = (data['cov_meas'][:,0] + data['cov_meas'][:,2])/2
    sigma2_rat = sigma2meas/sigma2

    Ttrue=data['cov'][:,2]+data['cov'][:,0]
    e1true = (data['cov'][:,2]-data['cov'][:,0])/Ttrue
    e2true = 2*data['cov'][:,1]/Ttrue
    etrue=sqrt(e1true**2 + e2true**2)

    Tmeas=data['cov'][:,2]+data['cov_meas'][:,0]
    e1m = (data['cov_meas'][:,2]-data['cov_meas'][:,0])/Tmeas
    e2m = 2*data['cov_meas'][:,1]/Tmeas
    em=sqrt(e1m**2 + e2m**2)

    radrat = data['image_radius']/sqrt(sigma2)
    tab=biggles.Table(2,1)
    splt=eu.plotting.bscatter(radrat,
                              sigma2_rat, 
                              yrange=[0,1.1],
                              xlabel=r'$image radius/\sigma$',
                              ylabel=r'$\sigma^2_{meas}/\sigma^2_{true}',
                              show=False)
    splt.add(biggles.Curve(radrat,ones(radrat.size),color='blue'))
    eplt=eu.plotting.bscatter(radrat,
                              em-etrue, 
                              xlabel=r'$image radius/\sigma$',
                              ylabel=r'$e_{meas}-e_{true}$',
                              show=False)
    eplt.add(biggles.Curve(radrat,zeros(radrat.size),color='blue'))
    tab[0,0] = splt
    tab[1,0] = eplt
    tab.title = r'$N_{sub}: %d \epsilon:%.2f \theta:%.2f' % (nsub,e,theta)
    tab.show()
    tab.write_eps(epsfile)

def test():

    fname='~/tmp/test-devprof.rec'
    if not os.path.exists(fname):

        dt=[('sigma2','f8'),('sigma2meas','f8')]
        data=zeros(n,dtype=dt)

        for i,sigma2 in enumerate(sigma2vals):
            # sigma2 would be equivalent of a half life radius,
            # and we expect <x^2> propto r_e^2
            sz=int(2*sigma2*sigfac)
            dims=array([sz,sz])
            if (dims[0] % 2) == 0:
                dims += 1
            cen = (dims-1)/2
            print sigma2,dims,cen
            im=fimage.model_image('dev',
                                  dims,
                                  cen,
                                  [sigma2,0.,sigma2],nsub=16)
            stat=fimage.statistics.fmom(im)

            data['sigma2'][i] = sigma2
            data['sigma2meas'][i] = stat['cov'][0]
            print data['sigma2meas'][i]/data['sigma2'][i]

            #images.multiview(im,title='%d' % i)

        eu.io.write(fname,data,verbose=True)
    else:
        data=eu.io.read(fname,verbose=True)


    s=data['sigma2']
    m=data['sigma2meas']

    deg=3
    plyf=polyfit(m, s, deg)
    ply=poly1d(plyf)
    print 'poly:'
    print ply

    plt=eu.plotting.bscatter(m, s, show=False,xlabel=r'$\sigma^2_{meas}$',ylabel=r'$\sigma^2_{req}$')
    eu.plotting.bscatter(m, ply(m),plt=plt,color='red',type='solid')

    return
    #ls = log10(data['sigma2'])
    #lm = log10(data['sigma2meas'])
    deg=1
    plyf=polyfit(lm, ls, deg)
    ply=poly1d(plyf)
    print 'poly:'
    print ply

    #plt=eu.plotting.bscatter(lm, ls, show=False,xlabel=r'$\sigma^2_{meas}$',ylabel=r'$\sigma^2_{req}$')
    #eu.plotting.bscatter(lm, ply(lm),plt=plt,color='red',type='solid')

def dopoly():
    c=[1.658e-06, -0.0005773, 0.1677, 0.6726]
    ply=poly1d(c)

    rat=zeros(n)
    for i,sigma2 in enumerate(sigma2vals):


        sigma2use = ply(sigma2)
        sz=2*sigma2*sigfac
        dims=array([sz,sz])
        if (dims[0] % 2) == 0:
            dims += 1
        cen = (dims-1)/2
        im=fimage.model_image('dev',
                              dims,
                              cen,
                              [sigma2use,0.,sigma2use],nsub=16)
        stat=fimage.statistics.fmom(im)

        rat[i] = stat['cov'][0]/sigma2

        print 'rat:',rat[i]


def dopow():

    #A = .7151
    #B = .3081
    #0.5026 x - 0.5806
    #A = 0.5026
    #B = 0.5806
    A = 0.6982 
    B = - 0.3023
    rat=zeros(n)
    for i,sigma2 in enumerate(sigma2vals):


        sigma2use = 10.**B * sigma2**A
        sz=2*sigma2*sigfac
        dims=array([sz,sz])
        if (dims[0] % 2) == 0:
            dims += 1
        cen = (dims-1)/2
        im=fimage.model_image('dev',
                              dims,
                              cen,
                              [sigma2use,0.,sigma2use],nsub=16)
        stat=fimage.statistics.fmom(im)

        rat[i] = stat['cov'][0]/sigma2

        print 'rat:',rat[i]

        #images.multiview(im)
        #stop


    #eu.plotting.bscatter(sigma2vals, rat, title='rat')

#test()
#dopow()
#dopoly()
test_flux_ap()
#test_ellip()
