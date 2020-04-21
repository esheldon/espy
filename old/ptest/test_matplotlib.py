import matplotlib
import pylab
import numpy

def t_mplib():

    n = 20
    x = numpy.arange(n)

    y = x*x


    # control font size on axis titles
    pylab.rc('axes',labelsize=18)

    pylab.clf()
    #pylab.plot(x, y, marker='o',linestyle='.')
    pylab.plot(x, y, marker='o',linestyle='None')

    
    pylab.xlabel(r'$\alpha$')
    pylab.ylabel(r'$\Delta\Sigma~~[$M$_{\odot} $pc$^{-2}]$')

    #pylab.show()
