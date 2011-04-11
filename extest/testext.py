import numpy
import cosmolib
import esutil
import time

def main():
    c=cosmolib.cosmo(2997.9245799999999, 1, 0.3, 0.7, 0.0)
    cold = esutil.cosmology.Cosmo()

    zmin =numpy.array([0.1,0.2], dtype='f8')
    zmax =numpy.array([0.3,0.4], dtype='f8')

    print 'Dc scalar'
    print c.cdist(0.2, 0.3)
    print cold.Dc(0.2, 0.3)

    print 'Dc vec1'
    print c.cdist_vec1(zmin, 0.3)
    print cold.Dc(zmin, 0.3)

    print 'Dc vec2'
    print c.cdist_vec2(0.1, zmax)
    print cold.Dc(0.1, zmax)

    print 'Dc 2vec'
    print c.cdist_2vec(zmin, zmax)
    print cold.Dc(zmin, zmax)


    print '\nDm scalar'
    print c.tcdist(0.2, 0.3)
    print cold.Dm(0.2, 0.3)

    print 'Dm vec1'
    print c.tcdist_vec1(zmin, 0.3)
    print cold.Dm(zmin, 0.3)

    print 'Dm vec2'
    print c.tcdist_vec2(0.1, zmax)
    print cold.Dm(0.1, zmax)

    print 'Dm 2vec'
    print c.tcdist_2vec(zmin, zmax)
    print cold.Dm(zmin, zmax)


    print '\nDa scalar'
    print c.angdist(0.2, 0.3)
    print cold.Da(0.2, 0.3)

    print 'Da vec1'
    print c.angdist_vec1(zmin, 0.3)
    print cold.Da(zmin, 0.3)

    print 'Da vec2'
    print c.angdist_vec2(0.1, zmax)
    print cold.Da(0.1, zmax)

    print 'Da 2vec'
    print c.angdist_2vec(zmin, zmax)
    print cold.Da(zmin, zmax)

    print '\nDl scalar'
    print c.lumdist(0.2, 0.3)
    print cold.Dl(0.2, 0.3)

    print 'Dl vec1'
    print c.lumdist_vec1(zmin, 0.3)
    print cold.Dl(zmin, 0.3)

    print 'Dl vec2'
    print c.lumdist_vec2(0.1, zmax)
    print cold.Dl(0.1, zmax)

    print 'Dl 2vec'
    print c.lumdist_2vec(zmin, zmax)
    print cold.Dl(zmin, zmax)


    print '\ndV scalar'
    print c.dV(0.2)
    print cold.dV(0.2)

    print 'dV vec'
    print c.dV_vec(zmin)
    print cold.dV(zmin)

    print '\nvolume'
    print c.V(0.2, 0.3)
    print cold.V(0.2, 0.3)

    print '\nsigmacritinv scalar'
    print c.scinv(0.2, 0.3)
    print cold.sigmacritinv(0.2, 0.3)

    print 'sigmacritinv vec1'
    print c.scinv_vec1(zmin, 0.3)
    print cold.sigmacritinv(zmin, 0.3)

    print 'sigmacritinv vec2'
    print c.scinv_vec2(0.1, zmax)
    print cold.sigmacritinv(0.1, zmax)

    print 'sigmacritinv 2vec'
    print c.scinv_2vec(zmin, zmax)
    print cold.sigmacritinv(zmin, zmax)



    return


    ntime=10
    zmin = numpy.linspace(0.1, 0.8, 1000000)

    print '\nDc timings'

    tm=time.time()
    for i in xrange(ntime):
        dold = cold.Dc(zmin, 0.3)
    print 'time old:',time.time()-tm

    tm=time.time()
    for i in xrange(ntime):
        d = c.cdist_vec1(zmin, 0.3)
    print 'time:',time.time()-tm

    print '\nDm timings'

    tm=time.time()
    for i in xrange(ntime):
        dold = cold.Dm(zmin, 0.3)
    print 'time old:',time.time()-tm

    tm=time.time()
    for i in xrange(ntime):
        d = c.tcdist_vec1(zmin, 0.3)
    print 'time:',time.time()-tm


main()
