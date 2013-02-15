import numpy

class GPriorVsMagNoSplit(object):
    def __init__(self):
        self._load_gpriors()


    def sample2d(self, mags):

        ndone=0
        g1=numpy.zeros(mags.size)
        g2=numpy.zeros(mags.size)

        pdict=self.gpriors[0]
        w,=numpy.where(mags < pdict['minmag'])
        if w.size > 0:
            gprior=pdict['gprior']
            print '    %4d at mag < %s' % (w.size,pdict['minmag'])
            g1[w],g2[w]= gprior.sample2d(w.size)
            ndone += w.size

        pdict=self.gpriors[-1]
        w,=numpy.where(mags > pdict['maxmag'])
        if w.size > 0:
            gprior=pdict['gprior']
            print '    %4d at mag > %s' % (w.size,pdict['maxmag'])
            g1[w],g2[w]= gprior.sample2d(w.size)
            ndone += w.size


        for pdict in self.gpriors:
            w,=numpy.where(  (mags >= pdict['minmag'])
                           & (mags <= pdict['maxmag']) )
            if w.size > 0:
                gprior=pdict['gprior']
                print '    %4d at %5.2f < mag < %5.2f' % (w.size,pdict['minmag'],pdict['maxmag'])
                g1[w],g2[w]= gprior.sample2d(w.size)
                ndone += w.size

        if ndone < mags.size:
            raise ValueError("generate only %d/%d" % (ndone,mags.size))
        return g1,g2

    def _load_gpriors(self):
        import cluster_step
        prior_pars=cluster_step.files.read_gprior(nosplit=True)

        gpriors=[]
        for i in xrange(prior_pars.size):
            pdict={}
            p=cluster_step.prior.GPriorExp(prior_pars['pars'][i])
            pdict['gprior'] = p
            pdict['minmag'] = prior_pars['minmag'][i]
            pdict['maxmag'] = prior_pars['maxmag'][i]

            gpriors.append(pdict)

        self.gpriors=gpriors


