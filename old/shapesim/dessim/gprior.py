"""
TODO:
    - decide if we want to allow the discontinuity that occurs at mag=23, or
    just use the last bin from SDSS

    - future:
        - implement bulge-disk by putting in two entries?
        - allow generic sersic?
"""
import numpy

class GPriorVsMag(object):
    def __init__(self):
        self._load_gpriors()


    def sample2d(self, mags, models):
        """
        models ignored for now
        """
        umodels=numpy.unique(models)
        g1=numpy.zeros(mags.size)
        g2=numpy.zeros(mags.size)

        for model in umodels:
            w,=numpy.where(models==model)
            tg1,tg2 = self._sample2d_by_model(mags, model)

            g1[w] = tg1
            g2[w] = tg2

        return g1,g2

    def _sample2d_by_model(self, mags, model):
        if model=='gauss':
            gpriors=self.gpriors['exp']
        elif model in ['exp','gexp']:
            gpriors=self.gpriors[model]
        elif model in ['dev','gdev']:
            gpriors=self.gpriors[model]
        else:
            raise ValueError("bad model: '%s'" % model)

        ndone=0
        g1=numpy.zeros(mags.size)
        g2=numpy.zeros(mags.size)

        pdict=gpriors[0]
        w,=numpy.where(mags < pdict['minmag'])
        if w.size > 0:
            gprior=pdict['gprior']
            print '    %s %4d at mag < %s' % (model,w.size,pdict['minmag'])
            g1[w],g2[w]= gprior.sample2d(w.size)
            ndone += w.size


        pdict=gpriors[-1]
        w,=numpy.where(mags > pdict['maxmag'])
        if w.size > 0:
            gprior=pdict['gprior']
            print '    %s %4d at mag > %s' % (model,w.size,pdict['maxmag'])
            g1[w],g2[w]= gprior.sample2d(w.size)
            ndone += w.size


        for pdict in gpriors:
            w,=numpy.where(  (mags >= pdict['minmag'])
                           & (mags <= pdict['maxmag']) )
            if w.size > 0:
                gprior=pdict['gprior']
                print '    %s %4d at %5.2f < mag < %5.2f' % (model,w.size,pdict['minmag'],pdict['maxmag'])
                g1[w],g2[w]= gprior.sample2d(w.size)
                ndone += w.size

        if ndone < mags.size:
            raise ValueError("generate only %d/%d" % (ndone,mags.size))

        return g1,g2

    def _load_gpriors(self):
        import cluster_step

        exp_prior_pars=cluster_step.files.read_gprior(type='gexp')
        dev_prior_pars=cluster_step.files.read_gprior(type='gdev')

        gpriors={}

        exp_gpriors=[]
        for i in xrange(exp_prior_pars.size):
            pdict={}
            p=cluster_step.prior.GPriorExp(exp_prior_pars['pars'][i])
            pdict['gprior'] = p
            pdict['minmag'] = exp_prior_pars['minmag'][i]
            pdict['maxmag'] = exp_prior_pars['maxmag'][i]

            exp_gpriors.append(pdict)

        dev_gpriors=[]
        for i in xrange(dev_prior_pars.size):
            pdict={}
            p=cluster_step.prior.GPriorExp(dev_prior_pars['pars'][i])
            pdict['gprior'] = p
            pdict['minmag'] = dev_prior_pars['minmag'][i]
            pdict['maxmag'] = dev_prior_pars['maxmag'][i]

            dev_gpriors.append(pdict)

        gpriors['exp'] = exp_gpriors
        gpriors['dev'] = dev_gpriors

        self.gpriors=gpriors


