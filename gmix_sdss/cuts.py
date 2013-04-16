PSF_EMIN = -0.10
PSF_EMAX =  0.10
S2N_MIN = 20.0
S2N_MAX = 300.0
SRATIO_MIN = 1.0

S2N_NBIN = 40

import numpy

class Selector(object):
    def __init__(self, cols, do_sratio_cut=True):
        """
        send a Columns instance
        """

        self.cols=cols
        self.sratio_cut=do_sratio_cut


    def do_select(self, camcol=None):
        cols=self.cols

        print 'loading select columns'
        print '    sratio'
        sratio = cols['sratio'][:]
        print '    s2n'
        s2n    = cols['s2n'][:]

        logic=(  (s2n > S2N_MIN)
               & (s2n < S2N_MAX)
               & (sratio > SRATIO_MIN) )

        if camcol is not None:
            print '    selecting camcol',camcol
            camcols = cols['camcol'][:]
            logic = logic & (camcols==camcol)

        print 'first cut'
        w,=numpy.where(logic)
        s2n=s2n[w]
        sratio=sratio[w]

        print '    psf pars (subset only)'
        ppars_rec=cols['psf_pars'][w]
        psf_e1=ppars_rec[:,2]
        psf_e2=ppars_rec[:,3]

        print 'second cut'
        w2,=numpy.where(  (psf_e1 > PSF_EMIN)
                        & (psf_e1 < PSF_EMAX)
                        & (psf_e2 > PSF_EMIN)
                        & (psf_e2 < PSF_EMAX) )

        psf_e1 = psf_e1[w2]
        psf_e2 = psf_e2[w2]
        s2n=s2n[w2]
        sratio=sratio[w2]
        w=w[w2]

        # keep these since we already read them
        data={}

        data['psf_e1'] = psf_e1
        data['psf_e2'] = psf_e2
        data['sratio'] = sratio
        data['s2n'] = s2n

        self.data=data
        self.indices=w


