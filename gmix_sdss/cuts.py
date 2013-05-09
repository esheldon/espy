import numpy

# just the defaults, broader than we will use
PSF_EMIN = -0.10
PSF_EMAX =  0.10
S2N_MIN = 10.0
S2N_MAX = 300.0
SRATIO_MIN = 0.5

S2N_NBIN = 40


class Selector(object):
    def __init__(self, cols, psf_model, do_sratio_cut=True):
        """
        send a Columns instance
        """

        self.cols=cols
        self.psf_model=psf_model
        self.sratio_cut=do_sratio_cut

    def do_select(self,
                  camcol=None,
                  pzrun=None,
                  s2n_min=S2N_MIN,
                  s2n_max=S2N_MAX,
                  psf_emin=PSF_EMIN,
                  psf_emax=PSF_EMAX,
                  sratio_min=SRATIO_MIN):
        """
        Send pzrun to demand matches to the indicated
        photoz run
        """
        cols=self.cols

        ntot=cols['s2n'].size

        print 'loading select columns'
        print '    sratio'
        sratio = cols['sratio'][:]
        print '    s2n'
        s2n    = cols['s2n'][:]

        logic=(s2n > s2n_min)
        wt,=numpy.where(logic)
        print("s2n > %.2f: %d/%d" % (s2n_min, wt.size,ntot))

        logic=logic & (s2n < s2n_max)
        wt,=numpy.where(logic)
        print("s2n < %.2f: %d/%d" % (s2n_max, wt.size,ntot))

        logic=logic & (sratio > sratio_min)
        wt,=numpy.where(logic)
        print("sratio > %.2f: %d/%d" % (sratio_min, wt.size,ntot))


        if camcol is not None:
            print '    selecting camcol',camcol
            camcols = cols['camcol'][:]
            logic = logic & (camcols==camcol)

            wt,=numpy.where(logic)
            print("column cut: %d/%d" % (wt.size,ntot))

        if pzrun is not None:
            print '    selecting zphot matches'
            colname='match_zphot%s' % pzrun
            m=cols[colname][:]
            logic = logic & (m >= 0)

            wt,=numpy.where(logic)
            print("photoz match cut: %d/%d" % (wt.size,ntot))

        w,=numpy.where(logic)

        s2n=s2n[w]
        sratio=sratio[w]

        print '    psf pars (subset only)'
        ppars_rec=cols['psf_pars'][w]
        if 'em' in self.psf_model:
            # not quite right since they may be non-cocentric
            psf_e1,psf_e2=_extract_psf_e1e2_full(ppars_rec)
        else:
            psf_e1=ppars_rec[:,2]
            psf_e2=ppars_rec[:,3]

        print 'epsf cut [%.2f,%.2f]' % (psf_emin,psf_emax)
        w2,=numpy.where(  (psf_e1 > psf_emin)
                        & (psf_e1 < psf_emax)
                        & (psf_e2 > psf_emin)
                        & (psf_e2 < psf_emax) )

        print("epsf cut: %d/%d" % (w2.size,ntot))

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


def _extract_psf_e1e2_full(pars):
    """
    not exactly right for non-cocentric
    """
    nobj=pars.shape[0]
    ngauss=pars.shape[1]/6

    irrsum=numpy.zeros(nobj)
    ircsum=numpy.zeros(nobj)
    iccsum=numpy.zeros(nobj)
    #psum=numpy.zeros(nobj)

    # sum irr*p, irc*p, icc*p
    # e1=(icc-irr)/(icc+irr)
    for i in xrange(ngauss):
        start=i*6
        end=(i+1)*6

        pt=pars[:,start]
        irrt=pars[:,start+3]
        irct=pars[:,start+4]
        icct=pars[:,start+5]

        #psum += pt
        irrsum += pt*irrt
        ircsum += pt*irct
        iccsum += pt*icct

    Tsum=irrsum+iccsum

    e1=(iccsum-irrsum)/Tsum
    e2=2*ircsum/Tsum

    return e1,e2
