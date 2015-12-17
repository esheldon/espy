"""
The data files were actually generated by the IDL code.
"""
from __future__ import print_function
import os
import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1
import numpy
from numpy import pi as PI, cos, sin, unique

import lensing

try:
    import sdsspy
except:
    pass

def read_rotfile(run, type='eq'):
    f = rotfile(run, type)
    print("Reading rotation:",f)
    # read subset to avoid mwrfits bug writing fieldid
    columns=['camcol','field','cos2angle','sin2angle','angle']
    return eu.io.read(f, columns=columns, lower=True)
def rotfile(run, type='eq'):
    d = rotdir(type)
    f = '%srot-%06i-301.fits' % (type, run)
    f = path_join(d, f)
    return f
def rotdir(type='eq'):
    d=os.environ['LENSDIR']
    d = path_join(d, 'sdss-shape-rot', type)
    return d

class SDSSRotator:
    """
    Rotate objecs from CCD coordinates into either equatorial or survey coords. 
    """
    def __init__(self, system='eq'):
        self.system=system
        self.current_run=None

    def rotate_multi(self, runs, camcols, fields, filter, e1pix, e2pix, 
                     getrot=False, verbose=False):

        e1=numpy.zeros(runs.size, dtype='f8')
        e2=e1.copy()
        if getrot:
            angles = e1.copy()

        urun = unique(runs)
        for run in urun:
            if verbose:
                print("  run:",run)
            wrun=where1(runs == run)
            ucamcols = unique(camcols[wrun])
            for camcol in ucamcols:
                if verbose:
                    print("    camcol:",camcol)
                wc=where1( camcols[wrun] == camcol)
                wc=wrun[wc]
                h,rev=eu.stat.histogram(fields[wc], binsize=1, rev=True)
                for i in xrange(h.size):
                    if rev[i] != rev[i+1]:

                        wf = rev[rev[i]:rev[i+1]]
                        wf = wc[wf]

                        #uf = unique(fields[wf])
                        #if uf.size != 1:
                        #    raise ValueError("expected unique field")

                        field = fields[wf[0]]
                        #print("      field:",field)

                        res=self.rotate(run, camcol, field, filter, e1pix[wf], e2pix[wf], 
                                        getrot=getrot, verbose=verbose)
                        e1[wf] = res[0]
                        e2[wf] = res[1]
                        if getrot:
                            angles[wf] = res[2]

        if getrot:
            return e1,e2,angles
        else:
            return e1,e2




    def rotate(self, run, camcol, field, filter, e1pix, e2pix,
               getrot=False, verbose=False):
        if not numpy.isscalar(run):
            return self.rotate_multi(run,camcol,field,filter,e1pix,e2pix,
                                     getrot=getrot, verbose=verbose)

        self.load_run(run)

        if verbose:
            print("%06i-%i-%04i" % (run,camcol,field))

        w=where1((self.rotstruct['camcol'] == camcol) & (self.rotstruct['field'] == field) )
        if w.size != 1:
            print("expected single match for %06i-%i-%04i, got %i" % (run,camcol,field,w.size))
            fields = numpy.unique(self.rotstruct['field'])
            print("here are the fields:",fields)
            raise ValueError("stopping")

        filternum = sdsspy.FILTERNUM[filter]
        cos2angle = self.rotstruct['cos2angle'][w,filternum]
        sin2angle = self.rotstruct['sin2angle'][w,filternum]

        # hmm... this seems to be applying -2*angle rotation...
        e1 =  e1pix*cos2angle + e2pix*sin2angle
        e2 = -e1pix*sin2angle + e2pix*cos2angle

        if getrot:
            angle = self.rotstruct['angle'][w,filternum]
            angles = numpy.zeros(e1.size, dtype='f8') + angle[0]
            return e1, e2, angles
        else:
            return e1,e2

    def load_run(self, run):
        if self.current_run != run:
            self.rotstruct = read_rotfile(run, self.system)
            self.current_run = run

def make_rotation_test_data(run=4335, camcol=3, field=100):
    c=lensing.regauss.open_columns('04')
    # note logic returns indices for Columns columns!
    w=c['run'] == run
    if w.size == 0:
        raise ValueError("no objects found for run %d" % run)
    
    flags = c['corrflags_rg_r'][w]
    camcols = c['camcol'][w]
    fields = c['field'][w]
    w2=where1((camcols == camcol) & (flags == 0) & (fields==field))
    #w2=where1( camcols == camcol )
    if w2.size == 0:
        raise ValueError("not objects in camcol %s with flags==0" % camcol)
    w=w[w2]

    data = c.read_columns(['ra','dec','run','camcol','field','e1_rg_r','e2_rg_r'], rows=w)

    output = numpy.zeros(data.size, 
                          dtype=[('ra','f8'),('dec','f8'),
                                 ('g1eq','f8'),('g2eq','f8'),
                                 ('clambda','f8'),('ceta','f8'),
                                 ('g1survey','f8'),('g2survey','f8')])

    output['ra'] = data['ra']
    output['dec'] = data['dec']

    lam, eta = eu.coords.eq2sdss(data['ra'],data['dec'])
    output['clambda'] =  lam
    output['ceta'] = eta

    eq_rotator = SDSSRotator('eq')
    survey_rotator = SDSSRotator('survey')

    g1eq_alt,g2eq_alt = eq_rotator.rotate(data['run'],data['camcol'],data['field'],2,
                                          data['e1_rg_r'],data['e2_rg_r'],verbose=True)
    g1eq_alt /= 2 
    g2eq_alt /= 2 


    output_file=os.path.expanduser('~/tmp/test-rot/test-rot.rec')
    print("writing output file:",output_file)
    fobj = open(output_file,'w')
    num = numpy.array([output.size], dtype='i8')
    num.tofile(fobj)
    #robj = eu.recfile.Recfile(output_file, 'w', delim=' ')
    robj = eu.recfile.Recfile(fobj, 'r+')
    robj.write(output)
    robj.close()



class ComparePrinceton(dict):

    def compare(self, run):
        import biggles
        import pcolors

        pdata, mdata = self.load_data(run)
        if len(pdata) == 0:
            print("no princeton data found")
            return
        if len(mdata) == 0:
            print("no my data found")
            return

        tab = biggles.Table(2,2)

        pcos_plt = biggles.FramedPlot()
        psin_plt = biggles.FramedPlot()
        mcos_plt = biggles.FramedPlot()
        msin_plt = biggles.FramedPlot()

        pcos_plots=[]
        psin_plots=[]
        mcos_plots=[]
        msin_plots=[]

        colors=pcolors.rainbow(6, 'hex')

        for camcol in xrange(1,6+1):
            # first princeton
            wp = where1(pdata['camcol'] == camcol)

            bcos = eu.stat.Binner(pdata['field'][wp], cos(2*pdata['phi_offset'][wp]))
            bcos.dohist(binsize=1.0)
            bcos.calc_stats()


            wgood_cos = where1(bcos['hist'] > 0)
            pcos = biggles.Curve(bcos['xmean'][wgood_cos],bcos['ymean'][wgood_cos], color=colors[camcol-1])
            pcos.label = 'princ camcol %s' % camcol

            pcos_plt.add(pcos)
            pcos_plots.append(pcos)


            bsin = eu.stat.Binner(pdata['field'][wp], sin(2*pdata['phi_offset'][wp]))
            bsin.dohist(binsize=1.0)
            bsin.calc_stats()

            wgood_sin = where1(bsin['hist'] > 0)
            psin = biggles.Curve(bsin['xmean'][wgood_sin],bsin['ymean'][wgood_sin], color=colors[camcol-1])
            psin.label = 'princ camcol %s' % camcol

            psin_plt.add(psin)
            psin_plots.append(psin)



            # now mine
            wm = where1(mdata['camcol'] == camcol)
            mpcos = biggles.Curve(mdata['field'][wm], cos(2*mdata['angle'][wm,2]), color=colors[camcol-1])
            mpcos.label = 'mine camcol %s' % camcol
            mcos_plt.add(mpcos)
            mcos_plots.append(mpcos)

            wm = where1(mdata['camcol'] == camcol)
            mpsin = biggles.Curve(mdata['field'][wm], sin(2*mdata['angle'][wm,2]), color=colors[camcol-1])
            mpsin.label = 'mine camcol %s' % camcol
            msin_plt.add(mpsin)
            msin_plots.append(mpsin)



        # princeton stuff
        pcos_key = biggles.PlotKey(0.1,0.9,pcos_plots)
        pcos_plt.add(pcos_key)

        pcos_plt.xlabel = 'Field'
        pcos_plt.title = 'Run: %s' % run
        pcos_plt.ylabel = 'cos(2*angle)'

        psin_key = biggles.PlotKey(0.1,0.9,psin_plots)
        psin_plt.add(psin_key)

        psin_plt.xlabel = 'Field'
        psin_plt.title = 'Run: %s' % run
        psin_plt.ylabel = 'sin(2*angle)'

        tab[0,0] = pcos_plt
        tab[0,1] = psin_plt

        # my stuff
        mcos_key = biggles.PlotKey(0.1,0.9,mcos_plots)
        mcos_plt.add(mcos_key)

        mcos_plt.xlabel = 'Field'
        mcos_plt.title = 'Run: %s' % run
        mcos_plt.ylabel = 'cos(2*angle)'

        msin_key = biggles.PlotKey(0.1,0.9,msin_plots)
        msin_plt.add(msin_key)

        msin_plt.xlabel = 'Field'
        msin_plt.title = 'Run: %s' % run
        msin_plt.ylabel = 'sin(2*angle)'

        tab[1,0] = mcos_plt
        tab[1,1] = msin_plt


        tab.show()
            

    def load_data(self, run):
        import pgnumpy
        print("loading princeton data")
        query="""
        select 
            run,
            camcol,
            field,
            phi_offset
        from 
            scat_princeton
        where 
            run = {run}
        """.format(run=run)
        print(query)

        pg=pgnumpy.PgNumpy()
        princeton_data = pg.fetchall(query)


        print("loading my rotation struct")
        mydata = read_rotfile(run)

        return princeton_data, mydata
