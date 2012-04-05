
"""
    %prog [options] lensrun randrun

Description:

    Correct the delta sigma for 

        - subtract the signal from randoms
        - 1/ssh  the shear polarizability correction
        - clustering of sources with lenses
        
    Note the subtraction of randoms should happen before the other two.  You
    need to pick a radius outside of which the subtraction occurs if you don't
    like the default.

    The inputs to this program are the binned outputs from /bin/bin-lenses.py
    and the matched random outputs from /bin/match-randoms.py"""

from __future__ import print_function
import sys, os
import lensing
import numpy
from numpy import where
from optparse import OptionParser
import biggles
import esutil as eu
import converter

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")

parser.add_option("-s",dest="subtract_rand",action='store_true',default=False,
                  help="Subtract the randoms.  default %default")
parser.add_option("-r",dest="minrad",default=1.,
                  help="The minimum radius for the subtraction of random signal. "
                       "Set to a large number to turn off subtraction. Default %default Mpc")

parser.add_option("--remove",action='store_true',default=False,
                  help="Use results from hist_match_remove.  default %default")
parser.add_option("-p",dest="prompt",action='store_true',default=False,
                  help="Prompt between each plot.  default %default")
parser.add_option("--show",action='store_true',default=False,
                  help="Show plots on screen.  default %default")

def doplot(binned_data, corr_data, rand, label, show=False):
    tab = biggles.Table(1,2)
    arr=lensing.plotting.plot2dsig(binned_data['r'], 
                                   binned_data['dsig'], 
                                   binned_data['dsigerr'],
                                   rand['r'],
                                   rand['dsig'], 
                                   rand['dsigerr'],
                                   plot_label=label, 
                                   label1='data', 
                                   label2='random',
                                   range4var=[0.1,100],
                                   show=False)

    
    tab[0,0] = arr


    cplt = eu.plotting.bscatter(corr_data['r'], 
                                corr_data['clust_corr']-1, 
                                yerr=corr_data['clust_corr_err'],
                                xlabel=lensing.plotting.labels['rproj'], 
                                ylabel='Corr-1', show=False, xlog=True, ylog=True, 
                                size=2)

    cplt_label = biggles.PlotLabel(0.9,0.9,label,halign='right')
    cplt.add(cplt_label)
    cplt.aspect_ratio=1


    ddict = lensing.plotting.plot_dsig(corr_data,show=False)
    ddict['p'].label = 'corrected'
    dplt = ddict['plt']

    # just correct for 1/ssh
    dsig = binned_data['dsig']/corr_data['ssh']

    nocorrp = biggles.Points(binned_data['r'], dsig, color='red')
    nocorrp.label = 'not corrected'
    dplt.add(nocorrp)
    key = biggles.PlotKey(0.9,0.9,[ddict['p'],nocorrp],halign='right')
    dplt.add(key)


    ctab = biggles.Table(2,1)
    ctab[0,0] = cplt
    ctab[1,0] = dplt

    tab[0,1] = ctab

    if show:
        tab.show()

    return tab

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]
    randrun=args[1]

    bintype=options.bintype
    nbin=int(options.nbin)
    minrad=float(options.minrad)
    subtract_rand = bool(options.subtract_rand)
    remove=options.remove

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    if subtract_rand:
        print("Will subtract randoms")

    b = lensing.binning.instantiate_binner(bintype, nbin)

    binned_data = lensing.files.sample_read('binned', lensrun, name=b.name())
    if remove:
        extra='randmatch-rm-%s' % randrun
    else:
        extra='randmatch-%s' % randrun
    allrand = lensing.files.sample_read('binned', lensrun, name=b.name(), extra=extra)

    alldata = lensing.correct.correct(binned_data, allrand, 
                                      subtract_rand=subtract_rand, minrad=minrad)

    lensing.files.sample_write(alldata,'corrected',lensrun,name=b.name(),extra=extra) 



    # now some plots
    biggles.configure('screen','width', 1100)
    biggles.configure('screen','height', 1100)

    
    range4var = [0.1,100]
    for binnum in xrange(nbin):
        eps_corr_extra='correction-%02d' % binnum
        eps_corr=lensing.files.sample_file('corrected-plots',
                                           lensrun,
                                           name=b.name(),
                                           extra=eps_corr_extra, ext='eps')
        eps_rand_extra='randcomp-%02d' % binnum
        eps_rand=lensing.files.sample_file('corrected-plots',
                                           lensrun,
                                           name=b.name(),
                                           extra=eps_rand_extra, ext='eps')
        eps_dsigcorr_extra='dsigcorr-%02d' % binnum
        eps_dsigcorr=lensing.files.sample_file('corrected-plots',
                                               lensrun,
                                               name=b.name(),
                                               extra=eps_dsigcorr_extra, ext='eps')

        eps_all_extra='allcorr-%02d' % binnum
        eps_all=lensing.files.sample_file('corrected-plots',
                                          lensrun,
                                          name=b.name(),
                                          extra=eps_all_extra, ext='eps')



        lensing.files.make_dir_from_path(eps_corr)

        data = alldata[binnum]
        rand = allrand[binnum]

        label=b.bin_label(binnum)

        tab=doplot(binned_data[binnum], data, rand, label, show=options.show)

        # the corr(r)-1 plot
        tab[0,1][0,0].write_eps(eps_corr)
        converter.convert(eps_corr, dpi=90, verbose=True)

        # the corr(r)-1 plot
        tab[0,1][1,0].write_eps(eps_dsigcorr)
        converter.convert(eps_dsigcorr, dpi=90, verbose=True)

        # the rand comparison plot
        tab[0,0].write_eps(eps_rand)
        converter.convert(eps_rand, dpi=120, verbose=True)

        # all
        tab.write_eps(eps_all)
        converter.convert(eps_all, dpi=150, verbose=True)


        if options.prompt:
            key=raw_input("hit a key (q to quit): ")
            if key == 'q':
                return

    d=os.path.dirname(eps_corr)
    os.chdir(d)

    outfile = os.path.join('correction.html')
    pattern = eps_corr.replace('%02d.eps' % (nbin-1,), '*.png')
    pattern=os.path.basename(pattern)
    print("making correction html file:",outfile)
    os.system('im2html -p '+pattern+' > '+outfile)

    outfile = os.path.join('randcomp.html')
    pattern = eps_rand.replace('%02d.eps' % (nbin-1,), '*.png')
    pattern=os.path.basename(pattern)
    print("making rand compare html file:",outfile)
    os.system('im2html -p '+pattern+' > '+outfile)

    outfile = os.path.join('dsigcorr.html')
    pattern = eps_dsigcorr.replace('%02d.eps' % (nbin-1,), '*.png')
    pattern=os.path.basename(pattern)
    print("making disg corr html file:",outfile)
    os.system('im2html -p '+pattern+' > '+outfile)

    outfile = os.path.join('allcorr.html')
    pattern = eps_all.replace('%02d.eps' % (nbin-1,), '*.png')
    pattern=os.path.basename(pattern)
    print("making all corr html file:",outfile)
    os.system('im2html -p '+pattern+' > '+outfile)



main()
