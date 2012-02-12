import os,sys
import recfile
import esutil as eu
from esutil.numpy_util import where1
import numpy
import deswl
import des
import converter

def make_html(pd,merun,serun):
    html="""
<html>
    <body bgcolor=white>

        <p>
        <img src="%(merun)s-%(serun)s-match-shear1diff-vs-shear_s2n.png">
        <img src="%(merun)s-%(serun)s-match-shear2diff-vs-shear_s2n.png">
        <img src="%(merun)s-%(serun)s-match-sizediff-vs-shear_s2n.png">
        <p>
        <img src="%(merun)s-%(serun)s-match-shear1diff-vs-mag_model.png">
        <img src="%(merun)s-%(serun)s-match-shear2diff-vs-mag_model.png">
        <img src="%(merun)s-%(serun)s-match-sizediff-vs-mag_model.png">

    </body>
</html>

    """ % {'merun':merun,'serun':serun}

    html_file=os.path.join(pd,'sheardiff.html')
    print 'writing html file:',html_file
    with open(html_file,'w') as fobj:
        fobj.write(html)

if len(sys.argv) < 3:
    print 'usage: plot-match-sheardiff.py merun serun'
    sys.exit(1)

merun=sys.argv[1]
serun=sys.argv[2]
d=deswl.files.collated_dir(merun)
pd=os.path.join(d,'sheardiff')
if not os.path.exists(pd):
    os.makedirs(pd) 

f='%(dir)s/match-%(serun)s-%(merun)s.dat' % {'dir':d,'merun':merun,'serun':serun}
f=os.path.expandvars(f)

print 'reading:',f
with recfile.Open(f,delim=' ',dtype=[('sid','i4'),('mid','i4')]) as fobj:
    matches = fobj.read()

print 'found:',matches.size,'matches'

mc = des.collate.open_columns(merun)
sc = des.collate.open_columns(serun)

mw = des.flags.select_good_me_bycol(mc)
sw = des.flags.select_good_se_bycol(sc)

print 'reading mag_model,shear_s2n,shear1,shear2,size from me'
data={}
data['mag_model_one'] = mc['mag_model'][mw]
data['shear_s2n_one'] = mc['shear_s2n'][mw]
data['mag_model'] = data['mag_model_one'][matches['mid']]
data['shear_s2n'] = data['shear_s2n_one'][matches['mid']]

data['mshear1'] = mc['shear1'][mw]
data['mshear2'] = mc['shear2'][mw]
data['msize'] = mc['shapelets_sigma'][mw]

print 'reading shear1,shear2,size from se'
data['sshear1'] = sc['shear1'][sw]
data['sshear2'] = sc['shear2'][sw]
data['ssize'] = sc['shapelets_sigma'][sw]



names={'shear1diff':r'$\langle \gamma_{1}^{me}-\gamma_{1}^{se}\rangle$',
       'shear2diff':r'$\langle \gamma_{2}^{me}-\gamma_{2}^{se}\rangle$',
       'sizediff':r'$\langle \sigma^{me}-\sigma^{se}\rangle$'}

dpi=90

nperbin=100000
for xfield in ['mag_model','shear_s2n']:
    if xfield == 'mag_model':
        xlog=False
        xmax=25
    else:
        xlog=True
        xmax=None

    for yfield in ['shear1','shear2','size']:
        diffname = yfield+'diff'
        if diffname not in data:
            mname='m'+yfield
            sname='s'+yfield
            data[diffname] = \
                data[mname][matches['mid']]-data[sname][matches['sid']]

        print 'plotting std(radiff) vs mag_model'
        plt = eu.plotting.bhist_vs(data, xfield, diffname, 
                                   names=names,
                                   max=xmax,
                                   nperbin=nperbin, 
                                   xlog=xlog)
        epsfile='%s-%s-match-%s-vs-%s.eps' % (merun,serun,diffname,xfield)
        epsfile=os.path.join(pd,epsfile)
        plt[0].write_eps(epsfile)
        converter.convert(epsfile, dpi=dpi, verbose=True)


make_html(pd,merun,serun)
