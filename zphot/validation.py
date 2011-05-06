import os
from sys import stdout,stderr
import glob
import pprint

import numpy
from numpy import where

import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import ahelp
from esutil.numpy_util import where1
from esutil.stat import histogram

import biggles
from biggles import FramedPlot, PlotKey
import converter

import zphot


class Validation:
    def __init__(self, wrun):
        """
        type='all' by default, meaning the combined file
        """

        raise ValueError("fix to work with new photo_dtype")

        self.wrun=wrun
        self.conf = zphot.read_config('weights',self.wrun)


        self.wt = zphot.weighting.WeightedTraining(self.wrun)
        self.types = self.wt.types




    def read_all(self):
        return read_training(self.allfile)
    

    def validation_dir(self):
        tmpname = self.wt.filename('all')
        validation_dir = os.path.dirname(tmpname)
        validation_dir = path_join(validation_dir,'validation')
        return validation_dir

    def maglim_dir(self):
        dir = self.validation_dir()
        dir = path_join(dir,'maglim')
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir

    def all_other_fdict(self,type):
        """
        Trying to recover a given training set using
        the everything set
        """
        dir = self.all_other_dir()

        fdict={}
        fdict['dir'] = dir

        fdict['nnear1'] = 5
        fdict['nnear2'] = 100

        fdict['trainfile'] = self.wt.filename('all')

        # this is the original file, it has the z values
        origfile = self.wt.filename(type)
        fdict['origfile'] = origfile

        # we need to write out origfile in the photofile format
        photofile = os.path.basename(origfile).replace('.dat','-asphoto.dat')
        photofile = path_join(dir, photofile)
        fdict['photofile'] = photofile

        tbase = os.path.basename(fdict['trainfile'])
        pbase = os.path.basename(fdict['photofile'])

        fdict['wfile1'] = \
            path_join(dir,tbase.replace('.dat','-%s-weights-5.dat' % type))
        fdict['wfile_nozero1'] = \
            path_join(dir,tbase.replace('.dat','-%s-weights-nozero-5.dat' % type))
        fdict['wfile2'] = \
            path_join(dir,tbase.replace('.dat','-%s-weights-100.dat' % type))
        fdict['numfile1'] = path_join(dir,pbase.replace('.dat','-num-5.dat'))
        fdict['numfile2'] = path_join(dir,pbase.replace('.dat','-num-100.dat'))

        fdict['wscript']=path_join(dir,type+'-doweights.sh')
        fdict['pofzscript']=path_join(dir,type+'-dopofz.sh')
        
        fdict['zhistfile'] = path_join(dir, pbase.replace('.dat','-zhist.eps') )
        return fdict


    def all_other_dir(self):
        dir = self.validation_dir()
        dir = path_join(dir,'all-other')
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir


    def prepare_all_other(self, type):
        """
        Split a sample into two and try to recover itself
        """

        fdict=self.all_other_fdict(type)
        orig = zphot.weighting.read_training(fdict['origfile'])
        
        dtype = zphot.weighting.photo_dtype()
        testset = numpy.zeros(orig.size, dtype=dtype)
        # NOTE this isn't really a photoid!
        testset['photoid'][:] = numpy.arange(orig.size,dtype='i4')
        testset['cmodelmag_dered_r'][:] = orig['cmodelmag_dered_r']
        testset['model_umg'][:] = orig['model_umg']
        testset['model_gmr'][:] = orig['model_gmr']
        testset['model_rmi'][:] = orig['model_rmi']
        testset['model_imz'][:] = orig['model_imz']

        print 'Writing:',fdict['photofile']
        r = eu.recfile.Open(fdict['photofile'],'w',delim=' ')
        r.write(testset)
        r.close()

        self.write_scripts(fdict)

    def compare_all_other(self, type, show=True):
        
        fdict=self.all_other_fdict(type)

        # this is the original file.  It has the redshifts
        orig = zphot.weighting.read_training(fdict['origfile'])

        # this is the outputs
        num = zphot.weighting.read_num(fdict['numfile1'])

        # this is the weights file
        weights = zphot.weighting.read_training(fdict['wfile2'])

        # recoverable set
        w_recoverable = where1(num['num'] > 0)
        # this is actually the indexes back into the "orig" file
        w_keep = num['photoid'][w_recoverable]

        # get the z values for these validation objects
        zrec = orig['z'][w_keep]

        binsize=0.0314
        valid_dict = histogram(zrec, min=0, max=1.1, binsize=binsize, more=True)
        plt=FramedPlot()

        vhist = valid_dict['hist']/(float(valid_dict['hist'].sum()))
        pvhist=biggles.Histogram(vhist, x0=valid_dict['low'][0], binsize=binsize)
        pvhist.label = 'truth'

        weights_dict = histogram(weights['z'], min=0, max=1.1, binsize=binsize,
                                 weights=weights['weight'], more=True)
        whist = weights_dict['whist']/weights_dict['whist'].sum()
        pwhist=biggles.Histogram(whist, x0=weights_dict['low'][0], 
                                 binsize=binsize, color='red')
        pwhist.label = 'weighted train'

        key = PlotKey(0.6,0.6,[pvhist,pwhist])
        plt.add(pvhist,pwhist,key)

        plt.add( biggles.PlotLabel(.8, .9, type) )

        plt.write_eps(fdict['zhistfile'])
        converter.convert(fdict['zhistfile'],dpi=90,verbose=True)
        if show:
            plt.show()



    def same_same_dir(self):
        dir = self.validation_dir()
        dir = path_join(dir,'same-same-test')
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir

    def same_same_fdict(self, type):
        dir = self.same_same_dir()

        fdict={}
        fdict['dir'] = dir

        origfile = self.wt.filename(type)
        fdict['origfile'] = origfile
        fdict['nnear1'] = 5
        fdict['nnear2'] = 100

        trainfile = os.path.basename(origfile).replace('.dat','-train.dat')
        validfile = os.path.basename(origfile).replace('.dat','-valid.dat')

        fdict['trainfile']=path_join(dir,trainfile)
        fdict['photofile']=path_join(dir,validfile)

        fdict['wfile1'] = fdict['trainfile'].replace('.dat','-weights-5.dat')
        fdict['wfile_nozero1'] = \
            fdict['trainfile'].replace('.dat','-weights-nozero-5.dat')
        fdict['wfile2'] = fdict['trainfile'].replace('.dat','-weights-100.dat')
        fdict['numfile1'] = fdict['photofile'].replace('.dat','-num-5.dat')
        fdict['numfile2'] = fdict['photofile'].replace('.dat','-num-100.dat')

        fdict['wscript']=path_join(dir,type+'-doweights.sh')
        fdict['pofzscript']=path_join(dir,type+'-dopofz.sh')
        
        fdict['zhistfile'] = fdict['photofile'].replace('.dat','-zhist.eps')

        return fdict

    def prepare_same_same(self, type):
        """
        Split a sample into two and try to recover itself
        """

        fdict=self.same_same_fdict(type)
        self.split_training(fdict['origfile'],fdict['trainfile'],fdict['photofile'])
        self.write_scripts(fdict)

    def split_training(self, origfile, trainfile, validfile):

        """
        Split the training weighting input into a training and validation
        set
        """

        data = zphot.weighting.read_training(origfile)

        # split in two
        nhalf = data.size/2

        ind1 = eu.numpy_util.random_subset(data.size, nhalf)

        # the other indices
        allind = numpy.arange(data.size, dtype='i4')
        allind[ind1] = -1
        ind2,=where(allind != -1)

        print ind1.size
        print ind2.size

        trainset=data[ind1]

        valid_dtype = zphot.weighting.photo_dtype()
        validset = numpy.zeros(ind2.size, dtype=valid_dtype)
        # NOTE this isn't really a photoid!
        validset['photoid'][:] = ind2
        validset['cmodelmag_dered_r'][:] = data['cmodelmag_dered_r'][ind2]
        validset['model_umg'][:] = data['model_umg'][ind2]
        validset['model_gmr'][:] = data['model_gmr'][ind2]
        validset['model_rmi'][:] = data['model_rmi'][ind2]
        validset['model_imz'][:] = data['model_imz'][ind2]

        print 'Writing:',trainfile
        r = eu.recfile.Open(trainfile,'w',delim=' ')
        r.write(trainset)
        r.close()
        print 'Writing:',validfile
        r = eu.recfile.Open(validfile,'w',delim=' ')
        r.write(validset)
        r.close()

    def compare_same_same(self, type, show=True):
        """
        Use the id from the validation set to go back and get the
        z for those objects.  Then plot histograms for comparision.

        read in all file
        read in validation set
            take recoverable subset based on num file
        Get z info for these points from the all file

        plot the histgram of actual validation set redshifts
        overplot the histgram of weighted redshifts

        Then bin by true validation set redshift and plot the
            ztrue - <z>
        Where <z> is the expectation value of z based on the p(z)
            <z> = integral( z*p(z) )/integral( p(z) )
        That will be noisy
        """
        
        fdict=self.same_same_fdict(type)

        # this is the original file
        all = zphot.weighting.read_training(fdict['origfile'])

        # this is the validation set, for which the "photoid" field
        # is actually an id pointing back into "all"
        # we take version 1 and will demand num > 0
        valid = zphot.weighting.read_photo(fdict['photofile'])
        num = zphot.weighting.read_num(fdict['numfile1'])


        # this is the weights file
        weights = zphot.weighting.read_training(fdict['wfile2'])

        # recoverable set
        w_recoverable = where1(num['num'] > 0)
        # this is actually the indexes back into the "all" file
        w_keep = num['photoid'][w_recoverable]

        # get the z values for these validation objects
        zvalid = all['z'][w_keep]

        binsize=0.0314
        valid_dict = histogram(zvalid, min=0, max=1.1, binsize=binsize, more=True)
        plt=FramedPlot()

        vhist = valid_dict['hist']/(float(valid_dict['hist'].sum()))
        pvhist=biggles.Histogram(vhist, x0=valid_dict['low'][0], binsize=binsize)
        pvhist.label = 'validation'

        weights_dict = histogram(weights['z'], min=0, max=1.1, binsize=binsize,
                                 weights=weights['weight'], more=True)
        whist = weights_dict['whist']/weights_dict['whist'].sum()
        pwhist=biggles.Histogram(whist, x0=weights_dict['low'][0], 
                                 binsize=binsize, color='red')
        pwhist.label = 'weighted train'

        key = PlotKey(0.6,0.6,[pvhist,pwhist])
        plt.add(pvhist,pwhist,key)

        plt.add( biggles.PlotLabel(.8, .9, type) )

        plt.write_eps(fdict['zhistfile'])
        converter.convert(fdict['zhistfile'],dpi=90,verbose=True)
        if show:
            plt.show()




    def write_scripts(self, fdict):

        weights_script_string="""
setup weighting -r ~/exports/weighting-work

trainfile="{trainfile}"
photofile="{photofile}"

nnear1={nnear1}
nnear2={nnear2}

wfile1="{wfile1}"
wfile_nozero1="{wfile_nozero1}"
wfile2="{wfile2}"

numfile1="{numfile1}"
numfile2="{numfile2}"

# first run

echo "
Running calcweights

First run with nnear=$nnear1
"

calcweights5         \\
    "$trainfile"    \\
    "$photofile"    \\
    "$nnear1"       \\
    "$wfile1" \\
    "$numfile1"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi


# remove training set objects with zero weight at first
# nnear
echo "
removing zero weight training objects to
    "$wfile_nozero1"
"
awk '$3>0' < "$wfile1" > "$wfile_nozero1"

if [ "$?" != "0" ]; then
    echo Error running awk.  Halting
    exit 45
fi

# second run
echo "
Second run with nnear=$nnear2
"

calcweights5                \\
    "$wfile_nozero1" \\
    "$photofile"           \\
    "$nnear2"              \\
    "$wfile2"        \\
    "$numfile2"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi
        """.format(photofile=fdict['photofile'],
                   trainfile=fdict['trainfile'],
                   nnear1=fdict['nnear1'],
                   nnear2=fdict['nnear2'],
                   wfile1=fdict['wfile1'],
                   wfile_nozero1=fdict['wfile_nozero1'],
                   wfile2=fdict['wfile2'],
                   numfile1=fdict['numfile1'],
                   numfile2=fdict['numfile2'])
    
        stdout.write("Writing weights script: '%s'\n" % fdict['wscript'])
        wf = open(fdict['wscript'],'w')
        wf.write(weights_script_string)
        wf.close()

def prepare_all_other(wrun):
    wt = zphot.weighting.WeightedTraining(wrun)
    types = wt.types
    v = Validation(wrun)
    for type in types:
        v.prepare_all_other(type)
 
def compare_all_other(wrun):
    wt = zphot.weighting.WeightedTraining(wrun)
    types = wt.types
    v = Validation(wrun)
    for type in types:
        v.compare_all_other(type,show=False)
 

def prepare_all_same_same(wrun):
    wt = zphot.weighting.WeightedTraining(wrun)
    types = wt.types
    types.append('all')
    v = Validation(wrun)
    for type in types:
        v.prepare_same_same(type)
    
def compare_all_same_same(wrun):
    wt = zphot.weighting.WeightedTraining(wrun)
    types = wt.types
    types.append('all')
    v = Validation(wrun)
    for type in types:
        if type not in ['cfrs','tkrs-fix']:
            v.compare_same_same(type,show=False)
 

