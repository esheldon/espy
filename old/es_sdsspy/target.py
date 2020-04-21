import os
import glob
from sys import stdout

import numpy

import esutil as eu
from esutil.numpy_util import where1
import sweeps

import re

class ThingIdAdder:

    def load_data(self):
        print 'opening columns'

        cdict = {}

        cdict['gal']  = {}
        cdict['gal']['cols'] = sweeps.open_columns('gal')
        cdict['star']  = {}
        cdict['star']['cols'] = sweeps.open_columns('star')

        for type in ['gal','star']:
            cols = cdict[type]['cols']

            print '%s dir: %s' % (type, cols.dir)

            print 'getting ra'
            cdict[type]['ra'] = cols['ra'][:]
            print 'getting dec'
            cdict[type]['dec'] = cols['dec'][:]
            print 'getting thing_id'
            cdict[type]['thing_id'] = cols['thing_id'][:]

            print 'getting htmid10'
            htmid10 = cols['htmid10'][:]
            cdict[type]['htmid10'] = htmid10

            cdict[type]['minid'] = htmid10.min()
            cdict[type]['maxid'] = htmid10.max()

            stdout.write("Getting reverse indices\n");stdout.flush()
            hist, htmrev = eu.stat.histogram(htmid10-cdict[type]['minid'],rev=True)

            cdict[type]['htmrev']=htmrev

        self.cdict = cdict
        self.radius = 2.0/3600.0

        self.htm = eu.htm.HTM(10)

    def match_radec(self, ra, dec, type):
        cd = self.cdict[type]
        mdata,mall,d = self.htm.match(ra, dec,
                                      cd['ra'],
                                      cd['dec'],
                                      self.radius,
                                      htmrev2=cd['htmrev'], 
                                      minid=cd['minid'],
                                      maxid=cd['maxid'],
                                      verbose=True)
        return mdata, mall, d


    def match(self, filename, primary_type, clobber=False):
        if primary_type == 'gal':
            secondary_type = 'star'
        else:
            secondary_type = 'gal'

        if filename.find('dr8matchboth') != -1:
            raise ValueError("dr8matchboth is already in filename")


        new_filename = filename.replace('.fits','-dr8matchboth.fits')
        if os.path.exists(new_filename):
            #raise ValueError("Output already exists: '%s'" % new_filename)
            if not clobber:
                print "Output already exists:",new_filename
                print "returning"
                return
            else:
                print 'Removing existing file:',new_filename
                os.remove(new_filename)

        print "Will write to file:",new_filename

        data1 = eu.io.read(filename, ext=1, lower=True, verbose=True)
        data2 = eu.io.read(filename, ext=2, lower=True, verbose=True)

        if 'thing_id' in data2.dtype.names:
            psweep = data1['photo_sweep'][0]
            if psweep.find('dr8') != -1:
                print 'thing_id is already correct, not creating output'
                return

        output = eu.numpy_util.add_fields(data2, [('dr8_thing_id','i4'),
                                                  ('dr8_ra','f8'),
                                                  ('dr8_dec','f8'),
                                                  ('dr8_match_type','S4')])
        output['dr8_thing_id'][:] = -1
        output['dr8_ra'][:] = -9999
        output['dr8_dec'][:] = -9999
        output['dr8_match_type'] = 'None'

        if not hasattr(self, 'cdict'):
            self.load_data()
        
        mdata,mall,d = self.match_radec(data2['ra'],data2['dec'],primary_type)
        print "matched: %s/%s to '%s'" % (mdata.size,data2.size,primary_type)

        cd = self.cdict[primary_type]
        output['dr8_thing_id'][mdata]   = cd['thing_id'][mall]
        output['dr8_ra'][mdata]         = cd['ra'][mall]
        output['dr8_dec'][mdata]        = cd['dec'][mall]
        output['dr8_match_type'][mdata] = primary_type

        if mdata.size != data2.size:
            w=where1(output['dr8_match_type'] == 'None')
            if w.size != (data2.size - mdata.size):
                raise ValueError("Expected %s unmatched" % (data2.size-mdata.size))
            print "matching remaining %s to secondary set: '%s'" \
                    % (w.size, secondary_type)
            mdata,mall,d = \
                self.match_radec(data2['ra'][w],data2['dec'][w],secondary_type)
            print "matched: %s/%s to '%s'" % (mdata.size,w.size,secondary_type)

            if mdata.size != 0:
                cd2 = self.cdict[secondary_type]
                mdata = w[mdata]
                output['dr8_thing_id'][mdata]   = cd2['thing_id'][mall]
                output['dr8_ra'][mdata]         = cd2['ra'][mall]
                output['dr8_dec'][mdata]        = cd2['dec'][mall]
                output['dr8_match_type'][mdata] = secondary_type
    

        #raise ValueError("stopping for debugging")
        eu.io.write(new_filename, data1, verbose=True)
        eu.io.write(new_filename, output, verbose=True, append=True)



def add_thingid_all():

    tadd = ThingIdAdder()

    dir = os.environ['BOSS_ROOT']
    dir = os.path.join(dir, 'target')

    subdirs = ['vcat-2010-07-02','comm','comm2',
               'main001','main002','main003','main004',
               'main005','main006','main007']

    for d in subdirs:
        indir = os.path.join(dir, d)
        pattern = os.path.join(indir, '*collate*.fits')
        files = glob.glob(pattern)

        #files = [f for f in files if f.find('dr8match') == -1]

        for f in files:
            # make sure it isn't one of the collate-002 type files
            if f.find('dr8match') == -1 \
                    and re.match('.*collate-.*[0-9][0-9][0-9].*',f) == None \
                    and re.match('.*bitmask.*',f) == None \
                    and re.match('.*inchunk.*',f) == None:
                #print f

                if f.find('lrg') != -1:
                    primary_type = 'gal'
                elif f.find('qso') != -1:
                    primary_type = 'star'
                elif f.find('std') != -1:
                    primary_type = 'star'
                else:
                    raise ValueError("Expected 'lrg','qso','std' in name: '%s'" % f)

                print 'Processing file:',f
                tadd.match(f, primary_type)


def check_add_thingid_all(filename=None):

    if filename is None:
        fobj = stdout
    else:
        fobj = open(filename,'w')

    dir = os.environ['BOSS_ROOT']
    dir = os.path.join(dir, 'target')

    subdirs = ['vcat-2010-07-02','comm','comm2',
               'main001','main002','main003','main004',
               'main005','main006','main007']

    for d in subdirs:
        indir = os.path.join(dir, d)
        pattern = os.path.join(indir, '*dr8matchboth*.fits')
        files = glob.glob(pattern)

        for f in files:
            data = eu.io.read(f,ext=2)

            w=where1(data['dr8_thing_id'] >= 0)
            nmiss = data.size-w.size
            perc_miss = float(nmiss)/data.size*100
            fb = os.path.basename(f)
            fb = os.path.join(d, fb)
            fobj.write('%-75s miss: %5i / %7i (%0.2f%%)' % (fb,nmiss,data.size,perc_miss))

            thing_id = data['dr8_thing_id']
            w=where1(thing_id >= 0)
            thing_id = thing_id[w]

            uid = numpy.unique(thing_id)
            fobj.write("  dups: %s\n" % (thing_id.size-uid.size,) )

    if filename is not None:
        fobj.close()

