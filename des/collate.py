from __future__ import print_function
import os
import sys
from sys import stdout,stderr
from copy import copy

import esutil
from esutil import json_util
from esutil import numpy_util
from esutil.plotting import setuplot, bwhiskers
from esutil.ostools import path_join, getenv_check, expand_path

import deswl
from deswl import wlpipe

import numpy
from numpy import where

import columns


class ColumnCollator: 
    def __init__(self, serun, split=None):
        self.serun = serun
        self.split=split

    def collate(self):
        self.load_columns_for_writing()

        self.load_flist()

        ntot = len(self.flist)
        i=1
        for fdict in self.flist:
            stdout.write('-'*70)
            stdout.write("\nProcessing %d/%d\n" % (i,ntot))
            data={}
            data['exposurename']=fdict['exposurename']
            data['ccd']=int(fdict['ccd'])

            # this will add data to the dictionary
            self.load_data(data)
            self.write_data(data)
            i+=1

    def write_data(self, data):
        self.write_from_stars(data)
        self.write_from_shear(data)
        self.write_from_psf(data)



    def load_data(self, data):
        
        if self.split not in [None,1,2]:
            raise ValueError("split must be None,1 or 2, found %s" % self.split)

        if not hasattr(self,'uid_offset'):
            # start id counter
            self.uid_offset = 0

        tdict = {}
        tdict['stars'] = 'stars'
        tdict['psf'] = 'psf'
        tdict['shear'] = 'shear'
        
        if self.split is not None:
            tdict['psf'] += str(self.split)
            tdict['shear'] += str(self.split)

        for type in tdict:
            fname=tdict[type]
            data[type] = deswl.files.se_read(data['exposurename'], 
                                             data['ccd'], 
                                             fname, 
                                             serun=self.serun, 
                                             ext=1, 
                                             ensure_native=True,
                                             verbose=True)

        data['idvals'] = self.uid_offset + numpy.arange(data['stars'].size,dtype='i4')
        self.uid_offset += data['stars'].size


    def coldir(self, fits=False):
        coldir = deswl.files.se_coldir(self.serun, fits=fits)
        coldir = expand_path(coldir)

        if self.split is not None:
            if not os.path.exists(coldir):
                raise RuntimeError("Coldir must exist before doing split")
            main_cols = columns.Columns(coldir)

            # in this case we need the columns dir to exist
            splitname = 'split%s' % self.split
            if fits:
                coldir = path_join(coldir,splitname+'-fits')
            else:
                coldir = path_join(coldir,splitname+'.cols')
        
        coldir = expand_path(coldir)
        return coldir


    def load_columns_for_writing(self):
        coldir=self.coldir()

        stdout.write('Will write to coldir %s\n' % coldir)
        if os.path.exists(coldir):
            raise RuntimeError("Coldir exists. Please start from scratch\n")

        self.cols = columns.Columns(coldir)
        self.cols.create()

        # load the sub columns dir for PSF stars
        psf_subdir = os.path.join(coldir, 'psfstars.cols')
        self.cols.load_coldir(psf_subdir)


    def load_flist(self):
        flistfile=deswl.files.se_collated_path(self.serun,'goodlist')
        self.flist=esutil.io.read(flistfile, verbose=True)


    def write_from_stars(self, data):

        cols=self.cols

        # make copies of the data before writing to the columns.  This is
        # to force the data type, including native endianness

        id = numpy.array(data['stars']['id'], dtype='i2')
        cols.write_column('id',id)
        del id

        x=numpy.array(data['stars']['x'], dtype='f4')
        cols.write_column('x',x)
        del x

        y=numpy.array(data['stars']['y'], dtype='f4')
        cols.write_column('y',y)
        del y

        sky=numpy.array(data['stars']['sky'], dtype='f4')
        cols.write_column('sky',sky)
        del sky


        size_flags=numpy.array(data['stars']['size_flags'], dtype='i2')
        cols.write_column('size_flags',size_flags)
        del size_flags


        sigma0=numpy.array(data['stars']['sigma0'], dtype='f4')
        cols.write_column('sigma0',sigma0)
        del sigma0


        cols.write_column('imag',data['stars']['mag'])


        star_flag=numpy.array(data['stars']['star_flag'], dtype='i1')
        cols.write_column('star_flag',star_flag)
        del star_flag

        if 'sg' in data['stars'].dtype.names:
            sg = numpy.array(data['stars']['sg'],dtype='f4')
            cols.write_column('sg',sg)
            del sg


    def write_from_shear(self, data):

        cols = self.cols

        ra = numpy.array(data['shear']['ra'], dtype='f8')
        cols.write_column('ra',ra)
        dec = numpy.array(data['shear']['dec'], dtype='f8')
        cols.write_column('dec',dec)

        gal_order=numpy.array(data['shear']['gal_order'], dtype='i1')
        cols.write_column('shapelets_order',gal_order)
        del gal_order

        cols.write_column('shear_flags',data['shear']['shear_flags'])

        shear1=numpy.array(data['shear']['shear1'], dtype='f4')
        cols.write_column('shear1',shear1)
        del shear1
        shear2=numpy.array(data['shear']['shear2'], dtype='f4')
        cols.write_column('shear2',shear2)
        del shear2

        shear_cov11=numpy.array(data['shear']['shear_cov00'], dtype='f4')
        cols.write_column('shear_cov11',shear_cov11)
        del shear_cov11
        shear_cov12=numpy.array(data['shear']['shear_cov01'], dtype='f4')
        cols.write_column('shear_cov12',shear_cov12)
        del shear_cov12
        shear_cov22=numpy.array(data['shear']['shear_cov11'], dtype='f4')
        cols.write_column('shear_cov22',shear_cov22)
        del shear_cov22

        shapelet_sigma=numpy.array(data['shear']['shapelet_sigma'], dtype='f4')
        cols.write_column('shapelets_sigma',shapelet_sigma)
        del shapelet_sigma


        shapelets_prepsf=numpy.array(data['shear']['shapelets_prepsf'], dtype='f4')
        cols.write_column('shapelets_prepsf',shapelets_prepsf)


        e1 = shapelets_prepsf[:,3]*numpy.sqrt(2)
        e2 = -shapelets_prepsf[:,4]*numpy.sqrt(2)
        del shapelets_prepsf

        cols.write_column('e1',e1)
        del e1
        cols.write_column('e2',e2)
        del e2

        # shapelets info interpolated from PSF stars
        psf_order=numpy.array(data['shear']['interp_psf_order'], dtype='i1')
        cols.write_column('interp_psf_order',psf_order)
        del psf_order


        psf_sigma=numpy.array(data['shear']['interp_psf_sigma'], dtype='f4')
        cols.write_column('interp_psf_sigma',psf_sigma)
        del psf_sigma

        psf_coeffs=numpy.array(data['shear']['interp_psf_coeffs'], dtype='f4')
        cols.write_column('interp_psf_shapelets',psf_coeffs)

        # get e1/e2 from interpolated coeffs
        psf_e1 = psf_coeffs[:,3]*numpy.sqrt(2)
        psf_e2 = -psf_coeffs[:,4]*numpy.sqrt(2)

        cols.write_column('interp_psf_e1',psf_e1)
        cols.write_column('interp_psf_e2',psf_e2)

        del psf_coeffs
        del psf_e1
        del psf_e2



        # nu, the S/N.  Error on either component e1,e2 is sqrt(2)/nu
        if 'nu' in data['shear'].dtype.names:
            nu=numpy.array(data['shear']['nu'],dtype='f4')
        elif 'shear_signal_to_noise' in data['shear'].dtype.names:
            nu=numpy.array(data['shear']['shear_signal_to_noise'],dtype='f4')
        else:
            raise ValueError('neither nu or shear_signal_to_noise found')
        cols.write_column('shear_s2n',nu)



        # our id column
        cols.write_column('uid', data['idvals'])

        # Same for all objects in this set
        ccd_array = numpy.zeros(data['stars'].size, dtype='i1')
        ccd_array[:] = data['ccd']
        cols.write_column('ccd',ccd_array)

        exp_array = numpy.zeros(data['stars'].size, dtype='S20')
        exp_array[:] = data['exposurename']
        cols.write_column('exposurename',exp_array)

        del ccd_array
        del exp_array

    def write_from_psf(self, data):
        cols = self.cols

        psf = data['psf']
        stars = data['stars']

        cols['psfstars'].write_column('psf_flags', psf['psf_flags'])

        if 'nu' in psf.dtype.names:
            nu=numpy.array(psf['nu'],dtype='f4')
        elif 'psf_signal_to_noise' in psf.dtype.names:
            nu=numpy.array(psf['psf_signal_to_noise'],dtype='f4')
        else:
            raise ValueError('no nu or psf_signal_to_noise in psf file')
        cols['psfstars'].write_column('psf_s2n', nu)

        order=numpy.array(psf['psf_order'],dtype='i1')
        cols['psfstars'].write_column('shapelets_order', order)

        sigma_p = numpy.array(psf['sigma_p'], dtype='f4')
        cols['psfstars'].write_column('shapelets_sigma', sigma_p)

        psf_shapelets = numpy.array(psf['shapelets'], dtype='f4')
        cols['psfstars'].write_column('shapelets', psf_shapelets)


        e1 = psf_shapelets[:,3]*numpy.sqrt(2)
        e2 = -psf_shapelets[:,4]*numpy.sqrt(2)

        cols['psfstars'].write_column('e1',e1)
        cols['psfstars'].write_column('e2',e2)


        pind, sind = esutil.numpy_util.match(psf['id'], stars['id'])

        if pind.size != psf.size:
            raise ValueError("Some PSF stars did not match\n")
        p_idvals = data['idvals'][sind]
        cols['psfstars'].write_column('uid', p_idvals)

        del nu
        del order
        del sigma_p
        del psf_shapelets
        del p_idvals
        del e1
        del e2

    def convert2fits(self):
        """
        Just make fits versions of all the files
        """

        coldir = self.coldir()
        fitsdir = self.coldir(fits=True)

        psfstars_fitsdir = os.path.join(fitsdir,'psfstars-fits')
        if not os.path.exists(fitsdir):
            os.makedirs(fitsdir)
            os.makedirs(psfstars_fitsdir)

        print('coldir:',coldir)
        print('fitsdir:',fitsdir)
        print('psfstars_fitsdir:',psfstars_fitsdir)

        cols = columns.Columns(coldir)

        print(cols)

        for col in sorted(cols):
            if col in ['shapelets_prepsf','interp_psf_shapelets']:
                print("skipping large column:",col)
            else:
                print('column: ',col)

                if col == 'psfstars':
                    continue

                fname = cols[col].filename

                data = esutil.sfile.read(fname)

                fitsfile = os.path.join(fitsdir, col+'.fits')
                print(fitsfile)

                # don't make a copy!
                esutil.io.write(fitsfile, data, clobber=True,copy=False)

                del data

        psfstars_cols = cols['psfstars']
        for col in sorted(psfstars_cols):
            print('column: ',col)

            fname = psfstars_cols[col].filename

            data = esutil.sfile.read(fname)

            fitsfile = os.path.join(psfstars_fitsdir, col+'.fits')
            print(fitsfile)

            # don't make a copy!
            esutil.io.write(fitsfile, data, clobber=True,copy=False)

            del data





def create_se_shear_columns_indexes(serun, coldir=None):

    if coldir is None:
        coldir = deswl.files.se_coldir(serun)
    coldir = expand_path(coldir)

    cols = columns.Columns(coldir,verbose=True)
    if not cols.dir_exists():
        raise RuntimeError("You haven't added the data yet")

    # create some indexes
    # after this, data automatically gets added to the index
    cols['ccd'].create_index()
    cols['size_flags'].create_index()
    cols['star_flag'].create_index()
    cols['shear_flags'].create_index()
    cols['exposurename'].create_index()


    cols['psfstars']['psf_flags'].create_index()
    cols['psfstars']['uid'].create_index()



def collate_se_shear_columns(serun, split=False,
                             coldir=None, limit=None):
    """
    """

    # start id counter
    uid = 0

    rc=deswl.files.Runconfig(serun)
    header = rc.asdict()
    header['serun'] = header['run']
    del header['run']

    if coldir is None:
        coldir = deswl.files.se_coldir(serun, split=split)
    coldir = expand_path(coldir)

    stdout.write('Will write to coldir %s\n' % coldir)

    if os.path.exists(coldir):
        raise RuntimeError("Coldir exists. Please start from scratch\n")


    flistfile=deswl.files.se_collated_path(serun,'goodlist')
    stdout.write("Reading goodlist: %s\n" % flistfile)
    flist=json_util.read(flistfile)

    cols = columns.Columns(coldir)
    cols.create()

    # load the sub columns dir for PSF stars
    psf_subcoldir = os.path.join(coldir, 'psfstars.cols')
    cols.load_coldir(psf_subcoldir)

    ntot=len(flist)
    i=1
    for fdict in flist:

        stdout.write('-'*70)
        stdout.write('\nReading: %d/%d\n' % (i,ntot))
        data = {}

        exposurename=fdict['exposurename']
        ccd=int(fdict['ccd'])

        data['exposurename'] = exposurename
        data['ccd'] = ccd

        data['stars'] = deswl.files.se_read(exposurename, ccd, 'stars', 
                                              serun=serun, ext=1, verbose=True)

        data['idvals'] = uid + numpy.arange(data['stars'].size,dtype='i4')

        data['shear'] = deswl.files.se_read(exposurename, ccd, 'shear', 
                                              serun=serun, ext=1, verbose=True)

        if data['stars'].size != data['shear'].size:
            raise ValueError('Size file and shears file must be same length')
        w,=numpy.where(data['stars']['id'] != data['shear']['id'])
        if w.size > 0:
            raise ValueError("id fields don't match up")

        
        # write main columns
        _do_collate_se_shear_columns(cols, data)


        # write psf columns
        data['psf'] = deswl.files.se_read(exposurename, ccd, 'psf', 
                                            serun=serun, ext=1,verbose=True)

        _do_collate_se_shear_columns_psf(cols, data)

        stdout.flush()
        
        i += 1
        uid += stars.size

        if limit is not None:
            if i > limit:
                return


def _do_collate_se_shear_columns(cols, stars, shear, idvals, ccd, exposurename):
    # note I'm using some smaller types here to save space

    # data from the stars structure
    id = numpy.array(stars['id'], dtype='i2')
    cols.write_column('id',id)
    del id

    x=numpy.array(stars['x'], dtype='f4')
    cols.write_column('x',x)
    del x

    y=numpy.array(stars['y'], dtype='f4')
    cols.write_column('y',y)
    del y

    sky=numpy.array(stars['sky'], dtype='f4')
    cols.write_column('sky',sky)
    del sky


    size_flags=numpy.array(stars['size_flags'], dtype='i2')
    cols.write_column('size_flags',size_flags)
    del size_flags


    sigma0=numpy.array(stars['sigma0'], dtype='f4')
    cols.write_column('sigma0',sigma0)
    del sigma0


    cols.write_column('imag',stars['mag'])


    star_flag=numpy.array(stars['star_flag'], dtype='i1')
    cols.write_column('star_flag',star_flag)
    del star_flag


    # Data from the shear structure
    cols.write_column('shear_flags',shear['shear_flags'])
    cols.write_column('ra',shear['ra'])
    cols.write_column('dec',shear['dec'])

    shear1=numpy.array(shear['shear1'], dtype='f4')
    cols.write_column('shear1',shear1)
    del shear1
    shear2=numpy.array(shear['shear2'], dtype='f4')
    cols.write_column('shear2',shear2)
    del shear2

    shear_cov00=numpy.array(shear['shear_cov00'], dtype='f4')
    cols.write_column('shear_cov00',shear_cov00)
    del shear_cov00
    shear_cov01=numpy.array(shear['shear_cov01'], dtype='f4')
    cols.write_column('shear_cov01',shear_cov01)
    del shear_cov01
    shear_cov11=numpy.array(shear['shear_cov11'], dtype='f4')
    cols.write_column('shear_cov11',shear_cov11)
    del shear_cov11

    gal_order=numpy.array(shear['gal_order'], dtype='i1')
    cols.write_column('shapelets_order',gal_order)
    del gal_order


    shapelet_sigma=numpy.array(shear['shapelet_sigma'], dtype='f4')
    cols.write_column('shapelets_sigma',shapelet_sigma)
    del shapelet_sigma


    shapelets_prepsf=numpy.array(shear['shapelets_prepsf'], dtype='f4')
    cols.write_column('shapelets_prepsf',shapelets_prepsf)


    e1 = shapelets_prepsf[:,3]*numpy.sqrt(2)
    e2 = -shapelets_prepsf[:,4]*numpy.sqrt(2)
    del shapelets_prepsf

    cols.write_column('e1',e1)
    del e1
    cols.write_column('e2',e2)
    del e2

    # shapelets info interpolated from PSF stars
    psf_order=numpy.array(shear['interp_psf_order'], dtype='i1')
    cols.write_column('interp_psf_order',psf_order)
    del psf_order


    psf_sigma=numpy.array(shear['interp_psf_sigma'], dtype='f4')
    cols.write_column('interp_psf_sigma',psf_sigma)
    del psf_sigma

    psf_coeffs=numpy.array(shear['interp_psf_coeffs'], dtype='f4')
    cols.write_column('interp_psf_shapelets',psf_coeffs)

    # get e1/e2 from interpolated coeffs
    psf_e1 = psf_coeffs[:,3]*numpy.sqrt(2)
    psf_e2 = -psf_coeffs[:,4]*numpy.sqrt(2)

    cols.write_column('interp_psf_e1',psf_e1)
    cols.write_column('interp_psf_e2',psf_e2)

    del psf_coeffs
    del psf_e1
    del psf_e2



    # nu, the S/N.  Error on either component e1,e2 is sqrt(2)/nu
    nu=numpy.array(shear['nu'],dtype='f4')
    cols.write_column('nu',nu)



    # our id column
    cols.write_column('uid', idvals)

    # Same for all objects in this set
    ccd_array = numpy.zeros(stars.size, dtype='i1')
    ccd_array[:] = ccd
    cols.write_column('ccd',ccd_array)

    exp_array = numpy.zeros(stars.size, dtype='S20')
    exp_array[:] = exposurename
    cols.write_column('exposurename',exp_array)

    del ccd_array
    del exp_array


def _do_collate_se_shear_columns_psf(cols, stars, psf, idvals):
    cols['psfstars'].write_column('psf_flags', psf['psf_flags'])

    nu=numpy.array(psf['nu'],dtype='f4')
    cols['psfstars'].write_column('nu', nu)

    order=numpy.array(psf['psf_order'],dtype='i1')
    cols['psfstars'].write_column('shapelets_order', order)

    sigma_p = numpy.array(psf['sigma_p'], dtype='f4')
    cols['psfstars'].write_column('shapelets_sigma', sigma_p)

    psf_shapelets = numpy.array(psf['shapelets'], dtype='f4')
    cols['psfstars'].write_column('shapelets', psf_shapelets)


    e1 = psf_shapelets[:,3]*numpy.sqrt(2)
    e2 = -psf_shapelets[:,4]*numpy.sqrt(2)

    cols['psfstars'].write_column('e1',e1)
    cols['psfstars'].write_column('e2',e2)


    pind, sind = esutil.numpy_util.match(psf['id'], stars['id'])

    if pind.size != psf.size:
        raise ValueError("Some PSF stars did not match\n")
    p_idvals = idvals[sind]
    cols['psfstars'].write_column('uid', p_idvals)

    del nu
    del order
    del sigma_p
    del psf_shapelets
    del p_idvals
    del e1
    del e2




def collate_se_shear_array(num=1, fields=None, allfields=False):
    """
    This is just a minimal structure

    Note I'm using f4 instead of f8 here
    """
    dt=[('x','f4'),
        ('y','f4'),
        ('ra','f8'),
        ('dec','f8'),
        ('exposurename','S20'),
        ('ccd','i1'),
        ('size_flags','i4'),
        ('magi','f4'),
        ('sigma0','f4'),
        ('star_flag','i4'),
        ('shear_flags','i4'),
        ('shapelet_sigma','f4'),
        ('shear1','f4'),
        ('shear2','f4'),
        ('shear_cov00','f4'),
        ('shear_cov01','f4'),
        ('shear_cov11','f4'),
        ('shapelets_prepsf','f4',28)]

    dnames=[d[0] for d in dt]


    if fields is None:
        if allfields:
            fields = dnames
        else:
            fields = [f[0] for f in dt if f[0] != 'shapelets_prepsf']

    # preserve order
    newdt = []
    for f in fields:
        try:
            ind = dnames.index(f)
        except:
            raise ValueError,"Field '%s' not available" % f

        newdt.append(dt[ind]) 
    dt=newdt

    return numpy.zeros(num, dtype=dt)


def collate_se_shear(serun, objclass='all',
                     fields=None,
                     allfields=False,
                     clobber=True,
                     outdir=None,
                     convert_to_degrees=False,
                     limit=None):

    """
    """



    rc=deswl.files.Runconfig(serun)
    header = rc.asdict()
    header['serun'] = header['run']
    del header['run']

    flistfile=deswl.files.se_collated_path(serun,'goodlist')
    stdout.write("Reading goodlist: %s\n" % flistfile)
    flist=json_util.read(flistfile)

    if outdir is not None:
        outdir=esutil.ostools.expand_path(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    outfile = deswl.files.se_collated_path(serun, objclass, ftype='rec', 
                                             dir=outdir)
    stdout.write('Will write to %s\n' % outfile)

    odir=os.path.dirname(outfile)
    if not os.path.exists(odir):
        stdout.write("Creating output directory: %s\n" % odir)
        os.makedirs(odir)

    ename=esutil.ostools.expand_path(outfile)
    if os.path.exists(ename):
        if not clobber:
            raise RuntimeError("Outfile exists. Send clobber=True to "
                               "overwrite")

    sf = esutil.sfile.Open(outfile, mode='w')

    ntot=len(flist)
    i=1
    for fdict in flist:

        exposurename=fdict['exposurename']
        ccd=int(fdict['ccd'])

        #if True:
        if (i % 10) == 0:
            stars_file=deswl.files.se_path(exposurename,ccd,'stars',
                                             serun=serun)
            shear_file=deswl.files.se_path(exposurename,ccd,'shear',
                                             serun=serun)
            stdout.write('Reading: %d/%d\n\t%s\n' % (i,ntot,stars_file))

        stars = deswl.files.se_read(exposurename, ccd, 'stars', 
                                      serun=serun, ext=1)

        if (i % 10) == 0:
            stdout.write('\t%s\n' % shear_file)
        shear = deswl.files.se_read(exposurename, ccd, 'shear', 
                                      serun=serun, ext=1)

        # we output ra/dec in arcsec!?!?!
        if convert_to_degrees:
            deswl.files.convert_to_degrees(shear)

        if stars.size != shear.size:
            raise ValueError('Size file and shears file must be same length')
        w,=numpy.where(stars['id'] != shear['id'])
        if w.size > 0:
            raise ValueError("id fields don't match up")

        outarr=collate_se_shear_array(stars.size, 
                                      fields=fields, 
                                      allfields=allfields)

        numpy_util.copy_fields(stars, outarr)
        numpy_util.copy_fields(shear, outarr)

        # renaming
        outarr['magi'] = stars['mag']

        outarr['ccd'] = ccd
        outarr['exposurename'] = exposurename

        sf.write(outarr, header=header)
        i += 1

        if limit is not None:
            if i > limit:
                return

    sf.close()
    stdout.write('Wrote to %s\n' % outfile)


def write_se_collate_html(serun, showsplit=False):
    """

    Need to make collation band-aware!

    """

    rc = deswl.files.Runconfig(serun)

    dir=deswl.files.se_collated_dir(serun)
    html_dir = path_join(dir, 'html')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)

    fitsdir=path_join(dir,serun+'-fits')
    coldir=path_join(dir,serun+'.cols')
    if os.path.exists(fitsdir):
        table_comment='<i>Note: The collated FITS files do not include the large columns <code>shapelets_prepsf</code> and <code>interp_psf_shapelets</code></i>'
    else:
        table_comment='<i><b>FITS Not Yet Available</b></i>'
    if not os.path.exists(coldir):
        col_comment += ' <i><b>Columns database not yet available</b></i>'


    collate_dir= deswl.files.se_collated_dir(serun)
    plotdir=path_join(collate_dir,'plots')

    band='i'
    tfile= path_join(plotdir, '%s-%s-sizemag.png' % (serun,band))
    if os.path.exists(tfile):
        qatext = """
		<td>
			<p>
			<p>
			<h2>QA and Analysis</h2>
			Here are  some <a href="html/qa.html">QA plots  and tests</a>.
			<p>
			<p>
		</td>
        """
    else:
        qatext=""



    if showsplit:
        splitinfo="NOTE: <b>Validation splits 1/2 are in the %s-fits/split* and %s.cols/split* subdirectories</b>" % (serun,serun)
    else:
        splitinfo=""

    dc4_shear_byreg="""
	<tr>
		<td>
			<p>
			<p>
			<h2>Shear Values in DC4 Simulations</h2>
			The shear was input in pixel coordinates.  Note X->DEC and Y->--RA.  The following table
			shows shear1/shear2 values in pixel and RA/DEC coords.


			<pre>
Region RA          Dec        gamma1     gamma2   gamma1     gamma2    posangle
                               (xy)       (xy)   (rad/dec)  (ra/dec)   (ra/dec)
--------------------------------------------------------------------------------
1      >= 335 deg  >= -25 deg     0.05    0        0         0.05       45.0
2      >= 335      <  -25        -0.025   0.025    0.025    -0.025     -22.5
3      <  335      >= -25         0.05    0.05    -0.05     -0.05      -67.5
4      <  335      <  -25         0      -0.025    0         0.025      45.0
			</pre>


		</td>

	</tr>
    """

    shear_desc=""
    if rc['dataset'].lower() == 'dc4':
        shear_desc = dc4_shear_byreg

    download="""
<html>
<head>
	<link rel="stylesheet" href="html/table.css" type="text/css">
</head>

<body bgcolor=white>

<h1>{dataset} Single Epoch Catalogs: Run {serun}</h1>

<table width="500px">
	<tr>
		<td>

            <h2>Raw WL outputs</h2>
            Raw outputs can be found <a href="..">Here</a>.   If you plan to download
            <i><b>all the data</b></i>, please contact Erin Sheldon.
            <p>
            These files are located on the <a href="https://sites.google.com/site/bnlwlens/software-tutorials/software-setup">astro cluster</a> at
            <pre>
    /astro/u/astrodat/data/DES/wlbnl/{serun}
            </pre>

			<h2>Collated WL outputs</h2>
            The data are broken up into two "tables": A main table for every object and a sub-table for
            PSF stars.
            <ul>
                <li><a href="html/all-columns-descr.html">Description of main table</a>
                <li><a href="html/psfstars-columns-descr.html">Description of PSF stars table</a>
            </ul>

            The structure of these "tables" is actually quite simple: There is a file for each
            column.  The PSF star sub-table is just a sub-directory with more column files.
            This table can be connected back to the main table using the uid column, which points
            to the row number in the main table.

            <p>
            <table>
                <table class=simple>
                    <caption>Collated downloads</caption>
                    <tr><th width=50%>Column Database</th><th width=50%>FITS</th></tr>
                    <tr>
                        <td>
                            Files for each column in "recfile" format.  This is a 
                            <a href="http://code.google.com/p/pycolumns/">columns database</a>.
                            <ul>
                                <li><a href="{serun}.cols">{serun}.cols</a>
                                <li><a href="{serun}.cols/psfstars.cols">{serun}.cols/psfstars.cols</a>
                            </ul>
                        </td>

                        <td>
                            Files for each column in FITS format<b>*</b>
                            <ul>
                                <li><a href="{serun}-fits">{serun}-fits</a> 
                                <li><a href="{serun}-fits/psfstars-fits">{serun}-fits/psfstars-fits</a>
                            </ul>
                        </td>
                    </tr>
                </table>
                {table_comment}
            </table>
            {splitinfo}

		</td>
	</tr>
    <tr>
        <p>
        The collate directories are located on the <a href="https://sites.google.com/site/bnlwlens/software-tutorials/software-setup">astro cluster</a> under<br>
        <pre>
    /astro/u/astrodat/data/DES/wlbnl/{serun}/collated
        </pre>
    </tr>

	<tr>
    {qatext}
	</tr>
    {shear_desc}
</table>



</body>
</html>
    """.format(serun=serun, 
               splitinfo=splitinfo, 
               shear_desc=shear_desc, 
               dataset=rc['dataset'],
               table_comment=table_comment,
               qatext=qatext)

    download_file=path_join(dir,'download.html')
    stdout.write("Writing download file: %s\n" % download_file)
    fobj=open(download_file,'w')
    fobj.write(download)
    fobj.close()







    allcol = """
<html>
<head>
	<link rel="stylesheet" href="table.css" type="text/css">
</head>

<body bgcolor=white>

<h2>Description of Fields in {serun} Shear Table</h2>

Data types are described using this key:

<pre>
f=floating point  
i=integer  
S=string
</pre>

<p>
<table class=simple>
	<tr><th>Column Name</th><th>Data Type<br>[type][bytes]</th><th>Description</th></tr>

	<tr> <td>uid</td>               <td> i4</td>    <td>A unique id</td>    </tr>
	<tr> <td>exposurename</td>      <td> S20</td>   <td>e.g. decam--27--41-i-11</td>  </tr>
	<tr> <td>ccd</td>               <td> i2</td>    <td>ccd number</td>  </tr>
	<tr> <td>id</td>                <td> i2</td>    <td>SExtractor id in this field</td>    </tr>

	<tr> <td>x</td>                 <td> f4</td>    <td>x in CCD</td>    </tr>
	<tr> <td>y</td>                 <td> f4</td>    <td>y in CCD</td> </tr>

	<tr> <td>ra</td>                <td> f8</td>    <td>ra (degrees, J2000)</td>  </tr>
	<tr> <td>dec</td>               <td> f8</td>    <td>dec (degrees, J2000)</td>  </tr>

	<tr> <td>magi</td>              <td> f4</td>    <td>i-band mag, uncalibrated</td>  </tr>

	<tr> <td>size_flags</td>        <td> i2</td>    <td>size error flags</td>  </tr>
	<tr> <td>sigma0</td>            <td> f4</td>    <td>basic size used for <br>star-galaxy separation</td>  </tr>
	<tr> <td>sg</td>                <td> f4</td>    <td>Auxillary star-galaxy separation <br> parameter <b><code>spread_model</code></b></td>  </tr>

	<tr> <td>star_flag</td>         <td> i2</td>    <td>star finder error flags</td>  </tr>

	<tr> <td>shear_flags</td>       <td> i4</td>    <td>shear and shapelets error flags</td>  </tr>
	<tr> <td>shear_s2n</td>         <td> f4</td>    <td>shear S/N for this object. Error on e1/e2 is sqrt(2)/s2n</td> </tr>
	<tr> <td>shapelets_sigma</td>   <td> f4</td>    <td>size used for pre-psf <br>shapelet decompositions </td>  </tr>
	<tr> <td>shapelets_order</td>   <td> i2</td>    <td>Order used for pre-psf <br>shapelet decompositions </td>  </tr>
	<tr> <td>shapelets_prepsf</td>  <td> f4</td>    <td>pre-psf shapelet decompositions </td>  </tr>

	<tr> <td>e1</td>  <td> f4</td>    <td>pre-psf e1 in ra/dec coords</td>  </tr>
	<tr> <td>e2</td>  <td> f4</td>    <td>pre-psf e2 in ra/dec coords</td>  </tr>

	<tr> <td>shear1</td>            <td> f4</td>    <td>shear component 1 in ra/dec coords</td>  </tr>
	<tr> <td>shear2</td>            <td> f4</td>    <td>shear compoment 2 in ra/dec coords</td>  </tr>
	<tr> <td>shear_cov00</td>       <td> f4</td>    <td>err squared for first shear compoment</td>  </tr>
	<tr> <td>shear_cov01</td>       <td> f4</td>    <td>covariance between 1st and 2nd <br>shear components</td>  </tr>
	<tr> <td>shear_cov11</td>       <td> f4</td>    <td>err squared for 2nd shear component</td> </tr>

	<tr> <td>interp_psf_order</td>  <td> f4</td>    <td>order of the PSF shapelets decomposition</td> </tr>
	<tr> <td>interp_psf_sigma</td>  <td> f4</td>    <td>scale of the PSF interpolated to the position of this object</td> </tr>
	<tr> <td>interp_psf_e1</td>     <td> f4</td>    <td>e1 of the PSF interpolated to the position of this object</td> </tr>
	<tr> <td>interp_psf_e2</td>     <td> f4</td>    <td>e2 of the PSF interpolated to the position of this object</td> </tr>
    <tr> <td>interp_psf_shapelets</td>     
                                    <td> f4</td>    <td>PSF shapelets interpolated to the position of this object</td> </tr>

</table>

</body>
</html>
""".format(serun=serun)


    allcol_descr_file=path_join(html_dir, 'all-columns-descr.html')
    stdout.write("Writing allcols description file\n")
    fobj = open(allcol_descr_file,'w')
    fobj.write(allcol)
    fobj.close()



    psfcol="""
<html>
<head>
	<link rel="stylesheet" href="table.css" type="text/css">
</head>

<body bgcolor=white>

<h2>Description of Fields in {serun} PSF Stars Table</h2>

Data types are described using this key:

<pre>
f=floating point  
i=integer  
S=string
</pre>

<p>

<p>
<table class=simple>
	<tr><th>Column Name</th><th>Data Type<br>[type][bytes]</th><th>Description</th></tr>

	<tr> <td>uid</td>               <td> i4</td>    <td>A unique id. Can use to match to the main catalog.</td> </tr>

	<tr> <td>psf_s2n</td>           <td> f4</td>    <td>S/N for this star. Error on e1/e2 is sqrt(2)/s2n </td> </tr>

	<tr> <td>psf_flags</td>         <td> i4</td>    <td>error flags for shapelets measurement</td>    </tr>

	<tr> <td>e1</td>                <td> f4</td>    <td>e1 for this star</td>    </tr>
	<tr> <td>e2</td>                <td> f4</td>    <td>e2 for this star</td>    </tr>

	<tr> <td>shapelets_order</td>   <td> f4</td>    <td>Order for shapelets decomposition</td>    </tr>
	<tr> <td>shapelets_sigma</td>   <td> f4</td>    <td>scale for shapelets decomposition</td>    </tr>
	<tr> <td>shapelets</td>         <td> f4</td>    <td>shapelets decomposition</td>    </tr>

</table>


</body>
</html>
    """.format(serun=serun)


    psfcol_descr_file=path_join(html_dir, 'psfstars-columns-descr.html')
    stdout.write("Writing psfstars cols description file\n")
    fobj = open(psfcol_descr_file,'w')
    fobj.write(psfcol)
    fobj.close()





    qa="""
<html>

	<head>
		<link rel="stylesheet" href="table.css" type="text/css">
	</head>

	<body>

		<h1>QA Plots for single epoch weak lensing run {serun}</h1>

        <h2>Comparison of star shapes and interpolated PSF shapes</h2>
        <table width="500px">
            <tr>
                <td>
                    This plot shows the mean star and PSF ellipticity over all 
                    exposures as a function of position in the focal plane.  Each 
                    ccd has been binned 3x3
                </td>
            </tr>
        </table>
        <p>
        <img src="../plots/psfcheck/allbothbin/all-checkpsf-allbothbin.png">


		<h2>Size-magnitude Diagrams</h2>

		<table width="500px">
			<tr>
				<td>

					This plot shows the size-mag diagram and psf_size/obj_size
					diagram for two different cuts on the data.  The left
					column of plots shows those object that had good size
					measurements.  The right column demands that the
					shear_flags variable is also zero, meaning there were no
					issues whatsoever with the pre-psf shapelets or shear
					determination.  Note this may be too restrictive for
					shear-only studies since the shear measurement is somewhat
					more stable.

					<p>

					Also note that, as expected, the size ratio is better at
					distinguishing the stellar locus than just the FWHM.

					<p>
					These 2-dimensional histograms are shown in an
					<b>asinh</b> stretch to bring out low-level features.
				</td>
			</tr>
		</table>
		<p>
		<img src="../plots/{serun}-{band}-sizemag.png">

		<h2>Testing the consistency of shapelets scales</h2>

		<table width="500px">
			<tr>
				<td>
					This plot shows a comparison between the object scale sigma0
					and the pre-PSF scale sigma.  Also shown is the result
					of "de-convolving" the interpolated PSF, which is
					simply<br>sqrt( sigma0^2 - sigma_PSF^2 )
					<p>
					As expected this is a good approximation of what is
					happening except for objects smaller than the PSF.
				</td>
			</tr>
		</table>
		<p>
		<img src="../plots/{serun}-{band}-sigma0-vs-shapelets-sigma.png">



	</body>
</html>
    """.format(serun=serun,band='i')

    qa_file=path_join(html_dir, 'qa.html')
    stdout.write("Writing qa file\n")
    fobj = open(qa_file,'w')
    fobj.write(qa)
    fobj.close()





    table_css="""
table.simple {
	border-width: 1px;
	border-collapse: collapse;
}

table.simple th {
	background-color: #faf0e6;
	padding: 7px;
	border-width: 1px;
	border-style: inset;
}
table.simple td {
	padding: 7px;
	border-width: 1px;
	border-style: inset;
}

    """

    table_css_file=path_join(html_dir, 'table.css')
    stdout.write("Writing table.css file\n")
    fobj=open(table_css_file,'w')
    fobj.write(table_css)
    fobj.close()



