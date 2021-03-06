from __future__ import print_function
import os
from sys import stderr

import shutil

import esutil as eu
from esutil import numpy_util
from esutil.numpy_util import replicate, where1
from esutil.ostools import path_join, expand_path

import deswl

import numpy
from numpy import arange

try:
    import columns
except:
    print("could not import columns")

def open_columns(run):
    coldir=deswl.files.coldir(run)
    return columns.Columns(coldir)

def make_rid(expid, ccd, fid):
    return expid*10**7 + ccd*10**5 + fid

class SEColumnCollator: 
    def __init__(self, serun, small=False, njob=None, job=None, 
                 no_cleanup=False):
        self.serun = serun
        self.small=small

        self.njob=None if njob is None else int(njob)
        self.job=None if job is None else int(job)

        self.no_cleanup=no_cleanup

        self.set_suffix()
        self.set_names()


    def collate(self):
        try:
            self.load_columns_for_writing()

            self.load_flist()

            ntot = len(self.flist)
            for i,fdict in enumerate(self.flist,1):
                print('-'*70,file=stderr)
                print("Processing %d/%d" % (i,ntot), file=stderr)

                # this will add data to the dictionary
                data=self.load_data(fdict)
                print('writing',file=stderr)
                self.write_data(data)

            self.copy_to_final()
        finally:
            if not self.no_cleanup:
                self.cleanup()

    def write_data(self, data):
        self.write_from_stars(data)
        self.write_from_shear(data)
        self.write_from_psf(data)



    def load_data(self, fdict):
        
        data={}

        if not hasattr(self,'uid_offset'):
            # start id counter
            self.uid_offset = 0

        data['expname'] = fdict['expname']
        data['ccd'] = fdict['ccd']

        for k in ['stars','psf','shear']:
            kk=k
            data[k] = eu.io.read(fdict[kk], ensure_native=True, verbose=True)

        data['idvals'] = self.uid_offset + numpy.arange(data['stars'].size,
                                                        dtype='i4')
        self.uid_offset += data['stars'].size

        return data

    def set_names(self):
        self.temp_coldir  = self.coldir(temp=True)
        self.final_coldir = self.coldir()

    def set_suffix(self):
        if self.job is not None:
            if self.njob is None:
                raise ValueError("Send both job and njob")
            if self.job < 0 or self.job > (self.njob-1):
                raise ValueError("job must be in [0,%s], "
                                 "got %s" % (self.njob-1,self.job))
            self.suffix = '-%03d' % self.job
        else:
            self.suffix = ''



    def coldir(self, fits=False, temp=False):
        coldir = deswl.files.coldir(self.serun, fits=fits, suffix=self.suffix)
        coldir = expand_path(coldir)

        if temp:
            if fits:
                raise ValueError("don't use temp=True for fits")
            coldir=os.path.basename(coldir)
            coldir=os.path.join('/data/esheldon/tmp',coldir)
        return coldir


    def load_columns_for_writing(self):

        print('Will write to temporary coldir:',self.temp_coldir,file=stderr)
        print('Will copy to:',self.final_coldir,file=stderr)

        if os.path.exists(self.temp_coldir):
            raise RuntimeError("temp coldir %s exists. Please start "
                               "from scratch" % self.temp_coldir)
        if os.path.exists(self.final_coldir):
            raise RuntimeError("Final coldir %s exists. Please start "
                               "from scratch" % self.final_coldir)

        self.cols = columns.Columns(self.temp_coldir)
        self.cols.create()

        # load the sub columns dir for PSF stars
        psf_subdir = os.path.join(self.temp_coldir, 'psfstars.cols')
        self.cols.load_coldir(psf_subdir)

    def copy_to_final(self):
        print("Copying to final dest:",self.final_coldir,file=stderr)
        shutil.copytree(self.temp_coldir, self.final_coldir)

    def cleanup(self):
        print("Cleaning up temp dir:",self.temp_coldir,file=stderr)
        if os.path.exists(self.temp_coldir):
            shutil.rmtree(self.temp_coldir)

    def load_flist(self):
        flistfile=deswl.files.collated_path(self.serun,'goodlist')
        flist=eu.io.read(flistfile, verbose=True)

        if self.job is not None:

            nf=len(flist)
            nper=nf/self.njob

            if self.job == (self.njob-1):
                flist=flist[self.job*nper:]
            else:
                flist=flist[self.job*nper:(self.job+1)*nper]

        self.flist=flist

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


        shapelets_prepsf=numpy.array(data['shear']['shapelets_prepsf'], 
                                     dtype='f4')
        if not self.small:
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
        if not self.small:
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
        exp_array[:] = data['expname']
        cols.write_column('expname',exp_array)

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


        pind, sind = eu.numpy_util.match(psf['id'], stars['id'])

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

        print('coldir:',coldir,file=stderr)
        print('fitsdir:',fitsdir,file=stderr)
        print('psfstars_fitsdir:',psfstars_fitsdir,file=stderr)

        cols = columns.Columns(coldir)

        print(cols,file=stderr)

        for col in sorted(cols):
            if col in ['shapelets_prepsf','interp_psf_shapelets']:
                print("skipping large column:",col,file=stderr)
            else:
                print('column: ',col,file=stderr)

                if col == 'psfstars':
                    continue

                fname = cols[col].filename

                data = eu.sfile.read(fname)

                fitsfile = os.path.join(fitsdir, col+'.fits')
                print(fitsfile,file=stderr)

                # don't make a copy!
                eu.io.write(fitsfile, data, clobber=True,copy=False)

                del data

        psfstars_cols = cols['psfstars']
        for col in sorted(psfstars_cols):
            print('column: ',col,file=stderr)

            fname = psfstars_cols[col].filename

            data = eu.sfile.read(fname)

            fitsfile = os.path.join(psfstars_fitsdir, col+'.fits')
            print(fitsfile,file=stderr)

            # don't make a copy!
            eu.io.write(fitsfile, data, clobber=True,copy=False)

            del data

    def create_indexes(self):

        coldir=self.coldir()
        cols = columns.Columns(coldir,verbose=True)
        if not cols.dir_exists():
            raise RuntimeError("You haven't added the data yet")

        # create some indexes
        # after this, data automatically gets added to the index
        cols['ccd'].create_index()
        cols['size_flags'].create_index()
        cols['star_flag'].create_index()
        cols['shear_flags'].create_index()
        cols['expname'].create_index()

        cols['psfstars']['psf_flags'].create_index()
        cols['psfstars']['uid'].create_index()

    def write_html(self):
        write_se_collate_html(self.serun)


    def add_rid(self):
        import desdb

        print("getting all expnames",file=stderr)
        rc=deswl.files.Runconfig(self.serun)
        all_expnames = desdb.files.get_expnames(rc['dataset'],rc['band'])
        all_expnames = numpy.array(all_expnames)
        all_expnames.sort()

        print("reading data",file=stderr)
        cols=open_columns(self.serun)
        expnames    = cols['expname'][:]
        ccds        = cols['ccd'][:]
        ids         = cols['id'][:]

        # we want this because we need to search the expname list to get and id
        # and we don't want to do it too many times
        print("getting expnames argsort",file=stderr)
        s=expnames.argsort()
        rid = numpy.zeros(ids.size, dtype='i8')

        print("getting rids...",file=stderr)
        expold=None
        for i in xrange(ids.size):
            if (i % 10000) == 0:
                print('%d/%d' % (i, ids.size),file=stderr)
            si=s[i] 

            expname = expnames[si]
            ccd = ccds[si]
            id = ids[si]

            if expname != expold:
                expid=where1(all_expnames == expname)
                expold = expname

                if expid.size == 0:
                    raise ValueError("expname not found: '%s'" % expname)
                expid=expid[0]

            rid[si] = make_rid(expid, ccd, id)

        print("\nwriting rid",file=stderr)
        cols.write_column('rid', rid, create=True)

class CollateWQJob:
    def __init__(self, run, njob, hosts=None, priority='low'):
        self.run=run
        self.njob=njob
        self.priority=priority

        if hosts is None:
            hosts = ['astro%04i' % i for i in xrange(1,12)]
            hosts += ['astro%04i' % i for i in xrange(21,34)]

        self.hosts=hosts

    def job_template(self):
        text = """
command: |
    hostname
    source ~esheldon/.bashrc
    source /opt/astro/SL53/bin/setup.hadoop.sh
    source ~astrodat/setup/setup.sh
    source ~/.dotfiles/bash/astro.bnl.gov/modules.sh
    {esutil_load}
    {tmv_load}
    {wl_load}
    module unload espy && module load espy/work
    source ~esheldon/local/des-oracle/setup.sh

    script=$ESPY_DIR/des/bin/collate2columns.py
    python $script -s -n {njob} -j {job} {run}

mode: byhost
host: {host}
priority: {priority}
job_name: {job_name}\n""" 
        return text

    def write(self):

        rc=deswl.files.Runconfig(self.run)
        wl_load = deswl.files._make_load_command('wl',rc['wlvers'])

        tmv_load = ''
        if self.run[0:2] == 'se' or self.run[0:2] == 'me':
            tmv_load = deswl.files._make_load_command('tmv',rc['tmvvers'])
        esutil_load = deswl.files._make_load_command('esutil', rc['esutilvers'])

        job_template=self.job_template()
        outd=deswl.files.wq_dir(self.run, subdir='collate')
        if not os.path.exists(outd):
            os.makedirs(outd)

        nh=len(self.hosts)
        for job in xrange(self.njob):
            job_name='%s-collate-%03i' % (self.run,job)
            job_file=os.path.join(outd,job_name+'.yaml')

            host = self.hosts[job % nh]
            text=job_template.format(esutil_load=esutil_load,
                                     tmv_load=tmv_load,
                                     wl_load=wl_load,
                                     host=host,
                                     njob=self.njob,
                                     job=job,
                                     run=self.run,
                                     job_name=job_name,
                                     priority=self.priority)

            print("Writing job file:",job_file,file=stderr)
            with open(job_file,'w') as fobj:
                fobj.write(text)

        comb_name=os.path.join(outd,'%s-combine.py' % self.run)
        print("Writing combine script:",comb_name,file=stderr)

        collated_dir=deswl.files.collated_dir(self.run)
        with open(comb_name,'w') as fobj:
            text="""
import columns
import glob

coldir='%(dir)s/%(run)s.cols'

pattern='%(dir)s/%(run)s-*.cols'
f=glob.glob(pattern)
f.sort()

c=columns.Columns(coldir)
c.from_columns(f,create=True)\n""" % {'dir':collated_dir,
                                      'run':self.run}
            fobj.write(text)



class MEColumnCollator: 
    def __init__(self, run, small=False, no_cleanup=False, njob=None, job=None):
        self.run = run
        self.small=small

        self.njob=None if njob is None else int(njob)
        self.job=None if job is None else int(job)

        self.no_cleanup=no_cleanup

        self.set_suffix()
        self.set_names()

    def collate(self):
        try:
            self.load_columns_for_writing()

            self.load_flist()

            ntot = len(self.flist)
            uid0=0
            cat_cols=['mag_model','magerr_model','x_image','y_image',
                      'flags_weight']
            for i,fdict in enumerate(self.flist,1):
                print('-'*70,file=stderr)
                print("Processing %d/%d" % (i,ntot),file=stderr)

                # eu.io.read works correctly with hdfs
                fname=fdict['multishear']
                print("    ",fname,file=stderr)
                data = eu.io.read(fname)

                catname=fdict['cat']
                print("    ",catname,file=stderr)
                cat=eu.io.read(catname,columns=cat_cols,lower=True)

                if cat.size != data.size:
                    raise ValueError("cat and multishear sizes don't "
                                     "match: %d/%d" % (cat.size,data.size))

                uids = uid0 + arange(data.size,dtype='i4')
                tilenames = replicate(str(fdict['tilename']), data.size)

                self.write(data)
                self.write(cat)
                self.cols.write_column('tilename', tilenames)
                self.cols.write_column('uid', uids)

                uid0 += data.size
            self.copy_to_final()
        finally:
            if not self.no_cleanup:
                self.cleanup()

    def write(self,data):
        for c in data.dtype.names:
            if self.small and c == 'shapelets_prepsf':
                continue

            if c in ['nimages_found','nimages_gotpix','gal_order']:
                d=numpy.array(data[c], dtype='i1')
            elif c == 'input_flags':
                d=numpy.array(data[c], dtype='i2')
            elif c == 'shear_signal_to_noise':
                d=data[c]
                c = 'shear_s2n'
            elif c == 'shapelet_sigma':
                d=data[c]
                c = 'shapelets_sigma'
            elif c == 'shear_cov00':
                d = data[c]
                c = 'shear_cov11'
            elif c == 'shear_cov01':
                d = data[c]
                c = 'shear_cov12'
            elif c == 'shear_cov11':
                d = data[c]
                c = 'shear_cov22'
            else:
                d = data[c]

            numpy_util.to_native(d, inplace=True)
            self.cols.write_column(c, d)




    def load_flist(self):
        flistfile=deswl.files.collated_path(self.run,'goodlist')
        flist=eu.io.read(flistfile, verbose=True)

        if self.job is not None:

            nf=len(flist)
            nper=nf/self.njob

            if self.job == (self.njob-1):
                flist=flist[self.job*nper:]
            else:
                flist=flist[self.job*nper:(self.job+1)*nper]

        self.flist=flist



    def load_columns_for_writing(self):

        print('Will write to temporary coldir:',self.temp_coldir,file=stderr)
        print('Will copy to:',self.final_coldir,file=stderr)

        if os.path.exists(self.temp_coldir):
            raise RuntimeError("temp coldir %s exists. Please start "
                               "from scratch" % self.temp_coldir)
        if os.path.exists(self.final_coldir):
            raise RuntimeError("Final coldir %s exists. Please start "
                               "from scratch" % self.final_coldir)

        self.cols = columns.Columns(self.temp_coldir)
        self.cols.create()

    def copy_to_final(self):
        print("Copying to final dest:",self.final_coldir,file=stderr)
        shutil.copytree(self.temp_coldir, self.final_coldir)

    def cleanup(self):
        print("Cleaning up temp dir:",self.temp_coldir,file=stderr)
        if os.path.exists(self.temp_coldir):
            shutil.rmtree(self.temp_coldir)



    def convert2fits(self):
        """
        Just make fits versions of all the files
        """

        coldir=deswl.files.coldir(self.run)
        fitsdir = coldir.replace('.cols','-fits')

        if not os.path.exists(fitsdir):
            os.makedirs(fitsdir)

        print('coldir:',coldir,file=stderr)
        print('fitsdir:',fitsdir,file=stderr)

        cols = columns.Columns(coldir)

        print(cols,file=stderr)

        for col in sorted(cols):
            if col in ['shapelets_prepsf']:
                print("skipping large column:",col,file=stderr)
            else:
                print('column: ',col,file=stderr)

                fname = cols[col].filename

                # use sfile so we get a rec array
                data = eu.sfile.read(fname)

                fitsfile = os.path.join(fitsdir, col+'.fits')
                print("    ",fitsfile,file=stderr)

                # don't make a copy!
                eu.io.write(fitsfile, data, clobber=True)

                del data

    def set_names(self):
        self.temp_coldir  = self.coldir(temp=True)
        self.final_coldir = self.coldir()

    def coldir(self, fits=False, temp=False):
        coldir = deswl.files.coldir(self.run, fits=fits, suffix=self.suffix)
        coldir = expand_path(coldir)

        if temp:
            if fits:
                raise ValueError("don't use temp=True for fits")
            coldir=os.path.basename(coldir)
            coldir=os.path.join('/data/wlpipe',coldir)
        return coldir


    def set_suffix(self):
        if self.job is not None:
            if self.njob is None:
                raise ValueError("Send both job and njob")
            if self.job < 0 or self.job > (self.njob-1):
                raise ValueError("job must be in [0,%s], "
                                 "got %s" % (self.njob-1,self.job))
            self.suffix = '-%03d' % self.job
        else:
            self.suffix = ''


    def create_indexes(self):

        coldir=deswl.files.coldir(self.run)
        cols = columns.Columns(coldir,verbose=True)
        if not cols.dir_exists():
            raise RuntimeError("You haven't added the data yet")

        # create some indexes
        # after this, data automatically gets added to the index
        for c in ['id','input_flags','shear_flags','flags_weight','mag_model',
                  'shear_s2n','ra','dec','nimages_found','nimages_gotpix']:
            if cols[c].has_index():
                print("skipping '%s' which already has an index" % c,
                      file=stderr)
            else:
                print("creating index for column:",c,file=stderr)
                cols[c].create_index()


    def write_html(self):

        rc = deswl.files.Runconfig(self.run)

        collate_dir = deswl.files.collated_dir(self.run)
        html_dir = path_join(collate_dir, 'html')
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)

        coldir = deswl.files.coldir(self.run)
        fitsdir = coldir.replace('.cols','-fits')
        plotdir=path_join(collate_dir,'plots')

        if os.path.exists(fitsdir):
            table_comment=('<i>* The collated FITS files do not '
                           'include the large '
                           'column <code>shapelets_prepsf</code></i>')
        else:
            table_comment='<i><b>* FITS Not Yet Available</b></i>'

        if os.path.exists(plotdir):
            qatext = """
            <td>
                <p>
                <p>
                <h2>QA and Analysis</h2>
                Here are  some <a href="html/qa.html">QA plots and tests</a>.
                <p>
                <p>
            </td>
            """
        else:
            qatext=""

        table_css_file=path_join(html_dir, 'table.css')
        write_css(table_css_file)

        download="""
<html>
<head>
	<link rel="stylesheet" href="html/table.css" type="text/css">
</head>

<body bgcolor=white>

<h1>{dataset} Multi-Epoch Catalogs: Run {run}</h1>

<table width="500px">
	<tr>
		<td>

            <h2>Raw WL outputs</h2>
            Raw outputs can be found <a href="..">Here</a>.   If you plan to download
            <i><b>all the data</b></i>, please contact Erin Sheldon.
            <p>
            These files are located on the <a href="https://sites.google.com/site/bnlwlens/software-tutorials/software-setup">astro cluster</a> at
            <pre>
    /astro/u/astrodat/data/DES/wlbnl/{run}
            </pre>
			<h2>Collated WL outputs</h2>
            <ul>
                <li><a href="html/columns-descr.html">Description of the table</a>
            </ul>

            This "table" is actually just a file for each column.

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
                                <li><a href="{run}.cols">{run}.cols</a>
                            </ul>
                        </td>

                        <td>
                            Files for each column in FITS format<b>*</b>
                            <ul>
                                <li><a href="{run}-fits">{run}-fits</a> 
                            </ul>
                        </td>
                    </tr>
                </table>
                {table_comment}
            </table>

		</td>
	</tr>
    <tr>
        <p>
        The collate directories are located on the <a href="https://sites.google.com/site/bnlwlens/software-tutorials/software-setup">astro cluster</a> under<br>
        <pre>
    /astro/u/astrodat/data/DES/wlbnl/{run}/collated
        </pre>
    </tr>

	<tr>
    {qatext}
	</tr>
</table>



</body>
</html>
        """.format(run=self.run, 
                   dataset=rc['dataset'],
                   table_comment=table_comment,
                   qatext=qatext)

        download_file=path_join(collate_dir,'download.html')
        print("Writing download file:",download_file,file=stderr)
        with open(download_file,'w') as fobj:
            fobj.write(download)


        coldescr = """
<html>
<head>
	<link rel="stylesheet" href="table.css" type="text/css">
</head>

<body bgcolor=white>

<h2>Description of Fields in {run} Multi-epoch Shear Table</h2>

Data types are described using this key:

<pre>
f  -> floating point  
i  -> integer  
S  -> string
[] -> vector
</pre>

<p>
<table class=simple>
	<tr><th>Column Name</th><th>Data Type<br>[type][bytes]</th><th>Description</th></tr>

	<tr> <td>uid</td>               <td> i4</td>    <td>A unique id for this run</td>    </tr>
	<tr> <td>tilename</td>          <td> S12</td>   <td>Coadd tilename</td>    </tr>
	<tr> <td>id</td>                <td> i4</td>    <td>SExtractor id in this field</td>    </tr>
	<tr> <td>x_image</td>           <td> f4</td>    <td>SExtractor x position in image (in python im[y,x])</td></tr>
	<tr> <td>y_image</td>           <td> f4</td>    <td>SExtractor y position in image (in python im[y,x])</td></tr>

	<tr> <td>ra</td>                <td> f8</td>    <td>SExtractor ra (degrees, J2000)</td>  </tr>
	<tr> <td>dec</td>               <td> f8</td>    <td>SExtractor dec (degrees, J2000)</td>  </tr>

	<tr> <td>input_flags</td>       <td> i2</td>    <td>SExtractor catalog flags</tr>
	<tr> <td>flags_weight</td>      <td> i2</td>    <td>SExtractor coadd weights flags</tr>

	<tr> <td>mag_model</td>         <td> f4</td>    <td>SExtractor i-band </td>  </tr>
	<tr> <td>magerr_model</td>      <td> f4</td>    <td>SExtractor i-band </td>  </tr>

	<tr> <td>nimages_found</td>     <td> i1</td>    <td>Number of images found to overlap</tr>
	<tr> <td>nimages_gotpix</td>    <td> i1</td>    <td>Number of images that contributed pixels</tr>

	<tr> <td>shear_flags</td>       <td> i4</td>    <td>shear and shapelets error flags</td>  </tr>
	<tr> <td>shear_s2n</td>         <td> f8</td>    <td>shear S/N for this object. Error on e1/e2 is sqrt(2)/s2n</td> </tr>
	<tr> <td>shapelets_sigma</td>   <td> f8</td>    <td>size used for pre-psf <br>shapelet decompositions </td>  </tr>
	<tr> <td>shapelets_prepsf</td>  <td> f8[]</td>    <td>pre-psf shapelet decompositions </td>  </tr>
	<tr> <td>gal_order</td>         <td> i1</td>    <td>Order used for shear </td>  </tr>

	<tr> <td>shear1</td>            <td> f8</td>    <td>shear component 1 in ra/dec coords</td>  </tr>
	<tr> <td>shear2</td>            <td> f8</td>    <td>shear compoment 2 in ra/dec coords</td>  </tr>
	<tr> <td>shear_cov00</td>       <td> f8</td>    <td>err squared for first shear compoment</td>  </tr>
	<tr> <td>shear_cov01</td>       <td> f8</td>    <td>covariance between 1st and 2nd <br>shear components</td>  </tr>
	<tr> <td>shear_cov11</td>       <td> f8</td>    <td>err squared for 2nd shear component</td> </tr>

</table>

</body>
</html>
""".format(run=self.run)


        descr_file=path_join(html_dir, 'columns-descr.html')
        print("Writing allcols description file",descr_file,file=stderr)
        with open(descr_file,'w') as fobj:
            fobj.write(coldescr)



def write_se_collate_html(serun):
    """

    Need to make collation band-aware!

    """

    rc = deswl.files.Runconfig(serun)

    dir=deswl.files.collated_dir(serun)
    html_dir = path_join(dir, 'html')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)

    fitsdir=path_join(dir,serun+'-fits')
    if os.path.exists(fitsdir):
        table_comment=('<i>Note: The collated FITS files do not '
                       'include the large columns '
                       '<code>shapelets_prepsf</code> and <code>'
                       'interp_psf_shapelets</code></i>')
    else:
        table_comment='<i><b>FITS Not Yet Available</b></i>'


    collate_dir= deswl.files.collated_dir(serun)
    plotdir=path_join(collate_dir,'plots')
    szmagdir=path_join(collate_dir,'../test/checksg/byexp')

    band='i'
    tfile= path_join(plotdir, '%s-%s-sizemag.png' % (serun,band))
    if os.path.exists(tfile) or os.path.exists(szmagdir):
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
               shear_desc=shear_desc, 
               dataset=rc['dataset'],
               table_comment=table_comment,
               qatext=qatext)

    download_file=path_join(dir,'download.html')
    print("Writing download file:",download_file,file=stderr)
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
f  -> floating point  
i  -> integer  
S  -> string
[] -> vector
</pre>

<p>
<table class=simple>
	<tr><th>Column Name</th><th>Data Type<br>[type][bytes]</th><th>Description</th></tr>

	<tr> <td>uid</td>               <td> i4</td>    <td>A unique id for this run</td>    </tr>
	<tr> <td>expname</td>      <td> S20</td>   <td>e.g. decam--27--41-i-11</td>  </tr>
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
	<tr> <td>shapelets_prepsf</td>  <td> f4[]</td>    <td>pre-psf shapelet decompositions </td>  </tr>

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
    print("Writing allcols description file",file=stderr)
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
    print("Writing psfstars cols description file",file=stderr)
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
    """.format(serun=serun)

    if os.path.exists(szmagdir):
        qa += """
    <h2>Size-mag diagrams for each exposure</h2>

            Size-mag and spreadmodel-mag plots for each exposure are <a
            href="../../test/checksg/byexp/phpshow.php">here</a>.  ccds are shown in different
            colors.
        """
    if os.path.exists(plotdir):
        qa += """
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
    print("Writing qa file",file=stderr)
    fobj = open(qa_file,'w')
    fobj.write(qa)
    fobj.close()

    table_css_file=path_join(html_dir, 'table.css')

    write_css(table_css_file)

def write_css(fname):
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

    print("Writing table.css file",fname,file=stderr)
    fobj=open(fname,'w')
    fobj.write(table_css)
    fobj.close()


class AMColumnCollator: 
    def __init__(self, run, njob=None, job=None):
        self.run = run

        self.njob=None if njob is None else int(njob)
        self.job=None if job is None else int(job)

        self.ftypes = ['am']
        self.set_suffix()
        self.set_names()
        self.load_expnames()

    def collate(self):
        try:
            self.load_columns_for_writing()

            self.load_flist()

            ntot = len(self.flist)
            for i,fdict in enumerate(self.flist,1):
                print('-'*70,file=stderr)
                print("Processing %d/%d" % (i,ntot), file=stderr)

                # this will add data to the dictionary
                data=self.load_data(fdict)
                print('writing',file=stderr)
                self.write_data(data)

            self.copy_to_final()
        finally:
            self.cleanup()


    def load_expnames(self):
        if not hasattr(self,'expnames'):
            import desdb
            rc=deswl.files.Runconfig(self.run)
            expnames = desdb.files.get_expnames(rc['dataset'],rc['band'])
            expnames = numpy.array(expnames)
            expnames.sort()
            self.expnames=expnames

    def load_data(self, fdict):
        
        data={}

        #if not hasattr(self,'uid_offset'):
        #    # start id counter
        #    self.uid_offset = 0
    
        expname = fdict['expname']
        ccd     = fdict['ccd']
        data['expname'] = expname
        data['ccd']     = ccd

        data['am'] = eu.io.read(fdict['output_files']['am'], 
                                ensure_native=True, 
                                verbose=True)

        #data['idvals'] = \
        #        self.uid_offset + numpy.arange(data['am'].size,dtype='i4')
        
        ids = data['am']['number']
        data['rid'] = self.make_rids(expname, ccd, ids)

        return data

    def make_rids(self, expname, ccd, ids):

        rid = numpy.zeros(ids.size, dtype='i8')
        expid=where1(self.expnames == expname)
        if expid.size == 0:
            raise ValueError("expname not found: '%s'" % expname)
        for i in xrange(ids.size):
            rid[i] = make_rid(expid[0], ccd, ids[i])

        return rid

    def set_names(self):
        self.temp_coldir  = self.coldir(temp=True)
        self.final_coldir = self.coldir()

    def set_suffix(self):
        if self.job is not None:
            if self.njob is None:
                raise ValueError("Send both job and njob")
            if self.job < 0 or self.job > (self.njob-1):
                raise ValueError("job must be in [0,%s], "
                                 "got %s" % (self.njob-1,self.job))
            self.suffix = '-%03d' % self.job
        else:
            self.suffix = ''



    def coldir(self, fits=False, temp=False):
        coldir = deswl.files.coldir(self.run, fits=fits, suffix=self.suffix)
        coldir = expand_path(coldir)

        if temp:
            if fits:
                raise ValueError("don't use temp=True for fits")
            coldir=os.path.basename(coldir)
            coldir=os.path.join('/data/esheldon/tmp',coldir)
        return coldir


    def load_columns_for_writing(self):

        print('Will write to temporary coldir:',self.temp_coldir,file=stderr)
        print('Will copy to:',self.final_coldir,file=stderr)

        if os.path.exists(self.temp_coldir):
            raise RuntimeError("temp coldir %s exists. Please start "
                               "from scratch" % self.temp_coldir)
        if os.path.exists(self.final_coldir):
            raise RuntimeError("Final coldir %s exists. Please start "
                               "from scratch" % self.final_coldir)

        self.cols = columns.Columns(self.temp_coldir)
        self.cols.create()

    def copy_to_final(self):
        print("Copying to final dest:",self.final_coldir,file=stderr)
        shutil.copytree(self.temp_coldir, self.final_coldir)

    def cleanup(self):
        print("Cleaning up temp dir:",self.temp_coldir,file=stderr)
        if os.path.exists(self.temp_coldir):
            shutil.rmtree(self.temp_coldir)

    def load_flist(self):
        flistfile=deswl.files.collated_path(self.run,'goodlist')
        flist=eu.io.read(flistfile, verbose=True)

        if self.job is not None:

            nf=len(flist)
            nper=nf/self.njob

            if self.job == (self.njob-1):
                flist=flist[self.job*nper:]
            else:
                flist=flist[self.job*nper:(self.job+1)*nper]

        self.flist=flist

    def write_data(self, data):
        cols=self.cols

        # make copies of the data before writing to the columns.  This is
        # to force the data type, including native endianness


        f8cols = ['row','col','Irr','Irc','Icc','e1','e2',
                  'rho4','a4','s2','uncer','s2n','wrow','wcol',
                  'sky','sigsky']
        for col in f8cols:
            vals=numpy.array(data['am'][col], dtype='f8')
            cols.write_column(col,vals)

        i4cols = ['numiter','whyflag']
        for col in i4cols:
            vals=numpy.array(data['am'][col], dtype='i4')
            cols.write_column(col,vals)

        # our id column
        #cols.write_column('uid', data['idvals'])
        cols.write_column('rid', data['rid'])

        # Same for all objects in this set
        ccd_array = numpy.zeros(data['am'].size, dtype='i1')
        ccd_array[:] = data['ccd']
        cols.write_column('ccd',ccd_array)

        exp_array = numpy.zeros(data['am'].size, dtype='S20')
        exp_array[:] = data['expname']
        cols.write_column('expname',exp_array)


        del ccd_array
        del exp_array


        # from sextractor
        id = numpy.array(data['am']['number'], dtype='i2')
        cols.write_column('id',id)

        mag = numpy.array(data['am']['mag_model'], dtype='f4')
        cols.write_column('mag_model',mag)

        flags = numpy.array(data['am']['flags'], dtype='i2')
        cols.write_column('sxflags',flags)

        ra = numpy.array(data['am']['alphawin_j2000'], dtype='f8')
        cols.write_column('ra',ra)

        dec = numpy.array(data['am']['deltawin_j2000'], dtype='f8')
        cols.write_column('dec',dec)

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

        print('coldir:',coldir,file=stderr)
        print('fitsdir:',fitsdir,file=stderr)
        print('psfstars_fitsdir:',psfstars_fitsdir,file=stderr)

        cols = columns.Columns(coldir)

        print(cols,file=stderr)

        for col in sorted(cols):
            if col in ['shapelets_prepsf','interp_psf_shapelets']:
                print("skipping large column:",col,file=stderr)
            else:
                print('column: ',col,file=stderr)

                if col == 'psfstars':
                    continue

                fname = cols[col].filename

                data = eu.sfile.read(fname)

                fitsfile = os.path.join(fitsdir, col+'.fits')
                print(fitsfile,file=stderr)

                # don't make a copy!
                eu.io.write(fitsfile, data, clobber=True,copy=False)

                del data

        psfstars_cols = cols['psfstars']
        for col in sorted(psfstars_cols):
            print('column: ',col,file=stderr)

            fname = psfstars_cols[col].filename

            data = eu.sfile.read(fname)

            fitsfile = os.path.join(psfstars_fitsdir, col+'.fits')
            print(fitsfile,file=stderr)

            # don't make a copy!
            eu.io.write(fitsfile, data, clobber=True,copy=False)

            del data

    def create_indexes(self):

        coldir=self.coldir()
        cols = columns.Columns(coldir,verbose=True)
        if not cols.dir_exists():
            raise RuntimeError("You haven't added the data yet")

        # create some indexes
        # after this, data automatically gets added to the index
        cols['ccd'].create_index()
        cols['size_flags'].create_index()
        cols['star_flag'].create_index()
        cols['shear_flags'].create_index()
        cols['expname'].create_index()

        cols['psfstars']['psf_flags'].create_index()
        #cols['psfstars']['uid'].create_index()

    def write_html(self):
        write_se_collate_html(self.run)



