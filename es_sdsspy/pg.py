import os,sys
from sys import stdout,stderr
import time

try:
    import pgnumpy
except:
    pass


import esutil
import sdsspy

import numpy
from numpy import where

class WindowStuffer:
    def __init__(self):
        self._window = sdsspy.window.Window()
        self._pg=pgnumpy.connect('dbname=boss user=postgres')

    def _index_colspec(self, type):
        if type == 'flist':
            # note multi-column index
            return ['run,camcol,field','rerun','ra','dec','score']

    def stuff(self, type=None):
        """
        Stuff all the types into database tables.  We will proceed the
        name with window_, e.g. window_flist
        """

        # if no types are sent, run through them all and
        # call recursively
        if type is None:
            for type in self._window._types:
                self.stuff(type)
            return


        if type == 'flist':
            t=self._window.read(type)
            tdata = t[type]
            data = esutil.numpy_util.to_native(tdata)
            del tdata
        else:
            # nothing else implemented yet
            return
        
        table = self.type_table(type)
        pgnumpy.array2table(data, table, 
                            tmpdir='~/tmp',
                            conninfo='dbname=boss user=postgres')
        self.grant(type)
        self.add_indexes(type)


    def add_indexes(self, type):
        table = self.type_table(type)
        idx_colspecs = self._index_colspec(type)
        for colspec in idx_colspecs:
            self._pg.create_index(table, colspec)

    def grant(self, type):
        table = self.type_table(type)
        grant_query = 'grant select on %s to public' % table
        print grant_query
        self._pg.execute(grant_query)

    def type_table(self, type):
        table = 'window_%s' % type
        return table

class SweepStuffer:
    def __init__(self, type, update=False, verbose=False, tmpdir='.'):
        """
        Class:
            SweepStuffer
        Purpose:
            Stuff the datasweeps into a table
        Construction:
            ss = SweepStuffer(type, update=False)

                if update=True, then the table can already exist and we are
                just adding to it.  Be careful with that!

        Examples:
            ss = SweepStuffer()

            # stuff all runs
            ss.stuff_runs(runs=ALL)

            # stuff a single run/camcol
            ss.stuff_one(run,rerun,camcol)

        """

        self._verbose=verbose

        self._type = type
        self._table = self._type
        self._stuffall=all
        self._update = update

        self._conninfo='dbname=boss user=postgres'

        self._htm_depth = 10
           
        #self._primary_key = 'photoid'
        self._index_columns = ['ra',
                               'dec',
                               'run,camcol,field,id',  # note combined index
                               'photoid',
                               'isprimary',
                               'modelmag_dered[1]',
                               'modelmag_dered[2]',
                               'modelmag_dered[3]',
                               'modelmag_dered[4]',
                               'modelmag_dered[5]',
                               'htmid10',
                               'basic_maskflags',
                               'tycho_maskflags']

        self._tmpdir = tmpdir

        # for the resolve
        self._minscore = 0.1

        self._photo_sweep = os.environ['PHOTO_SWEEP']
        self._photo_resolve = os.environ['PHOTO_RESOLVE']

        stdout.write('type: %s\n' % self._type)
        stdout.write('PHOTO_SWEEP: %s\n' % self._photo_sweep)
        stdout.write('PHOTO_RESOLVE: %s\n' % self._photo_resolve)



    def stuff_runs(self, runs=None):
        """
        Loop over runs and stuff them into the database table.  If runs
        is not entered, they are gotten using the window runlist
        """

        tm00 = time.time()
        # if not updating, make sure the table does not already exist
        if not self._update:
            pg = pgnumpy.connect(self._conninfo)
            if pg.table_exists(self._table):
                raise ValueError("Table '%s' already exists" % self._table)

        runs2use,reruns2use = self.matchruns(runs=runs)

        nrun=len(runs2use)

        srun = runs2use.argsort()
        irun=1
        for i in srun:
            run = runs2use[i]
            rerun = reruns2use[i]

            stdout.write("run: %06i (%s/%s)\n" % (run,irun,nrun))

            tm0 = time.time()
            for camcol in [1,2,3,4,5,6]:
                self.stuff_one(run,rerun,camcol)
            tm1=time.time()

            stdout.write("  Time: %s sec\n" % (tm1-tm0,) )
            stdout.flush()

            # grant select privledges to everyone
            if irun == 1:
                grant_query = 'grant select on %s to public' % self._table
                stdout.write('    '+grant_query+'\n')
                pg.execute(grant_query)
            irun += 1

        tm11 = time.time()
        stdout.write("Total time: ")
        esutil.misc.ptime(tm11-tm00)

    def stuff_one(self, run, rerun, camcol):
        """
        Stuff a single run,rerun,camcol into the table
        """
        stdout.write("  %06i %s %s:  " % (run,rerun,camcol))
        stdout.flush()

        data = self.read_sweep(run,rerun,camcol)

        output = self.get_subset_and_derived(data)

        self.write_to_table(output)

        del output
        del data

    def read_sweep(self, run, rerun, camcol):
        """
        Read a calibObj file with the given id info
        """
        # we ensure native so that the ascii write will work correctly
        data = sdsspy.files.read('calibobj',run,rerun,camcol=camcol,
                                 type=self._type, ensure_native=True)
        return data

    def write_to_table(self,output):
        """
        Use our PgInput object to write data to the table
        """
        # create the pginput object if it does not yet exist
        self.create_pginput(output)

        self._pginput.write(output)
        stdout.flush()

    def create_pginput(self, array):
        """
        Create the PgInput object if it does not exist
        """
        if not hasattr(self, 'pginput'):
            self._pginput = pgnumpy.PgInput(array.dtype, self._table, 
                                            conninfo=self._conninfo,
                                            tmpdir=self._tmpdir,
                                            verbose=self._verbose)


    def create_indices(self):
        """
        Create all the indices in the _index_columns variable
        """
        tm00 = time.time()
        pg = pgnumpy.connect(self._conninfo)
        for colspec in self._index_columns:

            stdout.write('-'*70 + '\n')

            if colspec == 'photoid':
                unique=True
            else:
                unique=False

            tm0 = time.time()
            pg.create_index(self._table, colspec, unique=unique)
            tm1 = time.time()
            esutil.misc.ptime(tm1-tm0)

        tm11 = time.time()
        stdout.write('-'*70 + '\n')
        stdout.write("Total time: ")
        esutil.misc.ptime(tm11-tm00)

    def create_meta_tables(self, subtypes=['','prim']):
        """
        Before calling this you must create the primary view

        Note count will be a bigint
        """

        photo_sweep = os.path.basename(self._photo_sweep)
        photo_resolve = os.path.basename(self._photo_resolve)


        pg = pgnumpy.connect(self._conninfo)
        for subtype in subtypes:

            query="""
            CREATE TABLE 
                {subtype}{type}_meta
            AS
                (
                    SELECT
                        now() as creation_date,
                        '{photo_sweep}'::varchar as photo_sweep,
                        '{photo_resolve}'::varchar as photo_resolve,
                        count(*)
                    FROM
                        {subtype}{type}
                )
            """.format(type=self._type,
                       subtype=subtype,
                       photo_sweep=photo_sweep,
                       photo_resolve=photo_resolve)

            stdout.write(query+'\n')

            pg.execute(query)

            grant_query = 'grant select on {subtype}{type}_meta to public'
            grant_query=grant_query.format(type=self._type,
                                           subtype=subtype)
            stdout.write(grant_query+'\n')
            pg.execute(grant_query)






    def dered_fluxes(self, flux, ivar, extinction):
        """
        Adam says this is more accurate at the faint end.  Why?
        """
        exponent = 0.4*extinction
        flux_correct = (10.0**exponent)*flux

        exponent = -0.8*extinction
        ivar_correct = (10.0**exponent)*ivar

        return flux_correct, ivar_correct

    def get_subset_and_derived(self,data):
        """
        Add columns such as photoid, htmid, mask flags and dereddened model
        mags
        """
        # create the output with the tags we want
        new = numpy.zeros(data.size, dtype=self.get_dtype())
        esutil.numpy_util.copy_fields(data, new)

        new['photoid'] = sdsspy.photoid(data)


        primflag = sdsspy.flags.flagval('resolve_status','survey_primary')
        w,=where((new['resolve_status'] & primflag) != 0)
        if w.size > 0:
            new['isprimary'][w] = 1


        flux,ivar = self.dered_fluxes(data['modelflux'], 
                                      data['modelflux_ivar'], 
                                      data['extinction'])

        mag,err = sdsspy.util.nmgy2mag(flux, ivar=ivar)

        new['modelmag_dered'] = mag
        new['modelmag_dered_err'] = err

        self.load_stomp_maps()
        self.load_htm()

        new['htmid10'] = SweepStuffer._htm.lookup_id(data['ra'], data['dec'])

        new['basic_maskflags'] = SweepStuffer._basic_map.Contains(data['ra'],data['dec'],'eq') 
        new['tycho_maskflags'] = SweepStuffer._tycho_map.Contains(data['ra'],data['dec'],'eq') 


        return new

    def matchruns(self,runs=None, reload=False):
        """
        Return a list of runs and reruns.  If runs= is input, these are
        matched against those from the window_runlist
        """
        self.load_runlist(reload=reload)

        if runs is None:
            return SweepStuffer._runs, SweepStuffer._reruns

        uruns = numpy.unique1d(runs)

        m1, m2 = esutil.numpy_util.match(uruns, SweepStuffer._runs)

        if m1[0] != -1:
            match_runs = SweepStuffer._runs[m2]
            match_reruns = SweepStuffer._reruns[m2]
            return match_runs, match_reruns
        else:
            # just return arrays [-1]
            return m1,m2

    def load_runlist(self, reload=False):
        """
        Load the window runlist
        """
        if not hasattr(SweepStuffer, '_runs') or reload:
            stdout.write("Getting run list, minscore=%s\n" % self._minscore)
            win = sdsspy.window.Window()
            SweepStuffer._runs, SweepStuffer._reruns = win.runlist(self._minscore)

    def load_stomp_maps(self):
        """
        Load the basic and tycho stomp maps
        """
        if not hasattr(SweepStuffer,'_basic_map'):
            stdout.write("loading boss basic stomp map\n")
            SweepStuffer._basic_map = sdsspy.stomp_maps.load('boss','basic')

            stdout.write("loading boss tycho stomp map\n")
            SweepStuffer._tycho_map = sdsspy.stomp_maps.load('boss','tycho')

    def load_htm(self):
        """
        Instantiate an HTM object
        """
        if not hasattr(SweepStuffer,'_htm'):
            stdout.write('loading HTM at depth: %s\n' % self._htm_depth)
            SweepStuffer._htm = esutil.htm.HTM(self._htm_depth)

    def get_dtype(self):
        """

        Removed aperflux from gal

        not in gal:
          aperflux6
          colv
          colverr
          objc_rowc
          rowv
          rowverr
          skyflux
          tmass_bl_flg
          tmass_cc_flg
          tmass_gal_contam
          tmass_h
          tmass_h_ivar
          tmass_j
          tmass_j_ivar
          tmass_jdate
          tmass_k
          tmass_k_ivar
          tmass_matchdist
          tmass_mp_flg
          tmass_nmatch
          tmass_ph_qual
          tmass_rd_flg
          zhedflag
        not in star:
          aperflux
          dev_lnl
          devflux
          devflux_ivar
          exp_lnl
          expflux
          expflux_ivar
          fracpsf
          m_rr_cc
          m_rr_cc_psf
          petroflux
          petroflux_ivar
          petror50
          petror90
          r_dev
          r_exp
          star_lnl
        """

        if self._type == 'gal':
            dtype=[('photoid','i8'),
                   ('run', 'i2'),
                   ('rerun', 'i2'),
                   ('camcol', 'i2'),
                   ('field', 'i2'),
                   ('id', 'i2'),
                   ('thing_id','i4'),
                   ('isprimary','i2'),
                   ('objc_type', 'i4'),
                   ('objc_flags', 'i4'),
                   ('objc_flags2', 'i4'),
                   ('rowc', 'f4', 5),
                   ('colc', 'f4', 5),
                   ('petror50', 'f4', 5),
                   ('petror90', 'f4', 5),
                   ('m_rr_cc', 'f4', 5),
                   ('m_rr_cc_psf', 'f4', 5),
                   ('r_dev', 'f4', 5),
                   ('r_exp', 'f4', 5),
                   ('star_lnl', 'f4', 5),
                   ('exp_lnl', 'f4', 5),
                   ('dev_lnl', 'f4', 5),
                   ('fracpsf', 'f4', 5),
                   ('flags', 'i4', 5),
                   ('flags2', 'i4', 5),
                   ('psp_status', 'i4', 5),
                   ('ra', 'f8'),
                   ('dec', 'f8'),
                   ('psf_fwhm', 'f4', 5),
                   ('extinction', 'f4', 5),
                   ('psfflux', 'f4', 5),
                   ('psfflux_ivar', 'f4', 5),
                   ('fiberflux', 'f4', 5),
                   ('fiberflux_ivar', 'f4', 5),
                   ('fiber2flux', 'f4', 5),
                   ('fiber2flux_ivar', 'f4', 5),
                   ('modelflux', 'f4', 5),
                   ('modelflux_ivar', 'f4', 5),
                   ('modelmag_dered','f4',5),
                   ('modelmag_dered_err','f4',5),
                   ('petroflux', 'f4', 5),
                   ('petroflux_ivar', 'f4', 5),
                   ('devflux', 'f4', 5),
                   ('devflux_ivar', 'f4', 5),
                   ('expflux', 'f4', 5),
                   ('expflux_ivar', 'f4', 5),
                   ('calib_status', 'i2', 5),
                   ('nmgypercount', 'f4', 5),
                   ('resolve_status', 'i2'),
                   ('ifield', 'i4'),
                   ('balkan_id', 'i4'),
                   ('ndetect', 'i2'),
                   ('basic_maskflags','i2'),
                   ('tycho_maskflags','i2'),
                   ('htmid10','i4')]

        elif self._type == 'star':
            dtype=[('photoid','i8'),
                   ('run', 'i2'),
                   ('rerun', 'i2'),
                   ('camcol', 'i2'),
                   ('field', 'i2'),
                   ('id', 'i2'),
                   ('thing_id','i4'),
                   ('isprimary','i2'),
                   ('objc_type', 'i4'),
                   ('objc_flags', 'i4'),
                   ('objc_flags2', 'i4'),
                   ('objc_rowc', 'f4'),
                   ('rowv', 'f4'),
                   ('rowverr', 'f4'),
                   ('colv', 'f4'),
                   ('colverr', 'f4'),
                   ('rowc', 'f4', 5),
                   ('colc', 'f4', 5),
                   ('flags', 'i4', 5),
                   ('flags2', 'i4', 5),
                   ('psp_status', 'i4', 5),
                   ('ra', 'f8'),
                   ('dec', 'f8'),
                   ('psf_fwhm', 'f4', 5),
                   ('extinction', 'f4', 5),
                   ('skyflux', 'f4', 5),
                   ('psfflux', 'f4', 5),
                   ('psfflux_ivar', 'f4', 5),
                   ('fiberflux', 'f4', 5),
                   ('fiberflux_ivar', 'f4', 5),
                   ('fiber2flux', 'f4', 5),
                   ('fiber2flux_ivar', 'f4', 5),
                   ('modelflux', 'f4', 5),
                   ('modelflux_ivar', 'f4', 5),
                   ('modelmag_dered','f4',5),
                   ('modelmag_dered_err','f4',5),
                   ('calib_status', 'i2', 5),
                   ('nmgypercount', 'f4', 5),
                   ('tmass_j', 'f4'),
                   ('tmass_j_ivar', 'f4'),
                   ('tmass_h', 'f4'),
                   ('tmass_h_ivar', 'f4'),
                   ('tmass_k', 'f4'),
                   ('tmass_k_ivar', 'f4'),
                   ('tmass_ph_qual', '|S3'),
                   ('tmass_rd_flg', 'i2'),
                   ('tmass_bl_flg', 'i2'),
                   ('tmass_cc_flg', '|S3'),
                   ('tmass_gal_contam', '|u1'),
                   ('tmass_mp_flg', '|u1'),
                   ('tmass_jdate', 'f8'),
                   ('tmass_matchdist', 'f4'),
                   ('tmass_nmatch', 'i2'),
                   ('resolve_status', 'i2'),
                   ('ifield', 'i4'),
                   ('balkan_id', 'i4'),
                   ('ndetect', 'i2'),
                   ('aperflux6', 'f4', 5),
                   ('zhedflag', 'i2'),
                   ('basic_maskflags','i2'),
                   ('tycho_maskflags','i2'),
                   ('htmid10','i4')]

        return dtype


    def type_differences(self, star=None, gal=None):
        # read in example structures
        run=94
        rerun=301
        camcol=2
        if gal is None:
            gal = sdsspy.files.read('calibobj',run,rerun,camcol=camcol,
                                    type='gal', ensure_native=True)
        if star is None:
            star = sdsspy.files.read('calibobj',run,rerun,camcol=camcol,
                                     type='star', ensure_native=True)

        # find tags not in both
        galnames = list(gal.dtype.names)
        starnames = list(star.dtype.names)
        allnames = numpy.array(galnames + starnames, dtype='S25')
        allnames = numpy.unique1d( allnames )

        notingal = []
        notinstar = []
        for n in allnames:
            # is the name in gal?
            if galnames.count(n) == 0:
                notingal.append(n)
            if starnames.count(n) == 0:
                notinstar.append(n)
        stdout.write("not in gal:\n")
        stdout.write('-'*72 + '\n')
        for n in notingal:
            stdout.write('  %s\n' % n)
        stdout.write("not in star:\n")
        stdout.write('-'*72 + '\n')
        for n in notinstar:
            stdout.write('  %s\n' % n)

        return gal,star
