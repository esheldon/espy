import os
import sys
from sys import stdout, stderr

import numpy
from numpy import where

# this is unlikely to import unless someone
# compiles it!
try:
    import pgnumpy
except:
    pass

import esutil
from esutil.ostools import path_join, expand_path, getenv_check

import sdsspy
import sdssgal



def pg_create_input(run_name):
    ms = MaxbcgPgSelector(run_name)
    ms.process_runs()

# This class used the postgres database as the starting point
class MaxbcgPgSelector:
    def __init__(self, run_name,
                 max_rmag=22.0, 
                 max_imag=21.0, 
                 minscore=0.1):

        self.maskname = 'tycho'

        self.run_name=run_name
        
        self.max_rmag = max_rmag
        self.max_imag = max_imag
        self.minscore=minscore
        self.flags = sdsspy.flags.Flags()

        win = sdsspy.window.Window()
        stdout.write("Loading run list\n")
        self.runs,reruns = win.runlist(minscore=self.minscore)

        stdout.write("connecting to boss database\n")
        self.pg = pgnumpy.connect('dbname=boss')

    def process_runs(self, **keys):
        d=self.outdir()
        if not os.path.exists(d):
            os.makedirs(d)

        outfile = self.outfile()
        stdout.write("Will write to file: %s\n" % outfile)
        if os.path.exists(outfile):
            stdout.write("Removing existing file: %s\n" % outfile)
            os.remove(outfile)

        i=1
        ntot=len(self.runs)

        for run in sorted(self.runs):
            stdout.write("Run: %s  %s/%s\n" % (run,i,ntot))
            stdout.write('-'*70 + '\n')
            if i == 1:
                # more verbose the first time
                stdout.write(self.run_query(run))
                stdout.write('\n')
            objs = self.read_run(run)

            if objs.size > 0:

                self.cmodelmag_dered,self.cmodelmag_dered_err = \
                        sdsspy.make_cmodelmag(objs, dered=True)

                stdout.write("    Getting logic\n")
                mag_logic = self.mag_logic()
                binned_logic = self.binned_logic(objs)
                object1_logic = self.object1_logic(objs)

                logic = mag_logic & binned_logic & object1_logic

                stdout.write("    Selecting good galaxies\n")
                keep, = where(logic)

                nkeep = keep.size
                stdout.write("    Keeping %s/%s\n" % (nkeep,objs.size))

                if nkeep > 0:
                    stdout.write("        Extracting fields...")
                    output= self.extract_fields(objs, keep)
                    stdout.write("Writing data...")
                    esutil.sfile.write(output, outfile,append=True)
                    stdout.write("Done\n")

            i+=1


            

    def extract_fields(self, objs, keep):
        """
        """
        dtype=[('photoid','i8'),
               ('ra','f8'),
               ('dec','f8'),
               ('modelmag_dered','f4',5),
               ('modelmag_dered_err','f4',5),
               ('cmodelmag_dered','f4',5),
               ('cmodelmag_dered_err','f4',5)]

        output = numpy.zeros(keep.size, dtype=dtype)

        output['photoid'] = objs['photoid'][keep]
        
        output['ra'] = objs['ra'][keep]
        output['dec'] = objs['dec'][keep]

        f, ivar = sdsspy.dered_fluxes(objs['extinction'][keep],
                                      objs['modelflux'][keep],
                                      objs['modelflux_ivar'][keep])
        mag,err = sdsspy.nmgy2mag(f,ivar)
        output['modelmag_dered'] = mag
        output['modelmag_dered_err'] = err

        output['cmodelmag_dered'] = self.cmodelmag_dered[keep,:]
        output['cmodelmag_dered_err'] = self.cmodelmag_dered_err[keep,:]

        return output

    def read_run(self, run, verbose=False, usemask=True):
        q=self.run_query(run, usemask=usemask)
        if verbose:
            stdout.write(q)
        data = self.pg.fetch(q)
        if data.size == 0:
            stdout.write("    No primary galaxies found in %s window "
                         "for run: %s\n" % (self.maskname,run))
        return data

    def run_query(self,run,usemask=True):
        """
        Select primary and inside the tycho window
        """

        extra=''
        if usemask:
            extra="and {maskname}_maskflags=1".format(maskname=self.maskname)
        query="""
        SELECT
            photoid,
            ra,
            dec,
            flags,
            flags2,
            objc_flags,
            objc_flags2,
            modelflux,
            modelflux_ivar,
            devflux,
            devflux_ivar,
            expflux,
            expflux_ivar,
            fracpsf,
            extinction
        FROM
            primgal
        WHERE
            run={run}
            {extra}
        """.format(run=run,extra=extra)
        return query
 

    def mag_logic(self):
        rmag_logic = self.cmodelmag_dered[:,2] < self.max_rmag
        imag_logic = self.cmodelmag_dered[:,3] < self.max_imag
        mag_logic = rmag_logic & imag_logic
        return mag_logic


    def object1_logic(self,objs):
        satur = self.flags.val('object1','satur')
        bright = self.flags.val('object1','bright')
        too_many_peaks = self.flags.val('object1','deblend_too_many_peaks')

        blended = self.flags.val('object1','blended')
        nodeblend = self.flags.val('object1','nodeblend')

        peakcenter = self.flags.val('object1','peakcenter')
        notchecked = self.flags.val('object1','notchecked')
        noprofile = self.flags.val('object1','noprofile')

        oflags = objs['objc_flags']

        # famously worded as double negatives
        bdb_logic = \
            ( (oflags & blended) == 0) | ((oflags & nodeblend) != 0)

        # now combine logic
        logic = \
            ((oflags & satur) == 0) \
            & ((oflags & bright) == 0) \
            & ((oflags & too_many_peaks) == 0) \
            & ((oflags & peakcenter) == 0) \
            & ((oflags & notchecked) == 0) \
            & ((oflags & noprofile) == 0) \
            & bdb_logic

        return logic

    def binned_logic(self, objs):
        binned1 = self.flags.val('object1','binned1')
        binned2 = self.flags.val('object1','binned2')
        binned4 = self.flags.val('object1','binned4')

        rflags = objs['flags'][:,2]
        iflags = objs['flags'][:,3]

        r_binned_logic = \
            ((rflags & binned1) != 0 )  \
            | ((rflags & binned2) != 0 )  \
            | ((rflags & binned4) != 0 )

        i_binned_logic = \
            ((iflags & binned1) != 0 )  \
            | ((iflags & binned2) != 0 )  \
            | ((iflags & binned4) != 0 )

        return r_binned_logic & i_binned_logic




    def outdir(self):
        maxbcg_dir=os.environ['MAXBCG_INPUT']
        d = path_join(maxbcg_dir, self.run_name)
        return d
    def outfile(self):
        d=self.outdir()
        f='maxbcg-input-%s-%s-r%0.1f-i%0.1f.rec'
        f=f % (self.run_name,self.maskname, self.max_rmag,self.max_imag)
        f=path_join(d,f)
        return f

