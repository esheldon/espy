"""

Simple routines to externally call the C/fortran mangle routines.  For complex,
pixelized masks this is actually faster than using the mangle python routines 
due to some optimizations that have not yet been included in the python.

Calls my shell script mpolyid for the calculations
"""

import os
import esutil as eu
import tempfile
import subprocess
import numpy

class Mangle:
    """
    Do mangle calculations

    Just calls mangle externally

    parameters
    ----------
    mangle_url: string
        The location of the mangle file
    veto: bool, optional
        If True, this is a veto mask.  The contains() function will thus treat
        objects with polyid of -1 (masked) as "contained".
    """

    def __init__(self, mask_url, veto=False):
        self.mask_url=mask_url
        self.is_veto_mask=veto

        if not os.path.exists(mask_url):
            raise ValueError("mask does not exist: %s" % mask_url)

    def contains(self, ra, dec):
        """
        Determine of the input ra,dec points are masked.

        parameters
        ----------
        mangle_url: string
            The location of the mangle file
        ra: array or scalar
            right ascension in degrees
        dec: array or scalar
            declination in degrees

        output
        ------
        A numpy array, 1 if unmasked 0 if masked.  Note for veto
        masks, you will "keep" points that are 0 in this.
        """

        ids=self.polyid(ra, dec)
        if self.is_veto_mask:
            return numpy.where(ids >= 0, 0, 1)
        else:
            return numpy.where(ids >= 0, 1, 0)

    def polyid(self, ra, dec):
        """
        Calculate the polygon id for the input ra,dec and mask

        parameters
        ----------
        mangle_url: string
            The location of the mangle file
        ra: array or scalar
            right ascension in degrees
        dec: array or scalar
            declination in degrees

        output
        ------
        A numpy array with the polygon id for each point. -1 if
        the point is masked.  Note for veto masks, you will want
        to "keep" the -1 points.
        """

        if numpy.isscalar(ra):
            ra=numpy.array(ra, ndmin=1, dtype='f8')
        if numpy.isscalar(dec):
            dec=numpy.array(dec, ndmin=1, dtype='f8')

        if len(ra) != len(dec):
            raise ValueError("ra dec must be same length")
        
        return self._call_polyid(ra,dec)

    def _call_polyid(self, ra, dec):
        radecfile=tempfile.mktemp(prefix='tmp-mangle-', suffix='-radec.dat')
        tmpfile=tempfile.mktemp(prefix='tmp-mangle-', suffix='-ids-pre.dat')
        idfile=tempfile.mktemp(prefix='tmp-mangle-', suffix='-ids.dat')

        try:
            with open(radecfile,'w') as fobj:
                for i in xrange(ra.size):
                    fobj.write("%.16g %.16g\n" % (ra[i],dec[i]))

            command='polyid -p+16-16 {mask} {radecfile} {tmpfile}'.format(mask=self.mask_url,
                                                                          radecfile=radecfile,
                                                                          tmpfile=tmpfile)
            # we want the stdout/err to go to the terminal here
            code,sout,serr=exec_process(command, stdout_file=None, stderr_file=None)
            if code != 0:
                raise ValueError("error running\n\t%s" % command)

            command="awk '(NF>3)' {tmpfile} | wc -l".format(tmpfile=tmpfile)

            code,sout,serr=exec_process(command)
            if code != 0:
                raise ValueError("error running id check\n\t%s\n\terr:%s" % (command,serr))
            nbad = int(sout)
            if nbad != 0:
                raise ValueError("found %s objects with > 1 polygons")


            command="awk '(NR>1) {if (NF==3) print $3; else print -1;}' %s > %s"
            command=command % (tmpfile,idfile)
            # we want the stderr to go to the terminal here
            code,sout,serr=exec_process(command, stderr_file=None)
            if code != 0:
                raise ValueError("error running\n\t%s" % command)
            
            ids = numpy.fromfile(idfile, dtype='i8', sep=' ')

        finally:
            # don't want files left sitting around in /tmp
            for f in [radecfile,tmpfile,idfile]:
                if os.path.exists(f):
                    os.remove(f)

        return ids

def exec_process(cmd, 
                 stdout_file=subprocess.PIPE, 
                 stderr_file=subprocess.PIPE):

    pobj = subprocess.Popen(cmd, 
                            stdout=stdout_file, 
                            stderr=stderr_file, 
                            shell=True)

    stdout_ret, stderr_ret = pobj.communicate()
    exit_status = pobj.returncode

    return exit_status, stdout_ret, stderr_ret
