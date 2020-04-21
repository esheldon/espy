"""
Module:
    window
Classes:

    Window

    Purpose:
        Class for working with the SDSS window function.  This parallels the
        functions in the sdss PHOTOOP package.   See the docs for window.Window
        for full details.

"""

import os
import sys
from sys import stdout, stderr

import numpy
import esutil
#import mangle


class Window():
    """
    Class:
        Window
    Purpose:
        Class for working with the SDSS window function.  This parallels the
        fuwctions in the sdss PHOTOOP package.
    Methods:
        read: 
            Read information about the SDSS window.  The information that is
            returned is governed by the types argument.  
        runlist:
            Return the runs,reruns corresponding to the request.  The runs can
            be trimmed by their score.
        window_types:
            Return a list of available file types, those which can be read
            using the read() method.
        name:
            Return the path to a window file.
        
    """
    def __init__(self):
        self._types = ['flist','blist','balkans','bcaps','findx','bindx']
        resolve_dir = os.getenv('PHOTO_RESOLVE')
        if resolve_dir is None:
            raise ValueError('$PHOTO_RESOLVE is not set')
        self._resolve_dir = resolve_dir


    def read(self, types, verbose=False):
        """
        Class:
            Window
        Method Name:
            read
        Purpose:
            Read information about the SDSS window.  The information that is
            returned is governed by the types argument.  

        Calling Sequence:
            import sdsspy
            w=sdsspy.window.Window()
            data = w.read(types, verbose=False)

        Inputs:
            types:  
                The type of info to read.  Available types can be got with
                window_types() method, and are:
                    ['flist','blist','balkans','bcaps','findx','bindx']
            verbose=: If True, print info about files read.

        Outputs:
            A dictionary with keys corresponding to the requested types. e.g
                data = w.read(['flist','bcaps'])
                data.keys()
                ['bcaps','flist']
            The data in each key is a numpy array.

        Requirements:
            Requires the esutil package.  http://code.google.com/p/esutil/

        Examples:
            import sdsspy
            w=sdsspy.window.Window()
            data = w.read(['flist','bindx'])

        """
        output = {}
        if isinstance(types, (str,unicode)):
            types = [types]

        if 'flist_rescore' in types:
            rfile = self.name('flist_rescore')
            if not os.path.exists(rfile):
                stdout.write("Generating rescore file\n")
                # still need to write the rescore method
                raise ValueError("implement scoring")
                self.rescore()
                if not os.path.exists(rfile):
                    raise RuntimeError("Could not generate rescore")

            output['flist_rescore'] = \
                esutil.io.read(rfile, ext=1, lower=True, verbose=verbose)
        elif 'flist' in types:
            rfile = self.name('flist')
            output['flist'] = esutil.io.read(rfile, ext=1, lower=True, verbose=verbose)

        
        if 'blist' in types or 'balkans' in types:
            fname = self.name('blist')
            output['blist'] = esutil.io.read(fname, ext=1, lower=True, verbose=verbose)

        if 'bcaps' in types or 'balkans' in types:
            fname = self.name('bcaps')
            output['bcaps'] = esutil.io.read(fname, ext=1, lower=True, verbose=verbose)

        if 'findx' in types:
            fname = self.name('findx')
            output['findx'] = esutil.io.read(fname, ext=1, lower=True, verbose=verbose)

        if 'bindx' in types:
            fname = self.name('bindx')
            output['bindx'] = esutil.io.read(fname, ext=1, lower=True, verbose=verbose)

        if 'balkans' in types:
            if verbose:
                stdout.write('Creating balkans\n')
            dtype = output['blist'].dtype.descr
            dtype += [('use_caps','i8'),('caps',numpy.object)]

            balkans = numpy.zeros(output['blist'].size, dtype=dtype)

            esutil.numpy_util.copy_fields(output['blist'], balkans)

            for i in range(balkans.size):
                if verbose:
                    if ( (i+1) % 10000 ) == 0:
                        stdout.write("%s/%s\n" % ( (i+1), balkans.size))
                # no copy is made
                tmp_balkan = balkans[i]
                ncaps = tmp_balkan['ncaps']
                # 'caps' is object array
                #tmp_balkan['caps'] = mangle.construct_cap(ncaps)
                tmp_balkan['caps'] = \
                    output['bcaps'][tmp_balkan['icap']:tmp_balkan['icap']+ncaps]
                #for j in range(ncaps):
                    #tmp_balkan['caps'][j] = output['bcaps'][tmp_balkan['icap']+j]


                mangle.set_use_caps(tmp_balkan, 
                                    range(tmp_balkan['ncaps']), 
                                    allow_doubles=True)

                # gotta do this or the memory blows up, and the garbage collector
                # really slows things down
                del tmp_balkan

            output['balkans'] = balkans

        if len(output) == 0:
            mess="""None of the requested types are valid.  Must be one of \n"""
            mess += "%s "% self._types
            raise ValueError(mess)

        return output            

    def runlist(self, minscore=None, rescore=False):
        """
        Class:
            Window
        Method Name:
            runlist
        Purpose:
            Return the runs,reruns corresponding to the request.  The runs can
            be trimmed by their score.
        Calling Sequence:
            import sdsspy
            w=sdsspy.window.Window()
            runs, rerun = w.runlist(minscore=None, rescore=False)

        Optional Inputs:
            minscore:  
                Only runs,reruns with score > minscore will be returned.  Note
                the version in IDL PHOTOOP uses >=, so be careful when
                converting your code.

            rescore: If True, re-score the runs.
        """

        if rescore:
            type = 'flist_rescore'
        else:
            type = 'flist'
        data = self.read(type)

        flist = data[type]

        runs,reruns = flist['run'], flist['rerun']
        if minscore is not None:
            w, = numpy.where( flist['score'] > minscore)

            if w.size == 0:
                return None,None

            runs=runs[w]
            reruns=reruns[w]

        # using a bigid is not necessary since the same run never ends up in
        # there with two different reruns

        uid, uindex = numpy.unique1d(runs, return_index=True)
        runs = runs[uindex]
        reruns = reruns[uindex]

        return runs, reruns


    def window_types(self):
        """
        Class:
            Window
        Method Name:
            window_types
        Purpose:
            return the available file types for window.
        Calling Sequence:
            import sdsspy
            w=sdsspy.window.Window()
            name = w.window_types()

            This should return:
                ['flist','blist','balkans','bcaps','findx','bindx']
        """

        return self._types

    def name(self, type):
        """
        Class:
            Window
        Method Name:
            name
        Purpose:
            return a full path to a window file of the given type.
        Calling Sequence:
            import sdsspy
            w=sdsspy.window.Window()
            name = w.name(type)
        Inputs:
            type: 
                One of the window types.  See window_types() for a list, which
                should be one of
                    ['flist','blist','balkans','bcaps','findx','bindx']
        """
        if type not in self._types:
            raise ValueError("Bad window type: '%s'" % type)

        fname = 'window_%s.fits' % type

        fname = os.path.join(self._resolve_dir,fname)
        return fname



