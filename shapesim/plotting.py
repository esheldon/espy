import os
from . import shapesim
import esutil as eu
from esutil.misc import wlog

class SimPlotter(dict):
    def __init__(self, run):
        self['run'] = run
        c = shapesim.read_config(run)
        for k,v in c.iteritems():
            self[k] = v

        d = shapesim.get_plot_dir(run)
        if not os.path.exists(d):
            os.makedirs(d)
        
        self._data=None

    def doplots(self, s2max=None, yrange=None):
        import biggles
        data = self.read_data(s2max=s2max)
        epsfile = shapesim.get_plot_file(self['run'],
                                         s2max=s2max,
                                         yrange=yrange)
        wlog("will plot to:",epsfile)


    def read_data(self, s2max=None):
        if self._data is None:
            wlog("reading data")
            self._data = shapesim.read_all_outputs(self['run'],average=True)
        
        alldata = self._data
        ntot=len(alldata)
        nkeep=ntot
        if s2max is not None:
            keepdata = []
            for st in alldata:
                if numpy.median(st['s2_meas']) < s2max:
                    keepdata.append(st)
            nkeep = len(keepdata)
            wlog("kept %d/%d with s2 < %.3g" % (nkeep,ntot,s2max)) 
        else:
            keepdata = alldata

        return keepdata
