import numpy
from . import files
from . import pipe

def sweep_camcol(**keys):
    sm=SweepMaker(**keys)
    sm.make_sweep()

class SweepMaker(dict):
    def __init__(self, **keys):
        self.update(keys)

        self._url=files.get_sweep_url(**self)
        self.conf=files.read_config(self['gmix_run'])

    def make_sweep(self):
        import esutil as eu
        fields=self._get_fields()

        data_list=[]
        for field in fields:
            d=files.read_output(verbose=True, field=field, **self)
            
            if d.size==0:
                continue

            w=self._select(d)
            if w.size==0:
                continue

            d=d[w]
            out=self._get_output(d)
            data_list.append(out)

        if len(data_list) == 0:
            data=None
        else:
            data=eu.numpy_util.combine_arrlist(data_list)

        print 'writing:',self._url
        eu.io.write(self._url,data,clobber=True)
 
    def _get_output(self, data):
        """
        Should be trimmed to objects with successful measurements

        Returns only the best model
        """
        import esutil as eu
        npars_psf=data['psf_pars'].shape[1]
        front=self.conf['obj_models'][0]
        front += '_'
        npars=data[front+'pars'].shape[1]

        dt=pipe.get_dtype(self.conf['obj_fitter'],
                          npars,
                          npars_psf,
                          noflags=True,
                          models=[''])
        dt.append( ('model','S4') )

        st=numpy.zeros(data.size, dtype=dt)
        eu.numpy_util.copy_fields(data, st)


        for i in xrange(data.size):
            aic_min=9.999e40
            for model in self.conf['obj_models']:
                front='%s_' % model
                aic = data[front+'aic'][i]
                if aic < aic_min:
                    use_model=model
                    aic_min=aic

            st['model'][i] = use_model

            front='%s_' % use_model
            for n in st.dtype.names:
                nn=front+n

                if nn in data.dtype.names:
                    st[n][i] = data[nn][i]

        return st

    def _select(self, objs):

        logic = (objs['flags']==0)

        # we only keep if all model fits succeed
        for model in self.conf['obj_models']:
            front='%s_' % model
            n=front+'flags'
            logic = logic & (objs[n] == 0)

        keep,=numpy.where(logic)
        return keep


    def _get_fields(self):
        primfields=files.read_field_cache(gmix_run=self['gmix_run'])
        w,=numpy.where((primfields['run']==self['run']) 
                      & (primfields['camcol'] == self['camcol']) )

        if w.size==0:
            raise ValueError("found no fields for run=%s "
                             "camcol=%s" % (self['run'],self['camcol']) )

        fields=primfields['field'][w]

        return fields

