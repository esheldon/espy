from __future__ import print_function
import sdsspy
import es_sdsspy
import sdsspy
import numpy




class Selector:
    '''
    Select Methods:
        These return boolean logic arrays.

        mask_logic(masktype)
        cmodelmag_logic(filter, limit)
        resolve_logic()

        binned_logic = self.binned_logic()
        calib_logic = self.calib_logic()
        object1_logic = self.object1_logic()

        flag_logic(): combines binned_logic, calib_logic, object1_logic


    if just selecting primary:

        selector = Selector(struct)
        logic = selector.resolve_logic()

        keep, = numpy.where(logic)

        keepstruct = struct[keep]

    Multiple selection

        mag_logic = selector.cmodelmag_logic('r', 22.0)
        resolve_logic = selector.resolve_logic()
        flag_logic = selector.flag_logic()

        logic = mag_logic & resolve_logic & flag_logic

        keep, = numpy.where(logic)

    Spatial selections:

        basic_logic = selector.mask_logic('basic')
        good_logic = selector.mask_logic('good')
        tycho_logic = selector.mask_logic('tycho')

    '''
    def __init__(self, objs, **keys):

        self.objs = objs
        self.flags = sdsspy.flags.Flags()
        self.masks={}
        self.mags={}

    def mask_logic(self, masktype):
        self.load_mask(masktype)

        ra,dec=self.objs['ra'],self.objs['dec']
        m = self.masks[masktype]
        if hasattr(m, 'Contains'):
            inmask = m.Contains(ra,dec,'eq')
        else:
            inmask = m.contains(ra,dec)
        return (inmask == 1)

    def load_mask_stomp(self, masktype):
        if masktype not in self.masks:
            self.masks[masktype] = \
                es_sdsspy.stomp_maps.load('boss',masktype,verbose=True)
    def load_mask(self, masktype):
        if masktype not in self.masks:
            self.masks[masktype] = \
                es_sdsspy.mangle_masks.load('boss',masktype,verbose=True)



    def modelmag_logic(self, filter, limit, dered=True):
        fnum = sdsspy.FILTERNUM[filter]
        if dered:
            name='cmodel_dered'
        else:
            name='cmodel'
        
        # cache the mags
        if name not in self.mags:
            mag = sdsspy.nmgy2mag(self.objs['modelflux'])
            if dered:
                mag -= self.objs['extinction']
            self.mags[name] = mag

        mag_logic = self.mags[name][:,fnum] < limit
        return mag_logic


    def cmodelmag_logic(self, filter, limit, dered=True):
        fnum = sdsspy.FILTERNUM[filter]
        if dered:
            name='cmodel_dered'
        else:
            name='cmodel'
        
        # cache the mags
        if name not in self.mags:
            self.mags[name] = sdsspy.make_cmodelmag(self.objs, doerr=False, dered=dered)
        mag_logic = self.mags[name][:,fnum] < limit
        return mag_logic

    def resolve_logic(self):
        primary_flag = self.flags.val('resolve_status','survey_primary')
        return (self.objs['resolve_status'] & primary_flag) != 0


    def flag_logic(self):

        # make sure detected in r and i in binned1,2, or 4
        binned_logic = self.binned_logic()

        # make sure no non-photometric data
        calib_logic = self.calib_logic()

        # checking lots of object1 flags
        object1_logic = self.object1_logic()

        # combined logic
        logic = binned_logic & calib_logic & object1_logic

        return logic

    def object1_logic(self):
        if 'objc_flags' not in self.objs.dtype.names:
            raise ValueError("need objc_flags in struct")


        satur = self.flags.val('object1','satur')
        bright = self.flags.val('object1','bright')
        too_many_peaks = self.flags.val('object1','deblend_too_many_peaks')

        blended = self.flags.val('object1','blended')
        nodeblend = self.flags.val('object1','nodeblend')

        peakcenter = self.flags.val('object1','peakcenter')
        notchecked = self.flags.val('object1','notchecked')
        noprofile = self.flags.val('object1','noprofile')

        oflags = self.objs['objc_flags']

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


    def binned_logic(self):
        binned1 = self.flags.val('object1','binned1')
        binned2 = self.flags.val('object1','binned2')
        binned4 = self.flags.val('object1','binned4')

        if 'flags' in self.objs.dtype.names:
            rflags = self.objs['flags'][:,2]
            iflags = self.objs['flags'][:,3]
        elif 'flags_r' in self.objs.dtype.names:
            rflags = self.objs['flags_r']
            iflags = self.objs['flags_i']
        else:
            raise ValueError("need flags or flags_{band} in struct")


        r_binned_logic = \
            ((rflags & binned1) != 0 )  \
            | ((rflags & binned2) != 0 )  \
            | ((rflags & binned4) != 0 )

        i_binned_logic = \
            ((iflags & binned1) != 0 )  \
            | ((iflags & binned2) != 0 )  \
            | ((iflags & binned4) != 0 )

        return r_binned_logic & i_binned_logic




    def calib_logic(self):
        objs = self.objs

        if 'calib_status' in objs.dtype.names:
            calib_status_u = objs['calib_status'][:,0]
            calib_status_g = objs['calib_status'][:,1]
            calib_status_r = objs['calib_status'][:,2]
            calib_status_i = objs['calib_status'][:,3]
            calib_status_z = objs['calib_status'][:,4]
        elif 'calib_status_r' in objs.dtype.names:
            calib_status_u = objs['calib_status_u']
            calib_status_g = objs['calib_status_g']
            calib_status_r = objs['calib_status_r']
            calib_status_i = objs['calib_status_i']
            calib_status_z = objs['calib_status_z']
        else:
            raise ValueError("need calib_status or calib_status_{band} in struct")


        flagval = self.flags.val('calib_status', 'photometric')
        logic = \
            ((calib_status_u & flagval) != 0)   \
            & ((calib_status_g & flagval) != 0) \
            & ((calib_status_r & flagval) != 0) \
            & ((calib_status_i & flagval) != 0) \
            & ((calib_status_z & flagval) != 0)


        return logic
        



class SweepSelector(object):
    """
    Take a photoObj catalog and trim it to be identical to the objects in the
    sweeps
    """
    def __init__(self, objs, type, trim=False):
        self._type=type
        self._trim=trim
        self._objs=objs
        self._flags=sdsspy.flags.Flags()

        if self._type!='gal':
            raise ValueError("Only gal implemented currently")

        self._type_info={'gal':{'magcuts':[21., 22., 22., 20.5, 20.1],
                                'typenum':3}}

        self._do_select()

    def get_logic(self):
        return self._logic
    def get_indices(self):
        return self._indices
    
    def _do_select(self):
        objs=self._objs

        type_logic = self._get_type_logic()
        mag_logic = self._get_mag_logic()
        flag_logic=self._get_flag_logic()
        resolve_logic = self._get_resolve_logic()

        self._logic=type_logic & mag_logic & flag_logic & resolve_logic
        self._indices, = numpy.where(self._logic)

    def _get_type_logic(self):
        typenum=self._type_info[self._type]['typenum']

        logic = self._objs['objc_type'] == typenum
        return logic
        
    def _get_flag_logic(self):
        flags=self._flags
        objs=self._objs

        oflags = objs['objc_flags']

        bright = flags.val('object1','bright')
        blended = flags.val('object1','blended')
        nodeblend = flags.val('object1','nodeblend')


        bdb_logic = \
            ( (oflags & blended) == 0) | ((oflags & nodeblend) != 0)

        bright_logic= ((oflags & bright) == 0)

        logic = bdb_logic & bright_logic
        if self._trim:
            edge = flags.val('object1','edge')
            too_many_peaks = flags.val('object1','deblend_too_many_peaks')
            badsky=flags.val('object1','badsky')
            satur=flags.val('object1','satur')
            satcen=flags.val('object2','satur_center')

            oflags2 = objs['objc_flags2']
            trim_logic = (
                ((oflags & too_many_peaks) == 0) 
                & ((oflags & edge)==0)
                & ((oflags & badsky) == 0)
                & ((oflags & satur) == 0)
                & ((oflags2 & satcen) == 0) )

            logic = logic & trim_logic

        return logic

    def _get_resolve_logic(self):
        flags=self._flags
        objs=self._objs

        primary_flag = flags.val('resolve_status','run_primary')
        return (objs['resolve_status'] & primary_flag) != 0

    def _get_mag_logic(self):
        magcuts = self._type_info[self._type]['magcuts']

        logic = numpy.zeros(self._objs.size, dtype='bool')

        for i in xrange(5):

            flux=self._objs['modelflux'][:,i]
            maxflux=flux.max()
            flux=flux.clip(0.0001, maxflux)

            mag=22.5-2.5*numpy.log10(flux)
            mag -= self._objs['extinction'][:,i]

            logic = logic | (mag < magcuts[i])

        return logic

def test_sweep_selector(run=756, camcol=3, type='gal'):
    sweep=sdsspy.read('calibobj.'+type, run=run, camcol=camcol, lower=True)

    fields=numpy.unique(sweep['field'])

    for field in fields:
        pobj=sdsspy.read('photoObj', run=run, camcol=camcol, field=field, 
                         lower=True)

        sel=SweepSelector(pobj, type)
        ind=sel.get_indices()

        w,=numpy.where(sweep['field'] == field)
        ndiff=w.size-ind.size
        print('field:',field,'sweep count:',w.size,'photoObj trimmed count:',ind.size,'diff:',ndiff)

