class StructFits(object):
    def __init__(self, struct_def):
        self.struct_def = struct_def
        self.process_struct_def()

    def get_name_raw(self):
        return self.struct_name_raw
    def get_name_norm(self):
        return self.struct_name_norm

    def process_struct_def(self):
        struct_def=self.struct_def
        sbeg = struct_def.find('{')
        send = struct_def.find('}')
        if sbeg == -1 or send == -1:
            raise ValueError("bad struct def: '%s'" % struct_def)

        self.struct_name_raw = struct_def[0:sbeg].strip()
        self.struct_name_norm = self.struct_name_raw.replace(' ','_')
        
        field_defs = struct_def[sbeg+1:send].split(';')

        field_defs = [f.strip() for f in field_defs] 

        fields=[_field_def_to_field(f) for f in field_defs if f != '']

        self.fields=fields

    def get_names_c(self):
        names=['"'+f.get_name()+'"' for f in self.fields]
        names = '{' + ','.join(names)+'}'
        return names

    def get_tforms_c(self):
        tforms=['"'+f.get_tform()+'"' for f in self.fields]
        tforms = '{' + ','.join(tforms)+'}'
        return tforms


    def get_create_table_signature(self):
        s="""
int %(struct_name_norm)s_create_table(fitsfile* fits, LONGLONG nrows, const char *extname, int *status)
        """ % {'struct_name_norm':self.get_name_norm()}
        return s.strip()


    def get_create_table_func(self):
        names=self.get_names_c()
        tforms=self.get_tforms_c()
        struct_name_norm=self.get_name_norm()

        sig=self.get_create_table_signature()
        func="""
%(sig)s
{
    int nfields=%(nfields)d;
    char **names=%(names)s;
    char **tforms=%(tforms)s;
    char **tunits=NULL;
    fits_create_tbl(fits, BINARY_TBL, nrows, nfields, names, tforms, tunits,
                    extname, status);
}\n""" % {'struct_name_norm':struct_name_norm,
          'nfields':len(self.fields),
          'names':names,
          'tforms':tforms,
          'sig':sig}
        return func

    def get_col_write_statements(self):
        s=[]
        for i,f in enumerate(self.fields):
            s.append(f.get_write_statement(i))

        return '\n        '.join(s)

    def get_write_table_signature(self):
        s="""
int %(struct_name_norm)s_write_fits(fitsfile* fits, %(struct_name_raw)s *self, LONGLONG nrows, int *status)
        """ % {'struct_name_raw':self.get_name_raw(),
               'struct_name_norm':self.get_name_norm()}
        return s.strip()


    def get_write_table_func(self):
        struct_name_raw=self.get_name_raw()
        struct_name_norm=self.get_name_norm()
        col_writes = self.get_col_write_statements()

        sig = self.get_write_table_signature()

        goto_name='_%s_fits_write_bail' % (struct_name_norm)

        # doesn't work for arrays
        func="""
%(sig)s
{
    LONGLONG row=0;
    LONGLONG firstelem=1;
    %(struct_name_raw)s *str=NULL;

    for (row=0; row<nrows; row++) {
        str = &self[row];

        %(col_writes)s
        if (*status) {
            goto %(goto_name)s;
        }
    }

%(goto_name)s:
    if (*status) {
        fits_report_error(stderr,*status);
    }
    return *status;
}\n""" % {'struct_name_raw':struct_name_raw,
          'struct_name_norm':struct_name_norm,
          'col_writes':col_writes,
          'goto_name':goto_name,
          'sig':sig}
        return func



    def get_col_read_statements(self):
        s=[]
        for i,f in enumerate(self.fields):
            s.append(f.get_read_statement(i))

        return '\n        '.join(s)

    def get_read_table_signature(self):
        s="""
%(struct_name_raw)s *%(struct_name_norm)s_read_fits(fitsfile* fits, int *status)
        """ % {'struct_name_raw':self.get_name_raw(),
               'struct_name_norm':self.get_name_norm()}
        return s.strip()

    def get_read_table_func(self):
        struct_name_raw=self.get_name_raw()
        struct_name_norm=self.get_name_norm()
        col_reads = self.get_col_read_statements()

        goto_name='_%s_fits_read_bail' % (struct_name_norm)

        sig = self.get_read_table_signature()

        # doesn't work for arrays
        func="""
%(sig)s
{
    LONGLONG nrows=0, row=0;
    LONGLONG firstelem=1;
    %(struct_name_raw)s *self=NULL, *str=NULL;

    if (fits_get_num_rows(fits, &nrows, status)) {
        goto %(goto_name)s;
    }

    self=calloc(nrows, sizeof(%(struct_name_raw)s));
    if (!self) {
        fprintf(stderr,"could not allocate %%ld %(struct_name_raw)s\\n",nrows);
        exit(1);
    }

    for (row=0; row<nrows; row++) {
        str = &self[row];

        %(col_reads)s
        if (*status) {
            goto %(goto_name)s;
        }

    }

%(goto_name)s:
    if (*status) {
        fits_report_error(stderr,*status);
        free(self);
        self=NULL;
    }
    return self;
}\n"""
        func = func % {'struct_name_raw':struct_name_raw,
                       'struct_name_norm':struct_name_norm,
                       'goto_name':goto_name,
                       'col_reads':col_reads,
                       'sig':sig}
        return func


class FieldBase(object):
    def __init__(self, field_def):
        self.field_def=field_def

        self._set_field_info()

    def is_scalar(self):
        """
        Always true for now
        """
        return True

    def get_name(self):
        return self.f_name
    def get_type(self):
        return self.f_type
    def get_shape(self):
        return self.f_shape
    def get_tform(self):
        return self.tform
    def get_fits_type(self):
        return self.fits_type
    def get_short_name(self):
        return self.short_name
    def get_tdim(self):
        return self.tdim

class StringField(FieldBase):
    def __init__(self, field_def):
        super(StringField,self).__init__(field_def)

        if '][' in field_def:
            raise ValueError("arrays not supported: '%s'" % field_def)

    def get_write_func(self):
        return 'fits_write_col_str'
    def get_read_func(self):
        return 'fits_read_col_str'

    def get_write_statement(self, colnum):
        write_func=self.get_write_func()

        if self.is_scalar():
            s="""
        (fits, %(colnum)s, row+1, firstelem, %(nelem)s, &str->%(colname)s, status);
            """ % {'write_func':write_func,
                   'colnum':colnum,
                   'nelem':self.n_elem,
                   'colname':self.f_name}
            s=s.strip()
        else:
            raise ValueError("implement non-scalar field")

        return s

    def get_read_statement(self, colnum):
        read_func=self.get_read_func()

        if self.is_scalar():
            s="""
        %(read_func)s(fits, %(colnum)s, row+1, firstelem, %(nelem)s, " ", &str->%(colname)s, NULL, status);
            """ % {'read_func':read_func,
                   'colnum':colnum,
                   'nelem':self.n_elem,
                   'colname':self.f_name}
            s=s.strip()
        else:
            raise ValueError("implement non-scalar field")

        return s


    def _set_field_info(self):
        field=self.field_def

        fs = field.split(' ')
        self.f_type = fs[0]
        f_name = fs[1]

        if self.f_type != 'char':
            raise ValueError("expected char field, got '%s'" % self.field_def)

        self.f_name, self.f_slen = _extract_string_name_len(f_name)

        # for now no arrays
        self.f_shape = (self.f_slen,)
        self.n_elem=1

        self.tform, self.fits_type, self.short_name, self.tdim = \
                c_string2fits(self.f_shape)
  
class NumberField(FieldBase):
    def __init__(self, field_def):
        super(NumberField,self).__init__(field_def)

        if '[' in field_def:
            raise ValueError("arrays not supported: '%s'" % field_def)

    def get_write_func(self):
        write_func='fits_write_col_%s' % self.short_name
        return write_func
    def get_read_func(self):
        write_func='fits_read_col_%s' % self.short_name
        return write_func


    def get_write_statement(self, colnum):
        write_func=self.get_write_func()

        if self.is_scalar():
            s="""
        %(write_func)s(fits, %(colnum)s, row+1, firstelem, %(nelem)s, &str->%(colname)s, status);
            """ % {'write_func':write_func,
                   'colnum':colnum,
                   'nelem':self.n_elem,
                   'colname':self.f_name}
            s=s.strip()
        else:
            raise ValueError("implement non-scalar field")

        return s

    def get_read_statement(self, colnum):
        read_func=self.get_read_func()

        if self.is_scalar():
            s="""
        %(read_func)s(fits, %(colnum)s, row+1, firstelem, %(nelem)s, 0, &str->%(colname)s, NULL, status);
            """ % {'read_func':read_func,
                   'colnum':colnum,
                   'nelem':self.n_elem,
                   'colname':self.f_name}
            s=s.strip()
        else:
            raise ValueError("implement non-scalar field")

        return s


    def _set_field_info(self):
        field=self.field_def

        fs = field.split(' ')
        self.f_type = fs[0]
        self.f_name = fs[1]

        if self.f_type == 'struct':
            raise ValueError("struct fields not supported")

        if self.f_type[0:4] == 'char':
            raise ValueError("use a StringField for char types: '%s'" % self.field_def)

        # for now scalar
        self.f_shape = (1,)
        self.n_elem = 1

        self.tform, self.fits_type, self.short_name, self.tdim = \
                c_num2fits(self.f_type, self.f_shape)


def _extract_string_name_len(f_name):
    ns=f_name.split('[')
    name=ns[0]
    ls=ns[1].split(']')
    slen = int( ls[0] )

    return name, slen

 
def _field_def_to_field(field_def):
    print field_def
    fd=field_def.strip()
    if fd[0:4] == 'char':
        return StringField(fd)
    else:
        return NumberField(fd)

def c_num2fits(f_type, shape):
    tdim=None

    if f_type not in _table_C2fits_tform:
        raise ValueError("unsupported type '%s'" % f_type)

    tform = _table_C2fits_tform[f_type]['char']
    fits_type = _table_C2fits_tform[f_type]['tstr']
    short_name = _table_C2fits_tform[f_type]['short_name']

    count=reduce(lambda x, y: x*y, shape)
    tform = '%d%s' % (count,tform)

    if len(shape) > 1:
        tdim = list(reversed(shape))
        tdim = [str(e) for e in tdim]
        tdim = '(' + ','.join(tdim)+')'
    return tform, fits_type, short_name, tdim

def c_string2fits(shape):

    tdim = None

    # get the size of each string
    string_size = shape[0]

    # now the dimensions
    if len(shape) == 1:
        tform = '%dA' % string_size
    else:
        sdims = shape[1:]
        count=reduce(lambda x, y: x*y, sdims)
        count = string_size*count
        tform = '%dA' % count

        tdim = list(reversed(sdims))
        tdim = [string_size_str] + [str(e) for e in tdim]
        tdim = '(' + ','.join(tdim)+')'

    return tform, 'TSTRING', 'str', tdim


# for TFORM and integer data type

_table_C2fits_tform = {'uint8_t':{'char':'B', 'tstr':'TBYTE','short_name':'byt'},
                       'int8_t':{'char':'S','tstr':'TSBYTE','short_name':'sbyt'}, # gets converted to unsigned

                       'S' :{'char':'A', 'tstr':'TSTRING','short_name':'str'},

                       'uint16_t':{'char':'U','tstr':'TUSHORT','short_name':'usht'}, # gets converted to signed
                       'int16_t':{'char':'I','tstr':'TSHORT','short_name':'sht'},

                       'uint32_t':{'char':'V','tstr':'TUINT','short_name':'uint'}, # gets converted to signed
                       'int32_t':{'char':'J','tstr':'TINT','short_name':'int'}, # gets converted to signed

                       # 'uint64_t':{'char':'K','tstr':'TULONG','short_name':'ulng'}, # not supported
                       'int64_t':{'char':'K','tstr':'TLONG','short_name':'lng'},

                       'float32':{'char':'E','tstr':'TFLOAT','short_name':'flt'},

                       'float64':{'char':'D','tstr':'TDOUBLE','short_name':'dbl'} }

_table_C2fits_tform['uint8'] =  _table_C2fits_tform['uint8_t']
_table_C2fits_tform['int8'] =  _table_C2fits_tform['int8_t']

_table_C2fits_tform['uint16'] =  _table_C2fits_tform['uint16_t']
_table_C2fits_tform['ushort'] =  _table_C2fits_tform['uint16_t']
_table_C2fits_tform['int16'] =  _table_C2fits_tform['int16_t']
_table_C2fits_tform['short'] =  _table_C2fits_tform['int16_t']

_table_C2fits_tform['uint32'] =  _table_C2fits_tform['uint32_t']
_table_C2fits_tform['uint'] =  _table_C2fits_tform['uint32_t'] # not portable
_table_C2fits_tform['int32'] =  _table_C2fits_tform['int32_t']
_table_C2fits_tform['int'] =  _table_C2fits_tform['int32_t']   # not portable

#_table_C2fits_tform['uint64'] =  _table_C2fits_tform['uint64_t']
_table_C2fits_tform['int64'] =  _table_C2fits_tform['int64_t']
_table_C2fits_tform['long'] =  _table_C2fits_tform['int64_t'] # not portable

_table_C2fits_tform['float'] =  _table_C2fits_tform['float32']
_table_C2fits_tform['double'] =  _table_C2fits_tform['float64']
