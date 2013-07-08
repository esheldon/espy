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


    def get_create_table_func(self):
        names=self.get_names_c()
        tforms=self.get_tforms_c()
        struct_name=self.get_name_norm()

        func="""
int %s_create_table(fitsfile* fits, LONGLONG nrows, char *extname, int *status)
{
    int nfields=%d;
    char *names=%s;
    char *tforms=%s;
    char *tunits=NULL;
    fits_create_tbl(fits, BINARY_TBL, nrows, nfields, names, tforms, tunits,
                    extname, status);
}\n""" % (struct_name,len(self.fields), names, tforms)
        return func

class FieldBase(object):
    def __init__(self, field_def):
        self.field_def=field_def

        self._set_field_info()

    def get_name(self):
        return self.f_name
    def get_type(self):
        return self.f_type
    def get_shape(self):
        return self.f_shape
    def get_tform(self):
        return self.tform
    def get_tdim(self):
        return self.tdim

class StringField(FieldBase):
    def __init__(self, field_def):
        super(StringField,self).__init__(field_def)

        if '][' in field_def:
            raise ValueError("arrays not supported: '%s'" % field_def)

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

        self.tform, self.tdim = c_string2fits(self.f_shape)
  
class NumberField(FieldBase):
    def __init__(self, field_def):
        super(NumberField,self).__init__(field_def)

        if '[' in field_def:
            raise ValueError("arrays not supported: '%s'" % field_def)

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

        self.tform,self.tdim = c_num2fits(self.f_type, self.f_shape)

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

    tform = _table_C2fits_tform[f_type]

    count=reduce(lambda x, y: x*y, shape)
    tform = '%d%s' % (count,tform)

    if len(shape) > 1:
        tdim = list(reversed(shape))
        tdim = [str(e) for e in tdim]
        tdim = '(' + ','.join(tdim)+')'
    return tform, tdim

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

    return tform, tdim


# for TFORM
_table_C2fits_tform = {'uint8':'B',
                       'uint8_t':'B',

                       'int8':'S', # gets converted to unsigned
                       'int8_t':'S', # gets converted to unsigned

                       'S' :'A',

                       'int16_t':'U', # gets converted to signed
                       'uint16_t':'U', # gets converted to signed

                       'int16':'I',
                       'int16_t':'I',

                       'uint32':'V', # gets converted to signed
                       'uint32_t':'V', # gets converted to signed

                       'int':'J',  # not portable
                       'int32':'J',
                       'int32_t':'J',

                       'long':'K',  # not portable
                       'int64':'K',
                       'int64_t':'K',

                       'float':'E',
                       'float32':'E',

                       'double':'D',
                       'float64':'D'}


