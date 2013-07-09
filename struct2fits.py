"""
todo:
    how to get header file name?  Just give header and struct name
    and extract?

    write auto-generated to files
"""

class StructFits(object):
    """

    struct_def = '''
        struct particle {
            int x;
            double y[4];
            long z;
            char sfield[35];
            uint32_t uf;
        };
    '''

    Note the struct definition does not have to be complete; only
    read/write the elements you want

    types can be standard C number types, those from <stdint.h>
    or char[len].  But note size_t etc are not supported.  Strings
    cannot be arrays.

    sf = Structfits(struct_def, header_list)
    sf.write_header_file(filename)
    sf.write_c_file(filename)

    Then use in C like this
    struct_particle_create_table(fits, nrows, "extname", &status);
    struct_particle_write_fits(fits, particles, nrows, &status);
    particles=struct_particle_read_fits(fits, &nrows, &status);


    """
    def __init__(self, struct_def, headers):
        self.struct_def = struct_def

        if not isinstance(headers,(list,tuple)):
            headers = [headers]

        self.headers=headers
        self._process_struct_def()

    def get_headers(self):
        """
        List of headers to include
        """
        return self.headers

    def get_name_raw(self):
        """
        The raw name, e.g. 'struct mystruct'.  If a typedef was used,
        the raw name equals the normalized name
        """
        return self.struct_name_raw

    def get_name_norm(self):
        """
        The normalized name, e.g. 'struct_mystruct'.  If a typedef was used,
        the raw name equals the normalized name
        """
        return self.struct_name_norm

    def get_names_c(self):
        """
        Get names as a C string array value
        """
        names=['"'+f.get_name()+'"' for f in self.fields]
        names = '{' + ','.join(names)+'}'
        return names

    def get_tforms_c(self):
        """
        Get tforms as a C string array value
        """
        tforms=['"'+f.get_tform()+'"' for f in self.fields]
        tforms = '{' + ','.join(tforms)+'}'
        return tforms


    def _process_struct_def(self):
        """
        Get the struct name and field list
        """
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

    def get_prototypes(self):
        proto = [self.get_create_table_signature(),
                 self.get_write_table_signature(),
                 self.get_read_table_signature()]
        proto = ';\n'.join(proto)
        return proto

    def get_guard_name(self):
        gname="_%s_INCLUDE_HEADER_GUARD" % self.get_name_norm().upper()
        return gname

    def get_includes(self):
        headers = self.get_headers()
        headers = ['#include "%s"' % h for h in headers]
        headers = '\n'.join(headers)

        return headers

    def write_header_file(self, filename):
        """
        Write the full header text
        """
        with open(filename,'w') as fobj:
            text=self.get_header_text()
            fobj.write(text)

    def write_c_file(self, filename):
        """
        Write the full .c text
        """
        with open(filename,'w') as fobj:
            text=self.get_c_text()
            fobj.write(text)


    def get_header_text(self):
        """
        The full text for the .h file
        """
        prototypes=self.get_prototypes()

        includes = self.get_includes()

        gname = self.get_guard_name()

        hdr="""#ifndef %(gname)s
#define %(gname)s

#include <fitsio.h>

%(includes)s

%(prototypes)s

#endif\n"""

        hdr = hdr % {'gname':gname,
                     'includes':includes,
                     'prototypes':prototypes}

        return hdr

    def get_c_text(self):
        """
        The full text for the .c file
        """

        includes = self.get_includes()
        texts = self.get_function_texts()

        hdr="""
#include <stdlib.h>
#include <stdio.h>
#include <fitsio.h>

%(includes)s

%(texts)s
\n"""

        hdr = hdr % {'includes':includes,
                     'texts':texts}

        return hdr




    def get_create_table_signature(self):
        """
        This is the function signature for the table creator
        """

        s="""
int %(struct_name_norm)s_create_table(fitsfile* fits, LONGLONG nrows, const char *extname, int *status)
        """ % {'struct_name_norm':self.get_name_norm()}
        return s.strip()


    def get_create_table_func(self):
        """
        Get the full create table function as a string
        """
        names=self.get_names_c()
        tforms=self.get_tforms_c()
        struct_name_norm=self.get_name_norm()

        sig=self.get_create_table_signature()
        func="""
// This function was automatically generated
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
        """
        Get all the column write statements for a given row
        """
        s=[]
        for i,f in enumerate(self.fields):
            s.append(f.get_write_statement(i))

        return '\n        '.join(s)

    def get_write_table_signature(self):
        """
        function signature for the table writer
        """
        s="""
int %(struct_name_norm)s_write_fits(fitsfile* fits, %(struct_name_raw)s *self, LONGLONG nrows, int *status)
        """ % {'struct_name_raw':self.get_name_raw(),
               'struct_name_norm':self.get_name_norm()}
        return s.strip()


    def get_function_texts(self):
        texts = [self.get_create_table_func(),
                 self.get_write_table_func(),
                 self.get_read_table_func()]
        texts = '\n\n'.join(texts)
        return texts

    def get_write_table_func(self):
        """
        get the full write table function as a string
        """
        struct_name_raw=self.get_name_raw()
        struct_name_norm=self.get_name_norm()
        col_writes = self.get_col_write_statements()

        sig = self.get_write_table_signature()

        goto_name='_%s_fits_write_bail' % (struct_name_norm)

        # doesn't work for arrays
        func="""
// This function was automatically generated
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
        """
        Get all the read statements for a given row
        """
        s=[]
        for i,f in enumerate(self.fields):
            s.append(f.get_read_statement(i))

        return '\n        '.join(s)

    def get_read_table_signature(self):
        """
        The function signature for reading the table
        """

        s="""
%(struct_name_raw)s *%(struct_name_norm)s_read_fits(fitsfile* fits, long *nrows, int *status)
        """ % {'struct_name_raw':self.get_name_raw(),
               'struct_name_norm':self.get_name_norm()}
        return s.strip()

    def get_read_table_func(self):
        """
        get the full function for reading the table as a string
        """
        struct_name_raw=self.get_name_raw()
        struct_name_norm=self.get_name_norm()
        col_reads = self.get_col_read_statements()

        goto_name='_%s_fits_read_bail' % (struct_name_norm)

        sig = self.get_read_table_signature()

        # doesn't work for arrays
        func="""
// This function was automatically generated
%(sig)s
{
    LONGLONG row=0;
    LONGLONG firstelem=1;
    %(struct_name_raw)s *self=NULL, *str=NULL;

    if (fits_get_num_rows(fits, nrows, status)) {
        goto %(goto_name)s;
    }

    self=calloc(*nrows, sizeof(%(struct_name_raw)s));
    if (!self) {
        fprintf(stderr,"could not allocate %%ld %(struct_name_raw)s\\n",*nrows);
        exit(1);
    }

    for (row=0; row<*nrows; row++) {
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
        *nrows=-1;
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
        return self.f_shape is None

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
            raise ValueError("string array columns not supported: '%s'" % field_def)

    def get_write_func(self):
        return 'fits_write_col_str'
    def get_read_func(self):
        return 'fits_read_col_str'

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

        # for now no arrays
        self.f_name, self.f_slen = _extract_string_name_len(f_name)

        self.f_shape = None
        self.n_elem=1

        self.tform, self.fits_type, self.short_name, self.tdim = \
                c_string2fits(self.f_slen)
  
class NumberField(FieldBase):
    def __init__(self, field_def):
        super(NumberField,self).__init__(field_def)

    def get_write_func(self):
        write_func='fits_write_col_%s' % self.short_name
        return write_func

    def get_read_func(self):
        write_func='fits_read_col_%s' % self.short_name
        return write_func

    def get_data_readwrite_name(self):
        if self.is_scalar():
            data_name = '&str->%s' % self.f_name
        else:
            # not no &, since this is an array
            data_name = 'str->%s' % self.f_name
        return data_name

    def get_write_statement(self, colnum):
        write_func=self.get_write_func()
        data_name = self.get_data_readwrite_name()

        s="""
        %(write_func)s(fits, %(colnum)s, row+1, firstelem, %(nelem)s, %(data_name)s, status);
        """ % {'write_func':write_func,
               'colnum':colnum,
               'nelem':self.n_elem,
               'data_name':data_name}
        s=s.strip()

        return s

    def get_read_statement(self, colnum):
        read_func=self.get_read_func()
        data_name = self.get_data_readwrite_name()

        s="""
        %(read_func)s(fits, %(colnum)s, row+1, firstelem, %(nelem)s, 0, %(data_name)s, NULL, status);
            """ % {'read_func':read_func,
                   'colnum':colnum,
                   'nelem':self.n_elem,
                   'data_name':data_name}
        s=s.strip()
        return s

    def _set_field_info(self):
        field=self.field_def

        fs = field.split(' ')
        self.f_type = fs[0]
        f_name = fs[1]

        if self.f_type == 'struct':
            raise ValueError("struct fields not supported")

        if self.f_type[0:4] == 'char':
            raise ValueError("use a StringField for char types: '%s'" % self.field_def)

        self.f_name, self.tform, self.fits_type, self.short_name, self.tdim, self.f_shape, self.n_elem = \
                c_num2fits(self.f_type, f_name)


def _extract_string_name_len(f_name):
    ns=f_name.split('[')
    name=ns[0]
    ls=ns[1].split(']')
    slen = int( ls[0] )

    return name, slen

 
def _field_def_to_field(field_def):
    #print field_def
    fd=field_def.strip()
    if fd[0:4] == 'char':
        return StringField(fd)
    else:
        return NumberField(fd)

def extract_shape(name0):
    import re
    ii=name0.find('[')
    if ii==-1:
        return name0,None


    name = name0[0:ii]
    shape_str = name0[ii:]
    dimstrs = re.findall('\[[0-9]+\]', shape_str)
    dims = [int(d.replace('[','').replace(']','')) for d in dimstrs]

    return name,dims


def c_num2fits(f_type, name):
    name,shape = extract_shape(name)

    if f_type not in _table_C2fits_tform:
        raise ValueError("unsupported type '%s'" % f_type)

    tform = _table_C2fits_tform[f_type]['char']
    fits_type = _table_C2fits_tform[f_type]['tstr']
    short_name = _table_C2fits_tform[f_type]['short_name']

    tdim=None
    if shape is None:
        count=1
    else:
        count=reduce(lambda x, y: x*y, shape)

        if len(shape) > 1:
            tdim = list(reversed(shape))
            tdim = [str(e) for e in tdim]
            tdim = '(' + ','.join(tdim)+')'

    tform = '%d%s' % (count,tform)
    return name, tform, fits_type, short_name, tdim, shape, count

def c_string2fits(string_size):
    tdim = None

    tform = '%dA' % string_size
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
