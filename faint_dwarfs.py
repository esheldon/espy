import es_sdsspy
import esutil as eu
import cosmology

def nasa_sdss_file():
    return 'hdfs:///user/esheldon/boss/faint-dwarfs/nsa_v0_1_1_rwf.fits'
def nasa_sdss_edge_file():
    return 'hdfs:///user/esheldon/boss/faint-dwarfs/nsa_v0_1_1_rwf_edge.fits'
def read_nasa_sdss():
    f=nasa_sdss_file()
    data = eu.io.read(f,lower=True,verbose=True)
    return data
def write_nasa_sdss_edge(data, header=None):
    f=nasa_sdss_edge_file()
    eu.io.write(f, data, clobber=True, header=header, verbose=True)

def make_nasa_sdss_edge(stomp_map='boss', map_type='basic'):
    radii = [0.15, 0.25]
    data=read_nasa_sdss()

    out = eu.numpy_util.add_fields(data, [('basic_maskflags','i2',len(radii))])

    m = es_sdsspy.stomp_maps.load(stomp_map,map_type)
    cosmo = cosmology.Cosmo()


    header=[{'name':'map',     'value':stomp_map,        'comment':'stompe map'},
            {'name':'map_type','value':map_type,         'comment':'stomp map type'},
            {'name':'H0',      'value':cosmo.H0(),       'comment':'H0 used for distance calculations'},
            {'name':'flat',    'value':cosmo.flat(),     'comment':'Was a flat universe used?'},
            {'name':'omega_m', 'value':cosmo.omega_m(),  'comment':'omega_m used for distance calculations'},
            {'name':'omega_l', 'value':cosmo.omega_l(),  'comment':'omega_l used for distance calculations'},
            {'name':'omega_k', 'value':cosmo.omega_k(),  'comment':'omega_k used for distance calculations'}]

    write_nasa_sdss_edge(out, header=header)
