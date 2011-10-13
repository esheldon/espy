import es_sdsspy
import esutil as eu
from esutil.numpy_util import where1
import cosmology
import numpy

_radii = [0.15, 0.25]

def nasa_sdss_file():
    return 'hdfs:///user/esheldon/boss/faint-dwarfs/nsa_v0_1_1_rwf.fits'
def nasa_sdss_edge_file():
    return 'hdfs:///user/esheldon/boss/faint-dwarfs/nsa_v0_1_1_rwf_edge.fits'

def read_nasa_sdss():
    f=nasa_sdss_file()
    data = eu.io.read(f,lower=True,verbose=True)
    return data

def read_nasa_sdss_edge():
    f=nasa_sdss_edge_file()
    data = eu.io.read(f,lower=True,verbose=True)
    return data

def write_nasa_sdss_edge(data, header=None):
    f=nasa_sdss_edge_file()
    eu.io.write(f, data, clobber=True, header=header, verbose=True)

def make_nasa_sdss_edge(stomp_map='boss', map_type='good', radmax=10.):
    """
    radmax in degrees
    """
    data=read_nasa_sdss()

    flagname='boss%s_maskflags' % map_type
    edgename='boss%s_edgeok' % map_type
    out = eu.numpy_util.add_fields(data, [(flagname,'i2',len(_radii)),
                                          (edgename,'i2',len(_radii))])
    out[flagname] = -1
    out[edgename] = 0

    m = es_sdsspy.stomp_maps.load(stomp_map,map_type)
    cosmo = cosmology.Cosmo()

    for i,r in enumerate(_radii):
        w=where1(out['z'] > 0.0)

        d = cosmo.Da(0.0, data['z'][w])
        srad = r/d*180.0/numpy.pi

        w2=where1(srad < radmax)
        print '    Found',w2.size,'with search rad <',radmax,'degrees'

        srad=srad[w2]
        w=w[w2]

        print '    search rad range from:',srad.min(),srad.max(),'degrees'
        print '        Checking edges'
        maskflags=m.Contains(out['ra'][w], out['dec'][w],"eq",srad)

        good=es_sdsspy.stomp_maps.quad_check(maskflags, strict=True)

        out[flagname][w,i] = maskflags[:]
        out[edgename][w[good],i] = 1

    rstr = ' '.join( [str(r) for r in _radii] )

    header=[{'name':'map',     'value':stomp_map,        'comment':'stomp map'},
            {'name':'map_type','value':map_type,         'comment':'stomp map type'},
            {'name':'H0',      'value':cosmo.H0(),       'comment':'H0 used for distance calculations'},
            {'name':'flat',    'value':cosmo.flat(),     'comment':'Was a flat universe used?'},
            {'name':'omega_m', 'value':cosmo.omega_m(),  'comment':'omega_m used for distance calculations'},
            {'name':'omega_l', 'value':cosmo.omega_l(),  'comment':'omega_l used for distance calculations'},
            {'name':'omega_k', 'value':cosmo.omega_k(),  'comment':'omega_k used for distance calculations'},
            {'name':'radii',   'value':rstr,             'comment':'radii used for edge checks'}]

    write_nasa_sdss_edge(out, header=header)

def doplot(region=None, file=None):
    import biggles
    data = read_nasa_sdss_edge()
    
    lam,eta=eu.coords.eq2sdss(data['ra'],data['dec'])

    w0=where1( data['bossgood_edgeok'][:,0] == 1)
    w1=where1( data['bossgood_edgeok'][:,1] == 1)

    xlabel=r'$\lambda$'
    ylabel=r'$\eta$'
    plt=eu.plotting.bscatter(lam,eta, type='dot', xlabel=xlabel, ylabel=ylabel, show=False)
    plt=eu.plotting.bscatter(lam[w0],eta[w0], type='dot', color='red', show=False,plt=plt)
    plt=eu.plotting.bscatter(lam[w1],eta[w1], type='dot', color='blue', show=False,plt=plt)

    plt=es_sdsspy.stomp_maps.plot_boss_geometry(plt=plt, region=region, show=False)

    fake=eu.plotting.fake_filled_circles(['all','r=%0.2f' % _radii[0], 'r=%0.2f' % _radii[1]],
                                         colors=['black','red','blue'])

    # bug in biggles, ignores valign='bottom' for PlotKey
    valign='top'

    if region is None:
        k_yval=0.6
        l_yval=0.6
    elif region == 'sgc':
        k_yval=0.2
        l_yval=0.05
        valign='bottom'
    else:
        k_yval=0.95
        l_yval=0.95

    key=biggles.PlotKey(0.95,k_yval,fake,halign='right',valign=valign)
    plt.add(key)

    if region is not None:
        nt=biggles.PlotLabel(0.05,l_yval,region.upper(),halign='left',valign=valign)
        plt.add(nt)

    if file is None:
        plt.show()
    else:
        print 'Writing plot file:',file
        plt.write_eps(file)
