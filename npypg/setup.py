from distutils.core import setup, Extension
import numpy
import os
import sys

# A few extra include dirs
include_dirs=[numpy.get_include()]
lib_dirs=[]

if os.path.exists('/sw/lib'):
    lib_dirs.append('/sw/lib')

inc_trydirs = ['/usr/include/postgresql',
           '/usr/local/include/postgresql',
           '/usr/local/postgresql/include',
           '/sw/include/postgresql']
for trydir in inc_trydirs:
    if os.path.exists(trydir):
        include_dirs.append(trydir)


sw_altdir='/sw/opt'
if os.path.exists(sw_altdir):
    # look for a postgresql-* directory there
    l=os.listdir(sw_altdir)
    pg_incdirs=[]
    pg_libdirs=[]
    for d in l:
        if d.find('postgresql') != -1:
            pgdir=os.path.join(sw_altdir, d)

            incdir=os.path.join(pgdir, 'include')
            if os.path.exists(incdir):
                pg_incdirs.append(incdir)

            libdir=os.path.join(pgdir, 'lib')
            if os.path.exists(libdir):
                pg_libdirs.append(libdir)

    if len(pg_incdirs) > 0:
        pg_incdirs.sort()
        include_dirs.append(pg_incdirs[-1])

    if len(pg_libdirs) > 0:
        pg_libdirs.sort()
        lib_dirs.append(pg_libdirs[-1])


module1 = Extension('npypg.PGConn', 
                    sources=['npypg/src/PGConn.cpp'],
                    include_dirs=include_dirs,
                    library_dirs=lib_dirs,
                    libraries=['pq'])


# create the ups table
pyvers='%s.%s' % sys.version_info[0:2]
d1='lib/python%s/site-packages' % pyvers
d2='lib64/python%s/site-packages' % pyvers

if not os.path.exists('ups'):
    os.mkdir('ups')
tablefile=open('ups/npypg.table','w')
tab="""
# The default version of this file will be overwritten on setup to include
# paths determined from the python version.  This is useful to have in place
# though so that dependencies can be checked *before* installation.  Currently
# there are no required dependencies, so this is somewhat moot.

setupOptional("python")
setupOptional("numpy")
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
""" % (d1,d2)
tablefile.write(tab)
tablefile.close()



setup(name='npypg',
        version='0.10',
        description='Numerical python interface to postgresql',
        packages=['npypg'],
        data_files=[('ups',['ups/npypg.table'])],
        ext_modules=[module1]
    )
