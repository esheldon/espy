"""
    %prog [options] gmix_run
"""
import sys, os
import gmix_sdss

from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option('-g','--groups',default=None,
                  help='groups for wq, csv. Default is unconstrained')
parser.add_option('-n','--notgroups',default=None,
                  help='not groups for wq, csv. Default is unconstrained')
parser.add_option('-p','--priority',default='med',
                  help='priority for queue')
parser.add_option('--vers',default='work',
                  help="version for espy and gmix image")

parser.add_option('--bycamcol',action='store_true',
                  help="generate wq scripts by camcol instead of by field")

_wq_template="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/%(vers)s
    module unload gmix_image && module load gmix_image/%(vers)s

    python $ESPY_DIR/gmix_sdss/bin/run-gmix.py %(gmix_run)s %(run)s %(camcol)s %(field)s

%(groups)s
%(notgroups)s
job_name: %(job_name)s
priority: %(priority)s

"""

_wq_camcol_template="""
command: |
    source ~/.bashrc
    module unload espy && module load espy/%(vers)s
    module unload gmix_image && module load gmix_image/%(vers)s

    python $ESPY_DIR/gmix_sdss/bin/run-gmix.py %(gmix_run)s %(run)s %(camcol)s

%(groups)s
%(notgroups)s
job_name: %(job_name)s
priority: %(priority)s

"""


class WQWriter(dict):
    def __init__(self):
        options, args = parser.parse_args(sys.argv[1:])
        if len(args) < 1:
            parser.print_help()
            sys.exit(1)

        gmix_run=args[0]
        conf=gmix_sdss.files.read_config(gmix_run)

        self.update(conf)
        self['gmix_run'] = conf['run']

        self['priority']=options.priority
        self['vers'] = options.vers
        self['bycamcol'] = options.bycamcol

        groups=''
        if options.groups is not None:
            groups = 'group: [%s]' % options.groups
        notgroups=''
        if options.notgroups is not None:
            notgroups = 'notgroup: [%s]' % options.notgroups

        self['groups']=groups
        self['notgroups']=notgroups

        self.flist=gmix_sdss.files.read_field_cache(gmix_run=gmix_run)

    def go(self):
        if self['bycamcol']:
            self.write_bycamcol()
        else:
            self.write_byfield()

    def write_byfield(self):

        flist=self.flist
        runold=None
        for i in xrange(flist.size):
            self['run']=flist['run'][i]
            self['camcol']=flist['camcol'][i]
            self['field']=flist['field'][i]

            self['job_name']='%s-%06d-%s-%04d' % (self['gmix_run'],
                                                  self['run'],
                                                  self['camcol'],
                                                  self['field'])

            url=gmix_sdss.files.get_wq_url(**self)


            if self['run'] != runold:
                print url
                runold=self['run']

            text=_wq_template % self
            
            self._make_output_dir(url)
            with open(url,'w') as fobj:
                fobj.write(text)
        
    def write_bycamcol(self):
        import numpy
        flist=self.flist
        runs=numpy.unique(flist['run'])

        for run in runs:
            w,=numpy.where(flist['run']==run)
            camcols=numpy.unique(flist['camcol'][w])

            for camcol in camcols:

                self['run']=run
                self['camcol']=camcol

                self['job_name']='%s-%06d-%s' % (self['gmix_run'],
                                                      self['run'],
                                                      self['camcol'])

                url=gmix_sdss.files.get_wq_url(**self)

                text=_wq_camcol_template % self
                
                self._make_output_dir(url)
                print url
                with open(url,'w') as fobj:
                    fobj.write(text)
     

    def _make_output_dir(self, url):
        d=os.path.dirname(url)
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass

def main():
    writer=WQWriter()
    writer.go()

main()
