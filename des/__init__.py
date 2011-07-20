# Don't edit these svn properties by hand
_property_headurl='$HeadURL: svn+ssh://howdy.physics.nyu.edu/usr/local/svn/esheldon/trunk/python/des/__init__.py $'

def GetWlVersion():
    import re

    thisname='/python/des/__init__.py'
    headstrip=re.compile( '^\$HeadURL: .*/trunk/' )
    stripped_headurl = headstrip.sub('trunk/', _property_headurl)
    tag = stripped_headurl.replace(' $','')
    tag = tag.replace(thisname,'')

    return tag

version=GetWlVersion()

import checksg

import util
import collate
import select
import dc4
import dc5b
import archive
