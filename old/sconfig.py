import os
import sys

# don't copy these from python exec statements
#ignores = [os,sys]
ignores=[]

default_delim = '='
default_comment = '#'

class sconfig(dict):
    """
    Class: sconfig
        A class to load simple config files.  Inherits from dict

    Initialization:
        conf = sconfig(path=None,
                       python=False,
                       delim='=',
                       comment='#',
                       verbose=0)
        conf.load(path,
                  python=False,
                  delim='=',
                  comment='#',
                  verbose=0)


    examples:
        from sconfig import sconfig
        conf = sconfig('/path/to/config/file')
        conf.load('/another/config/file.py',python=True)
        conf.load('/and/another/config', delim=':', comment='!')

        data_dir = conf['data_dir']
        for k in conf:
            sys.stdout.write('%s = %s\\n' % (k,conf[k]))

    Config files can be either of these 2 basic forms:
        1 
            python scripts, in which case the local variables are loaded using
            the exec() statement.  All the limitations of exec apply.  If
            modules are loaded in the config file and these will appear in the
            output.  The modules os and sys are ignored since these will be
            commonly loaded.
        2 
            files formatted as keyword-delimeter-value, e.g. 
                x = 3 
                y = stuff
            or 
                # this is colon separated
                x: 3 
                y: 'hello world'
            or even
                x 3
                y hello world      # a comment 
                
            The default delimiter is '=' but this can be changed through the
            delim keyword.  The first appearance of the dlimiter on the line
            indicates separation from keyword and value.  Everything after the
            first occurance of delim becomes part of the value.  Whitespace is
            stripped from around the keyword and variable.  
            
            The keyword must be a valid key name for a dictionary.  A bare
            number value gets converted to the python equivalent, 3 as an int,
            3.0 as a float.  If strings are surrounded by single or double
            quotes such as 
                x = 'hello world ' 
            then everything within quotes gets copied straight to a python
            string.  If the string is bare such as
                x = hello world
            then the string is everything after the delimiter with white 
            space stripped.

            Comments in these files are by default '#' but this can be changed
            through the comment keyword option.


    """
    def __init__(self, path=None, python=False, delim=default_delim,
                 comment=default_comment, verbose=0):
        if path is not None:
            self.load(path, python=python, delim=delim, 
                      verbose=verbose,comment=comment)

    def load_python(self, path, verbose=0):
        """
        Load the python config file using an exec() statement.
        """

        glob={}
        loc={}

        exec(open(path),glob,loc)

        for key in loc:
            if loc[key] not in ignores:
                self[key] = loc[key]
                if verbose > 1:
                    sys.stdout.write('    %s = %s\n' % (key,self[key]))


    def load_delim(self, path, delim, comment=default_comment, verbose=0):
        """
        Load the python config file using an exec() statement.
        """
        
        fobj = open(path,'r')

        for l in fobj:
            l = l.strip()

            lc = l.split(comment)
            lc = lc[0]
            
            ls = lc.split(delim)
            if len(ls) >= 2:
                varstring = ls[0].strip()
                rest = delim.join(ls[1:]).strip()

                # now variable assignment using eval.  First try direct
                # assignment, which should work for number, lists, or strings
                # marked with single quites 'hello' or double quotes "hello".
                # If the string is bare such as hello then we will get an
                # exception in most cases, in which case we will try gain
                # wrapping it in quotations

                try:
                    pcomm = "%s = %s" % (varstring, rest)
                    self[varstring] = eval(rest)
                except (NameError,SyntaxError):
                    pcomm = "%s = '%s'" % (varstring, rest)
                    self[varstring] = rest

                if verbose >= 2:
                    sys.stdout.write(pcomm+'\n')

    def load(self, path, python=False, delim=default_delim,
             comment=default_comment, verbose=0):
        """
        load(config_file, python=False, delim='=', comment='#', verbose=0)

        Load the config file.  
        """

        if verbose > 0: 
            sys.stdout.write("Reading from config file: %s\n" % path)
        
        if python:
            self.load_python(path, verbose=verbose)
        else:
            self.load_delim(path, delim, comment=comment, verbose=verbose)


