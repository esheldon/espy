"""
Module:
    flags

Convenience Functions:
    Name:
        flagval
    Calling Sequence:
        >>> import sdsspy
        >>> flagval = sdsspy.flags.flagval(flagprefix, label)

    This is a convenience function that instantias a Flags object and
    calculates the required flag info.

Classes:

    Flags

    Purpose:
        Work with SDSS flag values, e.g. photometric processing flags, target
        flags, etc.  The flag info is read from $ESPY_DIR/sdsspy/data/sdssMaskbits.par

        See the docs for the Flags class and it's val() function for more info.

    Calling Sequence:
        >>> import sdsspy
        >>> f=sdsspy.flags.Flags()
        >>> f.val('bosstile_status','previous_chunk')
        32768

    Other Useful Methods:
        load(reload=False)
        reload()
        filename()

Modification History:
    2010-04-07: Erin Sheldon, BNL
"""

import sdsspy
import numpy


def flagval(flagprefix, label):
    f=Flags()
    return f.val(flagprefix, label)

class Flags():
    """
    Class:
        Flags
    Purpose:
        Work with SDSS flag values, e.g. photometric processing flags, target
        flags, etc.  The flag info is read from $ESPY_DIR/sdsspy/data/sdssMaskbits.par

    Calling Sequence:
        >>> import sdsspy
        >>> f=sdsspy.flags.Flags()
        >>> f.val('bosstile_status','previous_chunk')
        32768

    Other Useful Methods:
        load(reload=False)
        reload()
        filename()

    Modification History:
        2010-04-07: Erin Sheldon, BNL
    """
    def __init__(self, reload=False):
        self.load(reload=reload)

    def load(self, reload=False):
        """
        This is a *class* attribute, so will not be loaded
        by each instance unless you use reload.
        """
        if not hasattr(Flags,'_flagstruct') or reload:
            fname = self.filename()
            y = sdsspy.yanny.Yanny(fname)

            Flags._flagstruct = y.read(one=True)
    def reload(self):
        self.load(reload=True)
        
    def filename(self):
        filename='$ESPY_DIR/sdsspy/data/sdssMaskbits.par'
        filename=sdsspy.files.expand_sdssvars(filename)
        return filename

    def val(self, flagprefix, label):
        """
        Class:
            Flags
        Method Name:
            val
        Purpose:
            Get the flag value for the named flag.

        Calling Sequence:
            >>> import sdsspy
            >>> f=sdsspy.flags.Flags()
            >>> flagval = f.val(flag_prefix, flag_label)

        Inputs:
            flag_prefix: 
                The flag prefix, e.g. OBJECT1 or BOSS_TARGET1.  The names are
                not case sensitive.  This is the overall name associated with
                a particilar bitmask.

            flag_label: 
                The label for a particilar bit in the bitmask named by
                flag_prefix.  e.g. for flag_prefix OBJECT1 one might use
                'BRIGHT' or 'SATURATED'.  This can be a sequence, in which case
                the flags or logically ored.

            Outputs:
                flag value:  
                    2**bit, or if a sequence of flag labels is entered, the
                    logical or or all, equivalent to 2**bit1 + 2*bit2 + ....

        Modification History:
            2010-04-07: Erin Sheldon, BNL
        """
        # get ref to this for simpler notation
        flags = Flags._flagstruct

        prefix=flagprefix.upper()

        wprefix, = numpy.where( flags['flag'] == prefix )
        if wprefix.size == 0:
            raise ValueError("No flag called '%s'" % prefix)

        # now get the labels to "or" together, make sure upper case
        if not isinstance(label, (tuple,list,numpy.ndarray)):
            label=[label]

        label = [l.upper() for l in label]

        # get the associated bits
        flagval = 0
        for l in label:
            w,=numpy.where(flags['label'][wprefix] == l)
            if w.size == 0:
                raise ValueError("Flag type '%s' has no flag label '%s'" % (prefix,l))
            
            bit = flags['bit'][wprefix[w[0]]]

            flagval |= 2**bit

        return flagval

