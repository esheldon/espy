import subprocess
def convert(filename, type=None, dpi=None, crop=False, margins=None, noclobber=False,
            rotation=None, dosubdir=False, verbose=False, debug=False):

    """
    Call converter on the file.  
        def convert(filename, type=None, dpi=None, crop=False, margins=None, 
                    noclobber=True, rotation=None, dosubdir=False, 
                    verbose=False, debug=False)

    These mimic the converter options
    -d dpi    set the dpi for image conversions.  default 72
    -t otype  set the output type.  Supported types are all basic
            image types supported by convert.  Note for cropped
            pdf files the pdfcrop program is needed. default png
    -c        crop. default is no cropping
    -m        margins to use in the cropping.  Only PDF at this point
    -n        don't overwrite existing files.  Default is to overwrite
    -s        place outputs in a subdirectory otype/
    -r angle  Apply the rotation angle in degrees, counter-clockwise
    -v        be verbose
    -D        print debug info
    """

    command = ['converter']

    if type is not None:
        command += ['-t %s' % type]
    if dpi is not None:
        command += ['-d %s' % dpi]
    if crop:
        command += ['-c']
    if margins is not None:
        command += ['-m %s' % margins]
    if noclobber:
        command += ['-n']
    if dosubdir:
        command += ['-s']
    if rotation is not None:
        command += ['-r %s' % rotation]
    if verbose:
        command += ['-v']
    if debug:
        command += ['-D']

    command += [filename]
    subprocess.call(command)
