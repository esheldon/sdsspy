"""
Package:
    sdsspy
Modules:
    See these modules for more documentation.


    atlas:
        Read SDSS atlas images.
    files:
        A set of functions for dealing with SDSS files, including reading
        files, finding file lists, getting run metadata, converting SDSS
        id numbers to strings, etc.
    flags:
        Get SDSS flag values by name, e.g. target flags, etc.

    window:
        Tools for working with the SDSS window function, as well as getting
        lists of runs that have certain properties.

    yanny:
        tools for working with yanny parameter files.

    util:
        Utilities for working with SDSS data, such as flux to magnitude
        transformations, etc.

Modification History:
    Erin Sheldon, BNL
"""

import sys


try:
    import atlas
except:
    sys.stderr.write('sdsspy.atlas not loaded\n')


try:
    from . import files
    from files import filename
    from files import filedir
    from files import read
    from files import filespec
    from files import runlist
    from files import file_list
except:
    sys.stderr.write('sdsspy.files module not loaded\n')

try:
    from . import window
except:
    sys.stderr.write('sdsspy.window module not loaded\n')


try:
    import yanny
except:
    sys.stderr.write('sdsspy.yanny module not loaded\n')

try:
    import flags
    from flags import flagval
    from flags import Flags
except:
    sys.stderr.write('sdsspy.flags module not loaded\n')


try:
    import util
    from util import *
except:
    sys.stderr.write('sdsspy.util module not loaded\n')

