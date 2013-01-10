"""
Package:
    atlas
Purpose:

    Read SDSS atlas images and PSF KL decompositions.  This is a wrapper for
    Robert Luptons C code.

Modules:
    atlas:
        wrappers for the C atlas reader.  See docs for sdsspy.atlas.atlas for
        more info.
    psf:
        Contains a wrapper class PSFKL for the C PSF reader, and reconstructing
        the PSF at a row and column.  See docs for sdsspy.atlas.psf and
        sdsspy.atlas.psf.PSFKL for more info.

atlas.read_atlas and psf.PSFKL get pulled into this namespace as well, so you can access
them with
    sdsspy.atlas.read_atlas
    sdsspy.atlas.PSFKL

"""

import atlas
from atlas import read_atlas
from atlas import NoAtlasImageError
import psf
from psf import PSFKL
