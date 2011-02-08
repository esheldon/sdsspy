"""
Package:
    atlas
Purpose:

    Read SDSS atlas images and PSF KL decompositions.  This is a wrapper for
    Robert Luptons C code.

Modules:
    atlas:
        wrappers for the C atlas reader.
    psf:
        Contains a wrapper class PSFKL for the C PSF reader, and reconstructing
        the PSF at a row and column.

See docs for the read function for more details.
"""
import atlas

from atlas import read_atlas

import psf
from psf import PSFKL
