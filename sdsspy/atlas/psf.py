"""
Package:
    atlas
Module:
    psf

Purpose:
    Wrappers for Robert Lupton's PSF KL decomposition reader. Reads
    psField files and reconstruct PSF images.

Classes:
    PSFKL:
        Read PSF KL decompositions and reconstruct PSF images at a particlar
        row and column

"""

import numpy
import sdsspy
import _py_atlas


class PSFKL:
    """
    Class:
        PSFKL
    Purpose:
        Read PSF KL decomposition from a psField file.  Reconstruct the
        PSF at a given location.

    Construction:
        kl = PSFKL(psfield_filename, filter)
            filename: the name of a psField file.
            filter: should be a number 0-4 or letter u,g,r,i,z
                the file extension is the filter number plus one.

    Methods: 

        load:
            Load data from a psField file
        rec:
            reconstruct the PSF at a given location.  See docs for the .rec
            method for more info.

    Requires:
        The C extension for reading PSF KL decompositions.

    """
    def __init__(self, filename=None, filter=None):
        self.load(filename, filter)
    def load(self, filename, filter):
        """
        Load a PSF KL decomposition from a psField file.  See the docs for the
        main class PSFKL.
        """
        if filename is not None:
            if filter is None:
                raise ValueError("Send filename AND filter")
            if isinstance(filter, (str, unicode)):
                filternum = sdsspy.util.BANDS[filter]
            else:
                filternum = filter

            self.basis = _py_atlas.py_read_kl_basis(filename, filternum+1)

            self.filename=filename
            self.filter=filter
            self.filternum = filternum
        else:
            self.filename=None
            self.filter=None
            self.filternum=None
            self.basis = None

    def rec(self, rowc, colc, counts=1.0, ncomp=None, trim=False, more=False):
        """
        Class:
            PSFKL
        Method:
            rec
        Purpose:
            Reconstruct a PSF image a a given row and column.
        Calling Sequence:
            kl = PSFKL(filename, filter)
            image = kl.rec(row, col, counts=1.0, ncomp=, trim=False)
            dict = kl.rec(row, col, counts=1.0, ncomp=, trim=False, more=False)

        Inputs:
            row,col: 
                The row column to create the reconstruction.
        Optional Inputs:
            counts: 
                The total counts in the image, default 1.
            ncomp: The number of eigenimages to use.  Default is the
                number that exists in the file, usually 4.
            trim:
                If True, remove the boundary region.  Default False.
            more:
                If true, return a dict with 'image','coeffs','ecoeffs'
        Outputs:
            image: 
                A PSF image.  If more=True, a dictionary, see the more keyword.

        """
        if self.basis is None:
            raise ValueError("No basis loaded")

        nrow_b = self.basis['nrow_b']
        ncol_b = self.basis['ncol_b']
        nb = nrow_b*ncol_b
        coeffs = numpy.zeros(nb, dtype='f4')
        rcs = self.basis['RC_SCALE']

        if ncomp is None:
            ncomp = self.basis['ncomp']
        ecoeffs = numpy.zeros(ncomp, dtype='f4')

        for i in xrange(nb):
            coeffs[i]=(rowc*rcs)**(i % nrow_b) * (colc*rcs)**(i//nrow_b)

        for j in xrange(ncomp):
            cmat = self.basis['c'][j]
            #print "j=",j
            for i in xrange(nb):
                #print "    ",i % nrow_b, i//nrow_b
                ecoeffs[j] += cmat[i % nrow_b, i//nrow_b] * coeffs[i]

        image = ecoeffs[0]*self.basis['eigen'][0]
        for i in xrange(1,ncomp):
            image[:,:] += ecoeffs[i]*self.basis['eigen'][i][:,:]

        if trim:
            image = image[10:41, 10:41]

        '''
        image /= image.sum()
        if counts != 1.0:
            image *= counts
        '''

        if more:
            out={}
            out['coeffs'] = coeffs
            out['ecoeffs'] = ecoeffs
            out['image'] = image
            return out
        else:
            return image
