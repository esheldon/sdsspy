"""
Package:
    atlas
Module:
    atlas

Purpose:
    Wrapper for Robert Lupton's atlas reader.

Functions:
    read_atlas: 
        Read atlas images and metadata.  Calls the C routine.  See docs for
        that function for more details.

"""

import numpy
import _py_atlas

def read_atlas(filename, id, trim=False):
    """
    Module:
        atlas
    Name:
        read_atlas
    Purpose:
        Read atlas images and metadata from SDSS fpAtlas files.

    Calling Sequence:

        from sdsspy import atlas
        atlas_dict=atlas.read(filename, id, trim=False)

        Note: Typically you will use the sdsspy.files.read('fpatlas', run, camcol, ...)
            to read these files.

    Inputs:
            filename: 
                the path to an fpAtlas file.
            id: 
                The id of the object in the field.

    Optional Inputs:
        trim=: 
            If true, the atlas images will be trimmed to the region with actual
            data.  For example if an object was deblended from a large complex
            region, most of the area in the altas image is just set to the
            default background value.  This region will be trimmed.

    Outputs:
        A dictonary.  The dictionary has the following keys
            'images': 
                A list with images from each bandpass. 0->u, 1->g, 2->r, 3->i, 4->i
                The images are 2-d numpy arrays of type uint16.
            'row0', 'col0': 
                The offset of the lower left hand corner into the original
                frame.  The center in the atlas image is then [rowc-row0,
                colc-col0].  This entry is a list, one for each bandpass.
            'drow','dcol': 
                The difference in center relative to the cononical bandpass,
                the r-band.
            'rmin','rmax','cmin','cmax': 
                The bounding box of the atlas image.
            'run','rerun','camcol','field','id': 
                The sdss file information.
            'filename': 
                The atlas image file name.
            'SOFT_BIAS':
                The value of the background in the image.

    Exceptions:
        IOError is raised is the file cannot be read.
        NoAtlasImageError is raised if the requested object does not
            have an atlas image.  This is a benign error which you
            will most likely want to catch.

    """

    id=int(id)
    try:
        imdict = _py_atlas.py_read_atlas(filename, id)
    except RuntimeError as e:
        # The RuntimeError happens when an object has no atlas image, this
        # is no big deal.  Re-raise as NoAtlasImageError
        raise NoAtlasImageError(e) 

    if trim:
        for band in xrange(5):
            im = imdict['images'][band]

            # there is no way that I know of to get the exact region.
            #
            # This is not exact, since even noisy spots could end up zero, but
            # at least we know we would be in the noise anyway

            w1,w2 = numpy.where(im != imdict['SOFT_BIAS'])
            if w1.size > 0:
                minrow = w1.min()
                maxrow = w1.max()
                mincol = w2.min()
                maxcol = w2.max()

                if ((minrow > 0) 
                        or (maxrow < (im.shape[0]-1))
                        or (mincol > 0)
                        or (maxcol < (im.shape[1]-1))):
                    im = im[minrow:maxrow+1, mincol:maxcol+1]
                    imdict['images'][band] = im


    return imdict


class NoAtlasImageError(Exception):
    def __init__(self, value):
        self.value = str(value)
    def __str__(self):
        return self.value

