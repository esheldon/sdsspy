#include <Python.h>
#include "numpy/arrayobject.h"

#include "dervish.h"
#include "phFits.h"
#include "phConsts.h"

#include "py_atlas.h"

PyObject* py_read_atlas(const char* filename, int id) {
    ATLAS_IMAGE* ai=NULL;
    PyObject* dict=NULL;

    import_array();

    // This actually reads the entire object into memory
    // So getting mutliple bands should not be much overhead
    
    ai = py_load_atlas_images(filename, id);
    if (ai == NULL) {
        // you must return NULL to raise the exception
        return NULL;
    }
    dict = py_store_atlas_data(ai);
    phAtlasImageDel(ai,1);

    PyDict_SetItemString(dict, "filename", PyString_FromString(filename));
    PyDict_SetItemString(dict, "SOFT_BIAS", PyInt_FromLong(SOFT_BIAS));

    return dict;
}

/* Load the atlas images into an ATLAS_IMAGE struct */
ATLAS_IMAGE* py_load_atlas_images(const char* filename, int id) {
    FITS *fits=NULL;
    ATLAS_IMAGE* ai=NULL;

    if((fits = open_fits_table(filename, 1)) == NULL) {
        /* Memory was cleaned up in open_fits_table() */
        PyErr_Format(PyExc_IOError, "Error reading file: %s", filename);
        return NULL;
    }

    if((ai = read_atlas_image(fits,id)) == NULL) {
        phFitsDel(fits);  
        PyErr_Format(PyExc_IOError, "In file %s unknown error reading id %d", filename, id);
        return NULL;
    }
    phFitsDel(fits);  
  
    // no atlas image for this object
    if((ai)->id < 0) {			
        phAtlasImageDel(ai,1);
        PyErr_Format(PyExc_IOError, "In file %s object %d has no atlas image", filename, id);
        return NULL;
    }

    return ai;
}


PyObject* py_store_atlas_data(ATLAS_IMAGE* ai) {
    PyObject* dict=NULL;


    dict = PyDict_New();

    py_copy_images(ai, dict);
    py_copy_stats(ai, dict);


    return dict;
}

void py_copy_images(ATLAS_IMAGE* ai, PyObject* dict) {
    int band;
    PyObject* imlist;
    PyObject* imobj=NULL;

    imlist = PyList_New(0);
    for (band=0; band<NCOLOR; band++) {
        imobj = py_copy_image(ai, band);
        PyList_Append(imlist, imobj);
    }
    PyDict_SetItemString(dict, "images", imlist);
}

// this assumes uint16 for numpy array
PyObject* py_copy_image(ATLAS_IMAGE* ai, int band) {

    U16* pix;               /* pointer to pixels for specified band */
    npy_uint16* image;
    OBJMASK* om;
    SPAN *sp;				/* == om->s[i] */
    int nrow, ncol;
    int x1, x2, y;			/* unpacked from sp */
    int npix;				/* number of pixels in OBJMASK */
    int i,j;
    PyObject* imobj;        /* output image */


    pix = ai->pix[band];
    om = ai->master_mask;

    nrow = om->rmax - om->rmin + 1;
    ncol = om->cmax - om->cmin + 1;

    imobj = py_create_uint16_image(nrow, ncol);
    image = (npy_uint16* ) PyArray_DATA(imobj);

    // row0 = rmin + drow
    // then sent drow - row0 to ai_reg_set, which is just -rmin
    // this gets added to om->s[i].y, so we should just subtract rmin

    for (i=0; i<om->nspan; i++) {
        sp = &om->s[i];
        y = sp->y; x1 = sp->x1; x2 = sp->x2;
        npix = x2 - x1 + 1;

        y  -= om->rmin; 
        x1 -= om->cmin; 
        x2 -= om->cmin;

        if(y >= 0 && y < nrow) {
            int col0 = x1;			/* starting column in image */
            int ncopy = npix;		/* number of pixels to copy */
            int ptr0 = 0;			/* starting index in ptr[] */

            if(x1 < 0) {
                col0 = 0;
                ptr0 = -x1;

                ncopy -= ptr0;
            }
            if(x2 >= ncol) {
                ncopy -= (x2 - ncol + 1);
            }

            if(ncopy > 0) {
                for (j=0; j<ncopy; j++) {
                    // note y is really rows
                    //
                    // Here we copy to a signed integer first so we can subtract
                    // the bias, which could result in a negative number.
                    *(image + y*ncol + col0 + j)  = pix[ptr0 + j];
                }
            } else {
                printf("ncopy not > 0: %d\n", ncopy);
            }
        }
        pix += npix;
    }
    shAssert(pix == ai->pix[band] + om->npix);

    return imobj;
}

PyObject* py_create_uint16_image(int nrow, int ncol) {
    int ndims=2;
    npy_intp dims[2];
    npy_intp i=0, n=0;
    int fortran=0;

    PyObject* imobj;
    npy_uint16* data;

    dims[0] = nrow;
    dims[1] = ncol;
    imobj = PyArray_EMPTY(ndims, dims, NPY_UINT16, fortran);

    data = (npy_int16* ) PyArray_DATA(imobj);

    n=PyArray_SIZE(imobj);
    for (i=0; i<n; i++) {
        *data = SOFT_BIAS;
        data++;
    }

    return imobj;
}

void py_copy_stats(ATLAS_IMAGE* ai, PyObject* dict) {

    int band;
    PyObject* row0list=NULL;
    PyObject* col0list=NULL;
    PyObject* drowlist=NULL;
    PyObject* dcollist=NULL;

    int tmp;

    row0list = PyList_New(0);
    col0list = PyList_New(0);
    drowlist = PyList_New(0);
    dcollist = PyList_New(0);

    for (band=0; band<NCOLOR; band++) {
        tmp = ai->master_mask->rmin + ai->drow[band];
        PyList_Append(row0list, PyInt_FromLong(tmp));

        tmp = ai->master_mask->cmin + ai->dcol[band];
        PyList_Append(col0list, PyInt_FromLong(tmp));

        tmp = ai->drow[band];
        PyList_Append(drowlist, PyInt_FromLong(tmp));

        tmp = ai->dcol[band];
        PyList_Append(dcollist, PyInt_FromLong(tmp));

    }

    PyDict_SetItemString(dict, "row0", row0list);
    PyDict_SetItemString(dict, "col0", col0list);
    PyDict_SetItemString(dict, "drow", drowlist);
    PyDict_SetItemString(dict, "dcol", dcollist);

    PyDict_SetItemString(dict, "run", PyInt_FromLong(ai->run));
    PyDict_SetItemString(dict, "rerun", PyInt_FromLong(ai->rerun));
    PyDict_SetItemString(dict, "camcol", PyInt_FromLong(ai->camCol));
    PyDict_SetItemString(dict, "field", PyInt_FromLong(ai->field));
    PyDict_SetItemString(dict, "id", PyInt_FromLong(ai->id));
    PyDict_SetItemString(dict, "parent", PyInt_FromLong(ai->parent));

    PyDict_SetItemString(dict, "rmin", PyInt_FromLong(ai->master_mask->rmin));
    PyDict_SetItemString(dict, "rmax", PyInt_FromLong(ai->master_mask->rmax));
    PyDict_SetItemString(dict, "cmin", PyInt_FromLong(ai->master_mask->cmin));
    PyDict_SetItemString(dict, "cmax", PyInt_FromLong(ai->master_mask->cmax));


}
