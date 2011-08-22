#include <Python.h>
#include "numpy/arrayobject.h"

#include "dervish.h"
#include "phVariablePsf.h"
#include "phFits.h"
#include "phConsts.h"

#include "py_atlas.h"

//#include "py_psf.h"

PyObject* py_read_kl_basis(const char* filename, int ext) {

    FITS *fits=NULL;
    PSF_BASIS *basis = NULL;		/* PSF_BASIS derived from kl_comps */
    PyObject* result;

    import_array();

    fits = py_load_psField(filename, ext);
    if ( (fits = py_load_psField(filename, ext)) == NULL) {
        return NULL;
    }

    if ( (basis = py_load_kl_basis(fits)) == NULL) {
        phFitsDel(fits);  
        return NULL;
    }

    result = py_copy_eigen(basis);

    phFitsDel(fits);  
    phPsfBasisDel(basis);

    //result = PyInt_FromLong(-1);
    return result;
}

FITS* py_load_psField(const char* filename, int ext) {
    FITS *fits;
    if((fits = open_fits_table(filename, ext)) == NULL) {
        PyErr_Format(PyExc_IOError, "Error opening ext %d from file: %s", ext, filename);
        return NULL;
    }

    return fits;
}

PSF_BASIS* py_load_kl_basis(FITS* fits) {
    PSF_KL_COMP *klc=NULL;			/* PSF_KL_COMP read from file */
    PSF_BASIS *basis = NULL;		/* PSF_BASIS derived from kl_comps */
    int row=0;

    for(row = 1; row <= fits->naxis2; row++) {
        if((klc = read_KLC(fits, row)) == NULL) {
            PyErr_Format(PyExc_IOError, "Error reading row %d", row);
            return NULL;
        }

        /*
         * Create the PSF_BASIS; we need the first PSF_KL_BASIS to do this
         */
        if(basis == NULL) {
            PSF_KERNEL *kern = phPsfKernelNew(fits->naxis2 - 1, 1, 1,
                    klc->nrow_b, klc->ncol_b);
            basis = phPsfBasisNew(kern, KL_BASIS, -1, -1, klc->reg->type);
            phPsfKernelDel(kern);		/* decrement reference counter */
        }

        phPsfBasisSetFromKLComp(basis, row - 1, klc, 0);

        klc->reg = NULL; phPsfKLCompDel(klc);
    }

    return basis;

}

PyObject* py_copy_eigen(PSF_BASIS* basis) {
    const PSF_KERNEL *kern;		/* == basis->kern */
    int i;
    int ncomp;
    int nrow, ncol;


    PyObject* dict;
    PyObject* clist;
    PyObject* eigenlist;
    PyObject* imobj;

    kern = basis->kern;

    // this is the number of eigen images
    ncomp = kern->nsigma + 1;


    dict = PyDict_New();

    PyDict_SetItemString(dict, "nrow_b", PyInt_FromLong(kern->nrow_b));
    PyDict_SetItemString(dict, "ncol_b", PyInt_FromLong(kern->ncol_b));
    PyDict_SetItemString(dict, "ncomp", PyInt_FromLong(ncomp));
    PyDict_SetItemString(dict, "MAX_ORDER_B", PyInt_FromLong(MAX_ORDER_B));
    PyDict_SetItemString(dict, "RC_SCALE", PyFloat_FromDouble(RC_SCALE));

    //
    // the "c" matrix, holding the coefficients
    //

    clist = PyList_New(0);
    for (i=0;i<ncomp;i++) {
        imobj = py_copy_c(kern, i);
        PyList_Append(clist, imobj);
    }

    PyDict_SetItemString(dict, "c", clist);

    eigenlist = PyList_New(0);
    for (i=0;i<ncomp;i++) {
        imobj = py_copy_eigen_image(basis, i);
        PyList_Append(eigenlist, imobj);
    }

    PyDict_SetItemString(dict, "eigen", eigenlist);



    return dict;
}

PyObject* py_copy_c(const PSF_KERNEL* kern, int icomp) {
    int row,col;
    int ndims=2;
    npy_intp dims[2] = {MAX_ORDER_B, MAX_ORDER_B};
    int fortran=0;
    PyObject* imobj;
    float* image;

    imobj = PyArray_EMPTY(ndims, dims, NPY_FLOAT32, fortran);
    image = PyArray_DATA(imobj);

    for (row=0; row<MAX_ORDER_B; row++) {
        for (col=0; col<MAX_ORDER_B; col++) {
            *(image + row*MAX_ORDER_B + col)  = kern->a[icomp-1][0][0].c[row][col];
        }
    }

    return imobj;
}

PyObject* py_copy_eigen_image(PSF_BASIS* basis, int icomp) {
    int ndims=2;
    npy_intp dims[2];
    int fortran=0;

    PyObject* imobj;
    float* image;

    int nrow, ncol;
    int row, col;
    const REGION *comp;			/* component of KL decomposition */
    FL32 **rows_comp, *row_comp;		/* == comp->rows, comp->rows[] */

    comp = basis->regs[icomp-1][0][0]->reg;
    rows_comp = comp->rows_fl32;

    nrow = comp->nrow; 
    ncol = comp->ncol;

    dims[0] = nrow;
    dims[1] = ncol;

    imobj = PyArray_EMPTY(ndims, dims, NPY_FLOAT32, fortran);
    image = PyArray_DATA(imobj);

    for(row = 0; row < nrow; row++) {
        row_comp = rows_comp[row];
        for(col = 0; col < ncol; col++) {
            *(image + row*ncol + col) = row_comp[col];
        }
    }

    return imobj;
}
