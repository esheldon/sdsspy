#include "dervish.h"
#include "phFits.h"
#include "phConsts.h"
#include "phVariablePsf.h"

#ifndef _PY_ATLAS_H
#define _PY_ATLAS_H

// Atlas images
PyObject* py_read_atlas(const char* filename, int id);

ATLAS_IMAGE* py_load_atlas_images(const char* filename, int id);
PyObject* py_store_atlas_data(ATLAS_IMAGE* ai);
void py_copy_atlas_images(ATLAS_IMAGE* ai, PyObject* dict);
PyObject* py_copy_atlas_image(ATLAS_IMAGE* ai, int band);

PyObject* py_create_uint16_image(int nrow, int ncol);

void py_copy_atlas_stats(ATLAS_IMAGE* ai, PyObject* dict);


// PSF KL decompositions
PyObject* py_read_kl_basis(const char* filename, int ext);
FITS* py_load_psField(const char* filename, int ext);
PSF_BASIS* py_load_kl_basis(FITS* fits);

PyObject* py_copy_eigen(PSF_BASIS* basis);
PyObject* py_copy_c(const PSF_KERNEL* kern, int icomp);
PyObject* py_copy_eigen_image(PSF_BASIS* basis, int icomp);

#endif
