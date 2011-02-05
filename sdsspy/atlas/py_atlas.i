%module py_atlas

%{
#include "py_atlas.h"
%}

PyObject* py_read_atlas(char* filename, int id);
