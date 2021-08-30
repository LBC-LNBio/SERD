%module SERD

%{
    #define SWIG_FILE_WITH_INIT
    #include "SERD.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

/* Solvent-exposed surface grid */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* grid, int size)}

/* Origin coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims)}

/* PDB coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

/* Sine and Cossine */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

%include "SERD.h"
