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
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(int *grid, int nx, int ny, int nz)}

/* Origin coordinates */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *reference, int ndims)}

/* PDB coordinates */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *atoms, int natoms, int xyzr)}

/* Sine and Cossine */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *sincos, int nvalues)}

%include "typemaps.i"

/* Map array of strings (Python) to char** (C) */
%typemap(in) 
char ** 
{
    /* Check if is a list */
    if (PyList_Check($input)) 
    {
        int size = PyList_Size($input);
        Py_ssize_t i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++)
        {
            PyObject *o = PyList_GetItem($input,i);
            if (PyUnicode_Check(o))
                $1[i] = PyUnicode_AsUTF8(PyList_GetItem($input,i));
            else 
            {
                //PyErr_SetString(PyExc_TypeError,"list must contain strings");
                PyErr_Format(PyExc_TypeError, "list must contain strings. %d/%d element was not string.", i, size);
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    } 
    else 
    {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

/* Free char ** array (C) */
%typemap(freearg) 
char ** 
{
  free((char *) $1);
}

/* Map char ** (C) to array of strings (Python) */
%typemap(out) 
(char **)
{
    int nC, nPy, py_err;
    PyObject *tmp;

    /* Define list length */
    for (nC = 0; $1[nC] != NULL; nC++);

    /* Create Python list */
    $result = PyList_New(nC);
    if (!$result) 
        return NULL;

    /* Pass C list to Python list */
    for (nPy = 0; nPy < nC; nPy++) 
    {
        if ($1[nPy] == NULL)
            break;

        /* Convert C string to Python string */
        tmp = PyString_FromString( $1[nPy] );
        if (!tmp) 
            return NULL;

        /* Pass Python string to Python list */
        py_err = PyList_SetItem($result, nPy, tmp);
        if (py_err == -1) 
            return NULL;
    }

   /* Delete C list */
   free($1);

   return $result; 
}

%include "SERD.h"
