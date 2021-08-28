%module SESD

%{
    #define SWIG_FILE_WITH_INIT
    #include "SESD.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%include "SESD.h"
