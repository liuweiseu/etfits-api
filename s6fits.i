%module s6fits
%include "pyabc.i"
%include "std_vector.i"
%{
#include <stdio.h>
#include <time.h>
#include "s6fits.h"
#include "fitsio.h"
%}

%typedef long time_t;
%include "s6fits.h"

%template(s6Vector) std::vector<s6hits_t>;


