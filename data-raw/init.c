#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _parasol_rcpp_hello();

static const R_CallMethodDef CallEntries[] = {
    {"_parasol_rcpp_hello", (DL_FUNC) &_parasol_rcpp_hello, 0},
    {NULL, NULL, 0}
};

void R_init_parasol(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
