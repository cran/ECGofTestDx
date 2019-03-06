#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void chebyshev_t_polynomials(void *, void *, void *);
extern void chebyshev_u_polynomials(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"chebyshev_t_polynomials", (DL_FUNC) &chebyshev_t_polynomials, 3},
    {"chebyshev_u_polynomials", (DL_FUNC) &chebyshev_u_polynomials, 3},
    {NULL, NULL, 0}
};

void R_init_ECGofTestDx(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

