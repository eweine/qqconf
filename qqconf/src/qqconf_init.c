#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void jointlevel_twosided(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"jointlevel_twosided", (DL_FUNC) &jointlevel_twosided, 4},
    {NULL, NULL, 0}
};

void R_init_qqconf(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}