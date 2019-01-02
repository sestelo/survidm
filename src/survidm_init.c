#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void KaplanMeierValueSort(void *, void *, void *, void *, void *);
extern void LLWeightsKernel(void *, void *, void *, void *, void *, void *);
extern void NWWeightsKernel(void *, void *, void *, void *, void *, void *);
extern void SurvBeranKernel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void WeightsKaplanMeierSort(void *, void *, void *, void *, void *);
extern void WeightsKaplanMeierSortEx(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"KaplanMeierValueSort",     (DL_FUNC) &KaplanMeierValueSort,      5},
    {"LLWeightsKernel",          (DL_FUNC) &LLWeightsKernel,           6},
    {"NWWeightsKernel",          (DL_FUNC) &NWWeightsKernel,           6},
    {"SurvBeranKernel",          (DL_FUNC) &SurvBeranKernel,          10},
    {"WeightsKaplanMeierSort",   (DL_FUNC) &WeightsKaplanMeierSort,    5},
    {"WeightsKaplanMeierSortEx", (DL_FUNC) &WeightsKaplanMeierSortEx,  5},
    {NULL, NULL, 0}
};

void R_init_survidm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
