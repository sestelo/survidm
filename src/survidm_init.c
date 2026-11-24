#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void KaplanMeierValueSort(double *, int *, const int *, const double *, double *);
extern void LLWeightsKernel(const double *, const int *, const double *, const double *, const char *const *, double *);
extern void NWWeightsKernel(const double *, const int *, const double *, const double *, const char *const *, double *);
extern void SurvBeranKernel(double *, int *, double *, int *, const int *, const double *, const double *, const double *, const char *const *, double *);
extern void WeightsKaplanMeierSort(double *, int *, const int *, const double *, double *);
extern void WeightsKaplanMeierSortEx(double *, double *, const int *, const double *, double *);

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
