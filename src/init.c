#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void DescendMin(void *, void *, void *, void *, void *);
extern void FindEqualGreaterM(void *, void *, void *, void *, void *);
extern void RectUnique(void *, void *, void *, void *, void *, void *, void *);



/* .Call calls */
extern SEXP rowcolttests(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"DescendMin",          (DL_FUNC) &DescendMin,           5},
  {"FindEqualGreaterM",          (DL_FUNC) &FindEqualGreaterM,           5},
  {"RectUnique",          (DL_FUNC) &RectUnique,           7},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"rowcolttests",          (DL_FUNC) &rowcolttests,           5},
    {NULL, NULL, 0}
};

 void R_init_MicrobiomeAnalystR(DllInfo *dll)
 {
   R_registerRoutines(dll, CEntries, CallEntries,NULL,NULL);
   R_useDynamicSymbols(dll, FALSE);
 }

//void R_init_MetaboAnalystR(DllInfo *info) {
//  R_RegisterCCallable("MetaboAnalystR", "add",  (DL_FUNC) &CEntries);
//}