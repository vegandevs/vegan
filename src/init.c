#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_stepacross(void *, void *, void *, void *);
extern void dykstrapath(void *, void *, void *, void *, void *);
extern void pnpoly(void *, void *, void *, void *, void *, void *, void *);
extern void primtree(void *, void *, void *, void *, void *);
extern void stepabyss(void *, void *, void *, void *);

/* .Call calls */
extern SEXP do_abuswap(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_chaoterms(SEXP);
extern SEXP do_curveball(SEXP, SEXP, SEXP);
extern SEXP do_decorana(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_getF(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_goffactor(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_backtrack(SEXP, SEXP, SEXP);
extern SEXP do_minterms(SEXP);
extern SEXP do_rcfill(SEXP, SEXP, SEXP);
extern SEXP do_qswap(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_swap(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_vegdist(SEXP, SEXP);
extern SEXP do_wcentre(SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(cepclose)();
extern void F77_NAME(cepcond)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cepfree)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cephead)(void *, void *, void *, void *, void *);
extern void F77_NAME(cepnames)(void *);
extern void F77_NAME(cepopen)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(monomds)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(orderdata)(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_stepacross", (DL_FUNC) &C_stepacross, 4},
    {"dykstrapath",  (DL_FUNC) &dykstrapath,  5},
    {"pnpoly",       (DL_FUNC) &pnpoly,       7},
    {"primtree",     (DL_FUNC) &primtree,     5},
    {"stepabyss",    (DL_FUNC) &stepabyss,    4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"do_abuswap",   (DL_FUNC) &do_abuswap,   4},
    {"do_chaoterms", (DL_FUNC) &do_chaoterms, 1},
    {"do_curveball", (DL_FUNC) &do_curveball, 3},
    {"do_decorana",  (DL_FUNC) &do_decorana,  7},
    {"do_getF",      (DL_FUNC) &do_getF,      8},
    {"do_goffactor", (DL_FUNC) &do_goffactor, 4},
    {"do_backtrack", (DL_FUNC) &do_backtrack, 3},
    {"do_minterms",  (DL_FUNC) &do_minterms,  1},
    {"do_qswap",     (DL_FUNC) &do_qswap,     4},
    {"do_rcfill",    (DL_FUNC) &do_rcfill,    3},
    {"do_swap",      (DL_FUNC) &do_swap,      4},
    {"do_vegdist",   (DL_FUNC) &do_vegdist,   2},
    {"do_wcentre",   (DL_FUNC) &do_wcentre,   2},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"cepclose",  (DL_FUNC) &F77_NAME(cepclose),   0},
    {"cepcond",   (DL_FUNC) &F77_NAME(cepcond),   11},
    {"cepfree",   (DL_FUNC) &F77_NAME(cepfree),    9},
    {"cephead",   (DL_FUNC) &F77_NAME(cephead),    5},
    {"cepnames",  (DL_FUNC) &F77_NAME(cepnames),   1},
    {"cepopen",   (DL_FUNC) &F77_NAME(cepopen),   10},
    {"monomds",   (DL_FUNC) &F77_NAME(monomds),   25},
    {"orderdata", (DL_FUNC) &F77_NAME(orderdata),  4},
    {NULL, NULL, 0}
};

void R_init_vegan(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
