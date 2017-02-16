#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void abuswap(void *, void *, void *, void *, void *);
extern void curveball(void *, void *, void *, void *, void *);
extern void data2hill(void *, void *, void *, void *, void *, void *, void *, void *);
extern void dykstrapath(void *, void *, void *, void *, void *);
extern void goffactor(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pnpoly(void *, void *, void *, void *, void *, void *, void *);
extern void primtree(void *, void *, void *, void *, void *);
extern void quasiswap(void *, void *, void *, void *);
extern void rswapcount(void *, void *, void *, void *);
extern void stepabyss(void *, void *, void *, void *);
extern void stepacross(void *, void *, void *, void *);
extern void swap(void *, void *, void *, void *);
extern void swapcount(void *, void *, void *, void *);
extern void trialswap(void *, void *, void *, void *);
extern void veg_distance(void *, void *, void *, void *, void *, void *);
extern void wcentre(void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(cepclose)();
extern void F77_NAME(cepcond)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cepfree)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cephead)(void *, void *, void *, void *, void *);
extern void F77_NAME(cepnames)(void *);
extern void F77_NAME(cepopen)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cutup)(void *, void *, void *, void *);
extern void F77_NAME(eigy)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(monomds)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(orderdata)(void *, void *, void *, void *);
extern void F77_NAME(yxmult)(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"abuswap",      (DL_FUNC) &abuswap,       5},
    {"curveball",    (DL_FUNC) &curveball,     5},
    {"data2hill",    (DL_FUNC) &data2hill,     8},
    {"dykstrapath",  (DL_FUNC) &dykstrapath,   5},
    {"goffactor",    (DL_FUNC) &goffactor,    10},
    {"pnpoly",       (DL_FUNC) &pnpoly,        7},
    {"primtree",     (DL_FUNC) &primtree,      5},
    {"quasiswap",    (DL_FUNC) &quasiswap,     4},
    {"rswapcount",   (DL_FUNC) &rswapcount,    4},
    {"stepabyss",    (DL_FUNC) &stepabyss,     4},
    {"stepacross",   (DL_FUNC) &stepacross,    4},
    {"swap",         (DL_FUNC) &swap,          4},
    {"swapcount",    (DL_FUNC) &swapcount,     4},
    {"trialswap",    (DL_FUNC) &trialswap,     4},
    {"veg_distance", (DL_FUNC) &veg_distance,  6},
    {"wcentre",      (DL_FUNC) &wcentre,       4},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"cepclose",  (DL_FUNC) &F77_NAME(cepclose),   0},
    {"cepcond",   (DL_FUNC) &F77_NAME(cepcond),   11},
    {"cepfree",   (DL_FUNC) &F77_NAME(cepfree),    9},
    {"cephead",   (DL_FUNC) &F77_NAME(cephead),    5},
    {"cepnames",  (DL_FUNC) &F77_NAME(cepnames),   1},
    {"cepopen",   (DL_FUNC) &F77_NAME(cepopen),   10},
    {"cutup",     (DL_FUNC) &F77_NAME(cutup),      4},
    {"eigy",      (DL_FUNC) &F77_NAME(eigy),      27},
    {"monomds",   (DL_FUNC) &F77_NAME(monomds),   25},
    {"orderdata", (DL_FUNC) &F77_NAME(orderdata),  4},
    {"yxmult",    (DL_FUNC) &F77_NAME(yxmult),     9},
    {NULL, NULL, 0}
};

void R_init_vegan(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
