\name{comparisons}
\alias{comparisons}
\alias{comparisons.default}
\alias{comparisons.cca}
\alias{comparisons.vegan_scenarios}
\alias{print.vegan_comparisons}
\title{Comparisons for Constrained Ordination Models}
\description{
Construct scientifically meaningful comparisons among fitted scenarios and
test their additional constrained inertia with permutation pseudo-\eqn{F}
statistics. Interactions and other fitted terms are represented through the
complete fitted model matrix.
}
\usage{
comparisons(object, ...)

\method{comparisons}{cca}(object, variable, comparison = "pairwise",
    reference = NULL, newdata = NULL, by = NULL, joint = FALSE,
    permutations = how(nperm = 999), strata = NULL,
    parallel = getOption("mc.cores"), adjust = "none", wts = NULL, ...)

\method{comparisons}{vegan_scenarios}(object, comparison = "pairwise",
    reference = NULL, joint = FALSE, permutations = how(nperm = 999),
    strata = NULL, parallel = getOption("mc.cores"), adjust = "none", ...)

\method{print}{vegan_comparisons}(x, ...)
}
\arguments{
\item{object}{A formula-fitted constrained ordination object, or a scenario
table returned by \code{comparison_scenarios}.}
\item{variable}{The name of one factor whose levels define scenarios.}
\item{comparison}{A built-in comparison name, character equations, a named
list, or a numeric comparison vector or matrix.}
\item{reference}{The reference or control factor level. The first fitted
factor level is used when this is omitted, including for
\code{comparison = "trt.vs.ctrl"}.}
\item{newdata}{The reference population over which scenarios are averaged.
The default uses fitted analysis data, \code{"balanced"} constructs a
balanced grid, \code{"mean"} constructs one mean-or-mode profile (or one per
\code{by} group), and a data frame supplies an explicit population or setting.}
\item{by}{Optional model variables defining groups within which scenarios
and recipe comparisons are computed.}
\item{joint}{Test all resulting rows as one joint comparison.}
\item{permutations}{A permutation count, matrix, or control object returned
by \code{\link[permute]{how}}.}
\item{strata}{An optional factor or integer vector restricting permutations.}
\item{parallel}{Number of parallel workers or a predefined socket cluster.}
\item{adjust}{A multiple-testing method accepted by
\code{\link[stats]{p.adjust}}.}
\item{wts}{Optional non-negative scenario-averaging weights, or the name of
a weight column in \code{newdata}.}
\item{x}{A result returned by \code{comparisons}.}
\item{\dots}{Arguments passed to methods or to \code{print.data.frame}.}
}
\details{
The model method is a convenience interface: it first calls
\code{comparison_scenarios} and then compares those rows. Constructing the scenario
table explicitly is useful when the meaning and ordering of the rows should be
checked, especially for interactions.

Built-in specifications are \code{"pairwise"}, \code{"sequential"}
(also \code{"consecutive"}), \code{"reference"},
\code{"trt.vs.ctrl"}, \code{"deviation"} (also \code{"effect"}),
\code{"helmert"}, and \code{"poly"}. Reversed variants are unnecessary
because reversing all signs does not change an inertia-based pseudo-\eqn{F}
test. For reference and treatment-versus-control comparisons, omission of
\code{reference} selects the first fitted factor level.

Character equations may use factor-level names in the one-call interface,
for example \code{"HF - BF = 0"} or
\code{"(NM + SF) / 2 - BF = 0"}. Alternatively, first construct scenarios
and use their displayed row identifiers, such as \code{"b2 - b1 = 0"}.
The latter syntax also permits cross-group questions such as a
difference-in-differences. Names containing special characters must be
backticked. Only linear comparisons with a zero null are supported.

When \code{by} is present, built-in recipes and factor-level equations are
repeated within every group. An expression containing scenario identifiers
instead addresses exactly those rows, even when they belong to different
groups. The identifiers refer to scenario rows, never to encoded model
coefficients, previously computed comparisons, inertias, or pseudo-\eqn{F}
statistics.

Each requested comparison is tested separately by default. With
\code{joint = TRUE}, their compiled design-space rows are stacked, reduced to
their independent rank, and tested once. Separate tests share one validated
permutation matrix. The \code{parallel} argument is passed to each individual
\code{permutest.cca} call, so comparisons are processed sequentially while
permutations within a test may run in parallel.

Internally, the fitted-model space allowed under the comparison null is added
to original \code{Condition} terms, and the complementary tested space is
passed to vegan's reduced-model permutation machinery. Inertia and
pseudo-\eqn{F} use the geometry of the fitted method. Coefficients are only an
internal coordinate system and are neither reported nor tested individually.

See \code{vignette("comparisons-vegan")} for a user-oriented guide and a
technical appendix.
}
\value{
A data frame of class \code{"vegan_comparisons"} containing the comparison,
degrees of freedom, inertia, proportion of total inertia, pseudo-\eqn{F}, and
permutation p-value. Grouping variables and adjusted p-values are included
when requested. Attributes retain the scenario table, compiled design-space
rows, model reference, permutation metadata, and diagnostics.
}
\examples{
data(dune)
data(dune.env)
mod <- rda(dune ~ Management + A1, data = dune.env)

comparisons(mod, "Management", permutations = 99)
comparisons(mod, "Management", comparison = "trt.vs.ctrl",
            permutations = 99)
comparisons(mod, "Management",
            comparison = c("HF - BF = 0", "(NM + SF) / 2 - BF = 0"),
            permutations = 99)

sc <- comparison_scenarios(mod, "Management")
sc
comparisons(sc, comparison = "b2 - b1 = 0", permutations = 99)
comparisons(sc, comparison = "(b3 + b4) / 2 - b1 = 0",
            permutations = 99)

set.seed(42)
idat <- expand.grid(treatment = factor(c("control", "low", "high"),
                                      levels = c("control", "low", "high")),
                    habitat = factor(c("grass", "heath")), replicate = 1:3)
iY <- matrix(rpois(nrow(idat) * 4, 5), ncol = 4)
imod <- rda(iY ~ treatment * habitat, data = idat)
isc <- comparison_scenarios(imod, "treatment", by = "habitat")
isc
comparisons(isc, comparison = "pairwise", permutations = 99)
}
\seealso{\code{\link{comparison_scenarios}}, \code{\link{anova.cca}},
\code{\link{rda}}, \code{\link{cca}}, \code{\link{dbrda}},
\code{\link[permute]{how}}}
\keyword{models}
\keyword{multivariate}
