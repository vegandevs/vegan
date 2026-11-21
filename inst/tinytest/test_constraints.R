## Constraint permutations share the same path for all constrained ordinations.
library("vegan")
library("permute")
set.seed(42)
sizes <- c(2, 5, 3, 4)
dat <- data.frame(plot = rep(LETTERS[1:4], sizes),
                   A = factor(rep(c('a','b','a','b'),sizes)),
                   z = rnorm(sum(sizes)),row.names = paste0('s',seq_len(sum(sizes))))
Y <- matrix(rexp(nrow(dat)*5),nrow(dat),5,dimnames=list(rownames(dat),NULL))
colnames(Y) <- paste0("sp",seq_len(ncol(Y)))
D <- vegdist(Y,'bray')
h <- how(plots=Plots(~plot,type='free',variables=~A),data=dat,
          observed=TRUE,complete=TRUE)
draws <- allPerms(nrow(dat),h)
indices <- permutationIndices(draws)

fit_constraints <- function(method, data, partial=FALSE) {
    response <- if (method %in% c('rda','cca')) Y else D
    formula <- if (partial) response ~ A + Condition(z) else response ~ A
    switch(method,
           rda=rda(formula,data),cca=cca(formula,data),
           dbrda=dbrda(formula,data),capscale=capscale(formula,data),
           adonis2=adonis2(formula,data=data,permutations=0))
}
stat_constraints <- function(fit) {
    if (inherits(fit,'anova')) return(unname(fit$F[1L]))
    if (is.null(fit$CCA) || fit$CCA$qrank == 0L) return(0)
    unname(anova(fit,permutations=0)$F[1L])
}
for (method in c('rda','cca','dbrda','capscale','adonis2')) {
    for (partial in c(FALSE,TRUE)) {
        references <- vapply(seq_len(nrow(indices)),function(i) {
            changed <- dat
            changed$A <- dat$A[indices[i,]]
            stat_constraints(fit_constraints(method,changed,partial))
        },numeric(1))
        fit <- fit_constraints(method,dat,partial)
        if (method == 'adonis2') {
            formula <- if (partial) D ~ A + Condition(z) else D ~ A
            result <- suppressWarnings(suppressMessages(adonis2(formula,data=dat,
                                model='constraints',permutations=h)))
        } else {
            result <- suppressWarnings(suppressMessages(anova(fit,model='constraints',permutations=h)))
        }
        expect_equal(unname(attr(result,'F.perm')[,1]),references,tolerance=1e-8,
                     info=paste(method,'partial',partial,'explicit refits'))
        observed <- stat_constraints(fit)
        expect_equal(unname(result$F[1]),observed,tolerance=1e-8)
        expect_equal(unname(result$`Pr(>F)`[1]),mean(references >= observed-sqrt(.Machine$double.eps)))
    }
}
fit <- rda(Y~A+Condition(z),dat)
expect_message(permutest(fit, model = "constraints", permutations = draws),
               pattern = "conditioning variables")
expect_silent(suppressMessages(permutest(fit, model = "constraints",
                                         permutations = draws)))
## A message must remain a message when warnings are promoted to errors.
local({
    previous <- options(warn = 2)
    on.exit(options(previous))
    expect_silent(suppressMessages(permutest(fit, model = "constraints",
                                             permutations = draws)))
})

## Model matrices suffice even after the original formula environment is gone.
fit <- local({ metadata <- dat; response <- Y; rda(response~A,data=metadata) })
fit$terminfo <- NULL
expect_equal(suppressMessages(anova(fit,model='constraints',permutations=h))$F[1],
             stat_constraints(fit))

## Plain constraint row shuffles, including numeric-matrix fits and workers.
fit <- rda(Y,model.matrix(~A,dat)[,-1,drop=FALSE])
set.seed(7)
p <- shuffleSet(nrow(dat),nset=7,check=FALSE)
result <- anova(fit,model='constraints',permutations=p)
expect_equal(nrow(attr(result,'F.perm')),7L)
reference <- vapply(seq_len(nrow(p)),function(i) {
    Xi <- model.matrix(~A,dat)[p[i,],-1,drop=FALSE]
    stat_constraints(rda(Y,Xi))
},numeric(1))
expect_equal(unname(attr(result,'F.perm')[,1]),reference)
if (.Platform$OS.type == 'unix') {
    expect_equal(anova(fit,model='constraints',permutations=p,parallel=2)$`Pr(>F)`,
                 result$`Pr(>F)`)
}
cluster <- parallel::makeCluster(2)
parallel::clusterCall(cluster,function(paths) .libPaths(paths),.libPaths())
worker_result <- anova(fit,model='constraints',permutations=p,parallel=cluster)
parallel::stopCluster(cluster)
expect_equal(worker_result$`Pr(>F)`,result$`Pr(>F)`)

## Two-factor whole-plot design, preserving complete factor/numeric profiles.
sizes2 <- c(1,2,2,3,3,4,4,5)
dat2 <- data.frame(plot=rep(letters[1:8],sizes2),
                    A=factor(rep(rep(c('a','b'),each=4),sizes2)),
                    B=factor(rep(rep(c('x','y'),each=2,times=2),sizes2)),
                    x=rep(rep(c(1,2,3,4),each=2),sizes2))
set.seed(8)
Y2 <- matrix(rexp(nrow(dat2)*5),nrow(dat2),5)
D2 <- vegdist(Y2,'bray')
h2 <- how(plots=Plots(~plot,type='free',variables=~A+B+x),data=dat2)
expect_equal(numPerms(nrow(dat2),h2),2520)
perms2 <- shuffleSet(nrow(dat2),nset=12,control=h2,check=FALSE)
idx2 <- permutationIndices(perms2)
for (method in c('rda','cca','dbrda','capscale','adonis2')) {
    construct <- function(data, testing=FALSE) {
        response <- if(method %in% c('rda','cca')) Y2 else D2
        if (method=='adonis2') {
            adonis2(response~A*B+x,data=data,permutations=if(testing) perms2 else 0,
                    model=if(testing) 'constraints' else 'reduced')
        } else {
            fn <- get(method,asNamespace('vegan'))
            fit <- fn(response~A*B+x,data=data)
            if(testing) anova(fit,model='constraints',permutations=perms2) else fit
        }
    }
    reference <- vapply(seq_len(nrow(idx2)),function(i) {
        changed <- dat2
        changed[c('A','B','x')] <- dat2[idx2[i,],c('A','B','x')]
        stat_constraints(construct(changed))
    },numeric(1))
    expect_equal(unname(attr(construct(dat2,TRUE),'F.perm')[,1]),reference,
                 tolerance=1e-8,info=paste('factorial',method))
}

## Subsets and missing rows align the design before generating permutations.
dat3 <- dat; dat3$A[2] <- NA
fit3 <- rda(Y~A,data=dat3,na.action=na.exclude,subset=seq_len(nrow(dat3))!=7)
h3 <- how(plots=Plots(~plot,type='free',variables=~A),data=dat3)
result3 <- suppressMessages(anova(fit3,model='constraints',permutations=h3))
expect_equal(nrow(getPlots(attr(result3,'control'))$variables),12L)

## Predictors and responses can have different row names.
rownames(Y) <- paste0("response", seq_len(nrow(Y)))
fit4 <- rda(Y ~ A, data = dat3, na.action = na.exclude,
            subset = seq_len(nrow(dat3)) != 7L)
result4 <- anova(fit4, model = "constraints", permutations = h3)
expect_identical(rownames(getPlots(attr(result4, "control"))$variables),
                 rownames(dat3)[-c(2L, 7L)])
expect_equal(attr(result4, "F.perm"), attr(result3, "F.perm"))

## A precomputed whole-plot design must have the fitted row order.
fit <- rda(Y ~ A, data = dat)
rev.dat <- dat[rev(seq_len(nrow(dat))), ]
rev.h <- how(plots = Plots(~plot, type = "free", variables = ~A),
             data = rev.dat, complete = TRUE)
rev.draws <- allPerms(nrow(rev.dat), rev.h)
expect_error(anova(fit, model = "constraints", permutations = rev.draws),
             pattern = "fitted order")

## Identical groups yield ties and a p-value of one.
two <- data.frame(plot=rep(c('a','b'),c(2,5)),A=rep(c(0,1),c(2,5)))
Ytwo <- matrix(seq_len(7),7)
htwo <- how(plots=Plots(~plot,type='free',variables=~A),data=two,observed=TRUE,complete=TRUE)
expect_equal(suppressMessages(adonis2(dist(Ytwo)~A,data=two,model='constraints',permutations=htwo))$`Pr(>F)`[1],1)

## Unsupported null hypotheses are explicit, including the pending contrast API.
fit <- rda(Y~A,dat)
expect_error(anova(fit,model='constraints',by='terms',permutations=h),pattern='overall')
expect_error(anova(fit,model='constraints',test='parametric'),pattern='overall')
expect_error(anova(fit,fit,model='constraints',permutations=0),pattern='single fitted')
expect_error(suppressMessages(permutest(fit,model='constraints',first=TRUE,permutations=h)),pattern='overall')
expect_error(anova(fit,permutations=h),pattern='model = "constraints"')
expect_error(permutest(fit,permutations=draws),pattern='model = "constraints"')
expect_error(anova(fit,model='constraints',permustat='modified-F',permutations=0),pattern='pseudo-F')
missing_variables <- how(plots=Plots(~plot,type='free',variables=~A),data=dat2)
expect_error(suppressMessages(anova(rda(Y2~A*B,dat2),model='constraints',permutations=missing_variables)),
             pattern='distinguish all distinct')

## Complete ordinary how() designs also count the observed draw once.
smallY <- matrix(c(1, 4, 2, 3, 7, 5), 3)
smallX <- c(1, 2, 4)
smallfit <- rda(smallY ~ smallX)
allrows <- allPerms(3, how(observed = TRUE))
explicit <- anova(smallfit, model = "constraints", permutations = allrows)
generated <- suppressMessages(anova(smallfit, model = "constraints",
    permutations = how(observed = TRUE, complete = TRUE)))
expect_equal(generated$`Pr(>F)`, explicit$`Pr(>F)`)

## A condition can use degrees of freedom while explaining zero inertia.
condition <- rep(c(-1, 1), each = 4)
constraint <- c(1, 4, 6, 2, 3, 8, 9, 5)
response <- cbind(c(-3, -1, 1, 3, -3, -1, 1, 3),
                  c(5, 3, -3, -5, 1, -1, -3, 3))
fit <- rda(response ~ constraint + Condition(condition))
take <- c(8, 3, 6, 7, 1, 4, 2, 5)
draw <- matrix(take, nrow = 1)
refit <- rda(response ~ I(constraint[take]) + Condition(condition))
result <- suppressMessages(permutest(fit, permutations = draw, model = "constraints"))
expect_equal(fit$pCCA$rank, 0)
expect_equal(result$df, c(1, 5))
expect_equal(result$F.perm[1, 1], stat_constraints(refit))

## PERMANOVA uses raw sums of squares, even for very small dissimilarities.
smallD <- dist(response * 1e-5)
result <- adonis2(smallD ~ constraint, permutations = draw, model = "constraints")
refit <- adonis2(smallD ~ I(constraint[take]), permutations = 0)
expect_equal(result$F[1], adonis2(smallD ~ constraint, permutations = 0)$F[1])
expect_equal(attr(result, "F.perm")[1, 1], stat_constraints(refit))

## Losing all residual degrees of freedom gives the same error in workers.
response <- cbind(c(3, 7, 1, 5), c(1, 2, 8, 0))
condition <- c(1, 2, 4, 8)
constraint <- c(2, 7, 5, 3)
fit <- suppressMessages(rda(response ~ condition + constraint + Condition(condition)))
draw <- rbind(1:4, c(2, 1, 4, 3))
expect_error(suppressMessages(permutest(fit, permutations = draw, model = "constraints")),
             pattern = "no residual degrees of freedom")
if (.Platform$OS.type == "unix") {
    expect_error(suppressMessages(permutest(fit, permutations = draw,
                                           model = "constraints", parallel = 2)),
                 pattern = "no residual degrees of freedom")
}
