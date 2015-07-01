"read.cep" <-
  function (file, maxdata = 10000, positive = TRUE, trace = FALSE,
            force = FALSE) 
{
  if (!force) {
    stop("R may crash: if you want to try, save your session and use 'force=TRUE'")
  }
  if (is.loaded("_gfortran_ioparm"))
      warning("It seems that you have used gfortran: the input may be corrupted\n")
  ftypes <- c("free", "open", "condensed")
  file <- path.expand(file)
  if (trace) 
    cat("File", file, "\n")
  if (file.access(file, 4) < 0) {
    stop("file does not exist or is not readable")
  }
  on.exit(.Fortran("cepclose", PACKAGE = "vegan"))
  cep <- .Fortran("cephead", file = file, kind = integer(1), 
                  nitem = integer(1), nst = integer(1), fmt = character(1), 
                  PACKAGE = "vegan")
  if (cep$kind > 3) 
    stop("Unknown CEP file type")
  if (trace) {
    cat("looks like", ftypes[cep$kind], "CEP file,\n")
    cat("with", cep$nitem, "items per record")
    if (cep$kind == 1) 
      cat(" and", cep$nst, "records")
    cat(".\n")
  }
  switch(cep$kind,
         cd <- .Fortran("cepfree",
                        nitem = as.integer(cep$nitem), 
                        axdat = as.integer(maxdata),
                        nsp = integer(1),
                        nst = as.integer(cep$nst),
                        i = integer(maxdata),
                        j = integer(maxdata),
                        y = double(maxdata),
                        w = double(cep$nitem),
                        ier = integer(1),
                        PACKAGE = "vegan"), 
         cd <- .Fortran("cepopen",
                        fmt = as.character(cep$fmt), 
                        nitem = as.integer(cep$nitem),
                        maxdat = as.integer(maxdata), 
                        nsp = integer(1),
                        nst = integer(1),
                        i = integer(maxdata), 
                        j = integer(maxdata),
                        y = double(maxdata),
                        w = double(cep$nitem), 
                        ier = integer(1),
                        PACKAGE = "vegan"),
         cd <- .Fortran("cepcond",
                        fmt = as.character(cep$fmt),
                        nitem = as.integer(cep$nitem),
                        maxdat = as.integer(maxdata),
                        nsp = integer(1),
                        nst = integer(1),
                        i = integer(maxdata),
                        j = integer(maxdata),
                        y = double(maxdata),
                        w = double(cep$nitem),
                        iw = integer(cep$nitem),
                        ier = integer(1),
                        PACKAGE = "vegan"))
  if (cd$ier) {
    if (cd$ier == 1) 
      stop("too many non-zero entries: increase maxdata")
    else stop("unknown and obscure error: I do not know what to do")
  }
  if (trace) 
    cat("Read", cd$nsp, "species, ", cd$nst, "sites.\n")
  d <- matrix(0, cd$nst, cd$nsp)
  for (i in seq_along(cd$i)) d[cd$i[i], cd$j[i]] <- cd$y[i]
  nlines <- ceiling(cd$nsp/10)
  names <- NULL
  for (i in seq_len(nlines)) {
    tmpnames <- .Fortran("cepnames", character(1), PACKAGE = "vegan")
    tmpnames <- substring(as.character(tmpnames), 1, 80)
    tmpnames <- substring(tmpnames, seq(1, 80, by = 8), seq(8, 
                                                 80, by = 8))
    names <- c(names, tmpnames)
  }
  names <- gsub(" ", "", names)
  names <- make.names(names, unique = TRUE)
  colnames(d) <- names[seq_len(ncol(d))]
  nlines <- ceiling(cd$nst/10)
  names <- NULL
  for (i in seq_len(nlines)) {
    tmpnames <- .Fortran("cepnames", character(1), PACKAGE = "vegan")
    tmpnames <- substring(as.character(tmpnames), 1, 80)
    tmpnames <- substring(tmpnames, seq(1, 80, by = 8), seq(8, 
                                                 80, by = 8))
    names <- c(names, tmpnames)
  }
  names <- gsub(" ", "", names)
  names <- make.names(names, unique = TRUE)
  rownames(d) <- names[seq_len(nrow(d))]
  if (positive) {
    rsum <- apply(d, 1, sum)
    csum <- apply(d, 2, sum)
    d <- d[rsum > 0, csum > 0]
  }
  as.data.frame(d)
}
