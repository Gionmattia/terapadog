# R/terapadog.R

#'
#' This function reads RNA and RIBO count files, checks input data validity and
#'  merges them into a single numerical matrix (expression.data.
#'  It also prepares the metatadata needed by padog (exp_de).
#' @importFrom KEGGREST keggLink keggList
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @importFrom utils combn
#' @importFrom stats var quantile sd na.omit
#' @param esetm A matrix containing the counts from RNA and RIBO samples.
#' Rownames must be ensembl GENEIDs, while column names must be sample names.
#' Refer to prepareTerapadogData.R to prepare input data.
#' @param exp_de A dataframe containing information regarding the samples.
#' It has number of rows equal to the columns of esetm.
#' It has a formatted vocabulary, but can be obtained by running
#' prepareTerapadogData.R.
#' @param paired Logical. Specify is the study has a paired design or not.
#' If it does, be sure that the pairs are specified in the "Block" column of the
#' exp_de dataframe.
#' @param gslist A list of named character vectors. Each vector is named after a
#' KEGG pathway ID and each element within the vector is an ENSEMBL gene ID for
#' a gene part of said pathway.
#' @param organism A three letter string giving the name of the organism
#' supported by the "KEGGREST" package.
#' @param gs.names Character vector with the names of the gene sets.
#' If specified, must have the same length as gslist.
#' @param NI Number of iterations allowed to determine the gene set score
#' significance p-values.
#' @param Nmin The minimum size of gene sets to be included in the analysis.
#' @param verbose Logical. If true, shows number of iterations done.
#' @param parallel Logical. Allows for parallel execution to speed up.
#' @param dseed Optional initial seed for random number generator (integer).
#' @param ncr The number of CPU cores used when parallel set to TRUE.
#' Default is to use all CPU cores detected
#' @return A dataframe with the PADOG score for each pathway in exam.
#' @export
#'
#' @example
#'
terapadog2 <- function (esetm = NULL, exp_de = NULL, paired = FALSE,
                       gslist = "KEGGRESTpathway", organism = "hsa" ,
                       gs.names = NULL, NI = 1000, Nmin = 3, verbose = TRUE,
                       parallel = FALSE, dseed = NULL, ncr = NULL) {

  # validity checks of the data
  # Initial checks on the data (as PADOG would do). Some have been modified/removed to fit new kind of data.
  if (length(gslist) == 1 && gslist == "KEGGRESTpathway") {
    stopifnot(nchar(organism) == 3)
    res <- KEGGREST::keggLink("pathway", organism)
    a <- data.frame(path = gsub(paste0("path:", organism),"", res),
                    gns = gsub(paste0(organism, ":"),"", names(res)))
    gslist <- tapply(a$gns, a$path, function(x) {
      as.character(x)
    })
    gs.names <- KEGGREST::keggList("pathway", organism)[paste0(organism, names(gslist))]
    names(gs.names) <- names(gslist)
    stopifnot(length(gslist) >= 3)
    rm(res, a)
  }
  stopifnot(is.matrix(esetm))
  stopifnot(all(dim(esetm) > 4))
  stopifnot(mode(gslist) == "list")
  stopifnot(length(gslist) >= 3)

  if (!is.null(gs.names)) {
    stopifnot(length(gslist) == length(gs.names))
  }
  stopifnot(class(NI) == "numeric")
  stopifnot(NI > 5)
  stopifnot(sum(rownames(esetm) %in% as.character(unlist(gslist))) >
                10 & !any(duplicated(rownames(esetm))))

  # Initial operations on the gene set lists. Unchanged from PADOG.
  gf <- table(unlist(gslist))
  if (!(stats::var(gf) == 0)) {
    if (stats::quantile(gf, 0.99) > mean(gf) + 3 * stats::sd(gf)) {
      gf[gf > stats::quantile(gf, 0.99)] <- stats::quantile(gf, 0.99)
    }
    gff <- function(x) {
      1 + ((max(x) - x)/(max(x) - min(x)))^0.5
    }
    gf <- gff(gf)
  }
  else {
    fdfd <- unique(unlist(gslist))
    gf <- rep(1, length(fdfd))
    names(gf) <- fdfd
  }
  allGallP <- unique(unlist(gslist))

  restg <- setdiff(rownames(esetm), names(gf))
  appendd <- rep(1, length(restg))
  names(appendd) <- restg
  gf <- c(gf, appendd)
  stopifnot(all(!duplicated(rownames(esetm))))
  stopifnot(sum(rownames(esetm) %in% allGallP) > 10)
  if (verbose) {
    cat(paste0("Starting with ", length(gslist), " gene sets!"))
    cat("\n")
  }
  gslist <- gslist[unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x)) >= Nmin
  }))]
  gs.names <- gs.names[names(gslist)]
  stopifnot(length(gslist) >= 3)
  if (verbose) {
    cat(paste0("Analyzing ", length(gslist), " gene sets with ", Nmin, " or more genes!"))
    cat("\n")
  }

  if (!is.null(dseed))
    set.seed(dseed)

  # Here starts what was modified more heavily.

  # This section groups matching RNA and RIBO samples under a single index
  # (to keep the pair together during permutations).

  # Orders matching RIBO and RNA samples by the SampleName value (which is provided by the user when submitting data)
  exp_de_ordered <- exp_de[do.call(order, exp_de["SampleName"]),]

  # Retrieves the indexes (from exp_de) for each RNA and matching RIBO Sample.
  ordered_indexes <- rownames(exp_de_ordered)
  grouped_indexes <- as.data.frame(matrix(ordered_indexes, ncol = 2, byrow = TRUE))

  # Adds to the dataframe info on the Group for each couple of indices
  grouped_indexes$Group <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Group

  block <- NULL
  if (paired) {
    # Retrieves the info for paired samples from the dataframe
    grouped_indexes$Block <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Block
    block <- factor(grouped_indexes$Block)
  }

  # Extracts the group vector from the dataframe, to use it for permutations as PADOG does
  group <- grouped_indexes$Group

  # Global variables

  G <- factor(group) # creates a vector of factors ("c" or "d" in riboNAV)
  Glen <- length(G) # number of samples (since each factor is tied to a sample)
  tab <- table(G) # counts frequency for each factor
  idx <- which.min(tab) # finds factors with least number of samples assigned
  minG <- names(tab)[idx] # retrieve its name
  minGSZ <- tab[idx] # retrieves number of samples assigned to it
  bigG <- rep(setdiff(levels(G), minG), length(G))

  # topSigNum <- dim(esetm)[1]

  # Creates all the possible "group" label permutations for the given data.
  combFun <- function(gi, countn = TRUE) {
    g <- G[gi]
    tab <- table(g)
    if (countn) {
      minsz <- min(tab)
      ifelse(minsz > 10, -1, choose(length(g), minsz))
    }
    else {
      dup <- which(g == minG)
      cms <- utils::combn(length(g), tab[minG])
      del <- apply(cms, 2, setequal, dup)
      if (paired) {
        cms <- cms[, order(del, decreasing = TRUE), drop = FALSE]
        cms[] <- gi[c(cms)]
        cms
      }
      else {
        cms[, !del, drop = FALSE]
      }
    }
  }


  if (paired) {
    bct <- tapply(seq_along(G), block, combFun, simplify = TRUE)
    nperm <- ifelse(any(bct < 0), -1, prod(bct))
    if (nperm < 0 || nperm > NI) {
      btab <- tapply(seq_along(G), block, `[`, simplify = FALSE)
      bSamp <- function(gi) {
        g <- G[gi]
        tab <- table(g)
        bsz <- length(g)
        minsz <- tab[minG]
        cms <- do.call(cbind, replicate(NI, sample.int(bsz,
                                                       minsz), simplify = FALSE))
        cms[] <- gi[c(cms)]
        cms
      }
      combidx <- do.call(rbind, lapply(btab, bSamp))
    }
    else {
      bcomb <- tapply(seq_along(G), block, combFun, countn = FALSE,
                      simplify = FALSE)
      colb <- expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1,
                                                                , drop = FALSE]
      combidx <- mapply(function(x, y) x[, y, drop = FALSE],
                        bcomb, colb, SIMPLIFY = FALSE)
      combidx <- do.call(rbind, combidx)
    }
  }
  else {
    # Checks number of permutations that would be created. Must be >0 and < 1000
    nperm <- combFun(seq_along(G))
    if (nperm < 0 || nperm > NI) {
      combidx <- do.call(cbind, replicate(NI, sample.int(Glen,
                                                         minGSZ), simplify = FALSE))
    } else {
      # If nperm is acceptable, then permutation matrix combidx is generated (is the product of "cms" of combFun)
      combidx <- combFun(seq_along(G), countn = FALSE)
    }
  }

  NI <- ncol(combidx)
  if (verbose) {
    cat("# of permutations used:", NI, "\n")
  }

  deINgs <- intersect(rownames(esetm), unlist(gslist))
  gslistINesetm <- lapply(gslist, match, table = deINgs, nomatch = 0)
  MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI +
                              1)

  gsScoreFun <- function(G, block) {
    force(G) # Forces evaluation of G
    force(block) # Forces evaluation of block

    # This bit re-assigns the "group" label to the samples according to the combidx matrix, if past the first iteration

    if (ite > 1) {
      G <- bigG # Need to check the original PADOG for bigG
      G[combidx[, ite - 1]] <- minG
      G <- factor(G)

      # Unpacks grouped_indexes, so to have "group" labels for each sample (RNA count and RIBO count) in a vector
      for (i in seq_along(G)) {

        # Get the new label from G
        value <- as.character(G[i])

        # Updates exp_de according to the retrieved indexes (from grouped_indexes) with the new "group" label
        exp_de[grouped_indexes[i,]$V1, ]$Group <- value
        exp_de[grouped_indexes[i,]$V2, ]$Group <- value
      }
    }

    if (paired) {
      design_TE <- ~ Block + Group + SeqType+ Group:SeqType # Paired designs do not have batch effect correction!
    }
    else { # Add if statement, if there is batch, then do this, else do not
      design_TE <- ~ Group + SeqType+ Group:SeqType
    }

    # Setup the ddsMat object
    ddsMat <- DESeq2::DESeqDataSetFromMatrix(
      countData = esetm,
      colData = exp_de,
      design = design_TE
    )

    # Calculate results (without printing messages on console)
    ddsMat <- suppressMessages(DESeq2::DESeq(ddsMat))

    # Extract specific comparison of interest
    res <- DESeq2::results(ddsMat, name = "Groupd.SeqTypeRIBO")

    # originally, moderated t-values (abs value) were retrived and stored in a df.
    # Now it just extracts adjusted p-values.
    de <- res$padj
    names(de) <- rownames(res)

    # The scaling happens.
    degf <- scale(cbind(de, de * gf[names(de)]))
    rownames(degf) <- names(de)
    degf <- degf[deINgs, , drop = FALSE]

    sapply(gslistINesetm, function(z) {
      X <- stats::na.omit(degf[z, , drop = FALSE])
      colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
    })
  }

  if (parallel && requireNamespace("doParallel", quietly = TRUE) &&
      requireNamespace("parallel", quietly = TRUE)) {
    ncores <- parallel::detectCores()
    if (!is.null(ncr))
      ncores <- min(ncores, ncr)
    if (verbose) {
      clust <- parallel::makeCluster(ncores, outfile = "")
    }
    else {
      clust <- parallel::makeCluster(ncores)
    }
    doParallel::registerDoParallel(clust)
    tryCatch({
      parRes <- foreach::foreach(ite = 1:(NI + 1), .combine = "c",
                        .packages = "DESeq2") %dorng% { # original: .packages = "limma"
                          Sres <- gsScoreFun(G, block)
                          tmp <- list(t(Sres))
                          names(tmp) <- ite
                          if (verbose && (ite%%10 == 0)) {
                            cat(ite, "/", NI, "\n")
                          }
                          tmp
                        }
      parRes <- do.call(cbind, parRes[order(as.numeric(names(parRes)))])
      evenCol <- (1:ncol(parRes))%%2 == 0
      MSabsT[] <- parRes[, !evenCol]
      MSTop[] <- parRes[, evenCol]
      rm(parRes)
    }, finally = parallel::stopCluster(clust))
  }
  else {
    if (parallel)
      message("Execute in serial! Packages 'doParallel' and 'parallel'\n needed for parallelization!")
    for (ite in 1:(NI + 1)) {
      Sres <- gsScoreFun(G, block)
      MSabsT[, ite] <- Sres[1, ]
      MSTop[, ite] <- Sres[2, ]
      if (verbose && (ite%%10 == 0)) {
        cat(ite, "/", NI, "\n")
      }
    }
  }
  meanAbsT0 <- MSabsT[, 1]
  padog0 <- MSTop[, 1]
  MSabsT <- scale(MSabsT)
  MSTop <- scale(MSTop)

  # mff is what checks the real score against the iterative ones.
  # Originally, PADOG compares t-scores (the bigger, the more significant).
  # Since we are using p-value, the smaller the more significant, so the comparison sign was changed.

  mff <- function(x) {
    mean(x[-1] < x[1], na.rm = TRUE) # Original mean(x[-1] > x[1], na.rm = TRUE)
  }
  PSabsT <- apply(MSabsT, 1, mff)
  PSTop <- apply(MSTop, 1, mff)
  PSabsT[PSabsT == 0] <- 1/NI/100
  PSTop[PSTop == 0] <- 1/NI/100

  # Removed part for plots
  if (!is.null(gs.names)) {
    myn <- gs.names
  }
  else {
    myn <- names(gslist)
  }
  SIZE <- unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x))
  }))
  res <- data.frame(Name = myn, ID = names(gslist), Size = SIZE,
                    meanAbsT0, padog0, PmeanAbsT = PSabsT, Ppadog = PSTop,
                    stringsAsFactors = FALSE)
  ord <- order(res$Ppadog, -res$padog0)
  res <- res[ord, ]
  return(res)
}
