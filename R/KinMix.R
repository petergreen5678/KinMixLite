aca2gt <-
function (mixture, aca) 
{
    nm <- length(names(aca))
    allele1 <- allele2 <- rep(0, nm)
    for (j in 1:nm) {
        m <- names(aca)[j]
        alleles <- mixture$data[[m]]$allele
        w <- which(aca[[m]] > 0)
        if (length(w) == 1) 
            allele1[j] <- allele2[j] <- alleles[w]
        else {
            allele1[j] <- alleles[w[1]]
            allele2[j] <- alleles[w[2]]
        }
    }
    data.frame(marker = names(aca), allele1, allele2)
}
add.child.meiosis.nodes <-
function (mixture, aca, ind = 1) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            nC <- aca[[m]]
            na <- length(nC)
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            for (a in 1:na) {
                Cpa <- paste("Cp", a, sep = "")
                Fna <- paste("n_", ind, "_", a, sep = "")
                ga <- paste("g", a, sep = "")
                gam1 <- paste("g", a - 1, sep = "")
                Cma <- paste("Cm", a, sep = "")
                CmSa <- paste("CmS", a, sep = "")
                CmSam1 <- ifelse(a == 2, "Cm1", paste("CmS", 
                  a - 1, sep = ""))
                Cna <- paste("Cn", a, sep = "")
                add.node(d, Cpa, states = 0:1, subtype = "numbered")
                add.edge(d, Cpa, Fna)
                if (a == 1) {
                  tab <- get.table(d, Cpa)
                  tab$Freq <- dbinom(tab[, Cpa], 1, tab[, Fna]/2)
                }
                else {
                  add.edge(d, Cpa, gam1)
                  tab <- get.table(d, Cpa)
                  probs <- (tab[, Fna] == 1) * tab[, gam1]/2 + 
                    (tab[, Fna] == 2)
                  tab$Freq <- dbinom(tab[, Cpa], 1, probs)
                }
                set.table(d, Cpa, tab, type = "cpt")
                if (a < na) {
                  add.node(d, ga, states = 0:2, subtype = "numbered")
                  if (a == 1) {
                    add.edge(d, ga, Cpa)
                    add.edge(d, ga, Fna)
                    tab <- get.table(d, ga)
                    tab$Freq <- ((tab[, ga] == 2) * (tab[, Cpa] == 
                      0) + (tab[, ga] == 0) * (tab[, Cpa] == 
                      1)) * (tab[, Fna] >= 1) + (tab[, ga] == 
                      1) * (tab[, Fna] == 0)
                  }
                  else {
                    add.edge(d, ga, gam1)
                    add.edge(d, ga, Cpa)
                    add.edge(d, ga, Fna)
                    tab <- get.table(d, ga)
                    tab$Freq <- ((tab[, ga] == 2) * (tab[, Cpa] == 
                      0) + (tab[, ga] == 0) * (tab[, Cpa] == 
                      1)) * (tab[, Fna] >= 1) * (tab[, gam1] == 
                      1) + (tab[, ga] == tab[, gam1]) * ((tab[, 
                      Fna] == 0) + (tab[, Fna] >= 1) * (tab[, 
                      gam1] != 1))
                  }
                  set.table(d, ga, tab, type = "cpt")
                }
                add.node(d, Cma, states = 0:1, subtype = "numbered")
                if (a == 1) {
                  tab <- get.table(d, Cma)
                  tab$Freq <- dbinom(tab[, Cma], 1, q[1])
                }
                else {
                  add.edge(d, Cma, CmSam1)
                  tab <- get.table(d, Cma)
                  s <- sum(tail(q, -(a - 1)))
                  tab$Freq <- dbinom(tab[, Cma], 1 - tab[, CmSam1], 
                    q[a]/ifelse(s > 0, s, 1))
                }
                set.table(d, Cma, tab, type = "cpt")
                if (a > 1 && a < na) {
                  add.node(d, CmSa, states = 0:1, subtype = "numbered")
                  add.edge(d, CmSa, Cma)
                  add.edge(d, CmSa, CmSam1)
                  tab <- get.table(d, CmSa)
                  tab$Freq <- ifelse(tab[, CmSa] == pmax(tab[, 
                    Cma], tab[, CmSam1]), 1, 0)
                  set.table(d, CmSa, tab, type = "cpt")
                }
                add.node(d, Cna, states = 0:2, subtype = "numbered")
                add.edge(d, Cna, Cpa)
                add.edge(d, Cna, Cma)
                tab <- get.table(d, Cna)
                tab$Freq <- ifelse(tab[, Cna] == tab[, Cma] + 
                  tab[, Cpa], 1, 0)
                set.table(d, Cna, tab, type = "cpt")
            }
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
add.motherchild.likd.node <-
function (mixture, Cgt, Mgt, db, ind = 1) 
{
    revid <- list()
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Cgt$marker) {
            add.node(d, "Rlikd", subtype = "boolean")
            datam <- mixture$data[[m]]
            q <- datam$freq
            q <- q/sum(q)
            ca1 <- match(Cgt[Cgt$marker == m, ]$allele1, datam$allele)
            ca2 <- match(Cgt[Cgt$marker == m, ]$allele2, datam$allele)
            ma1 <- match(Mgt[Mgt$marker == m, ]$allele1, datam$allele)
            ma2 <- match(Mgt[Mgt$marker == m, ]$allele2, datam$allele)
            w <- match(ca1, c(ma1, ma2))
            if (is.na(w)) {
                w <- match(ca2, c(ma1, ma2))
                if (!is.na(w)) {
                  t <- ca1
                  ca1 <- ca2
                  ca2 <- t
                }
            }
            if (is.na(w)) {
                tab <- get.table(d, "Rlikd")
                tab$Freq <- c(1, 0)
            }
            else {
                if (w == 2) {
                  t <- ma1
                  ma1 <- ma2
                  ma2 <- t
                }
                if (ca2 == ma2 && ca1 != ca2) {
                  Fna <- paste("n_", ind, "_", ca1, sep = "")
                  Fnb <- paste("n_", ind, "_", ca2, sep = "")
                  add.edge(d, "Rlikd", Fna)
                  add.edge(d, "Rlikd", Fnb)
                  tab <- get.table(d, "Rlikd")
                  lrs <- 0.5 * (tab[, Fna] + tab[, Fnb])/(q[ca1] + 
                    q[ca2])
                  mlr <- max(lrs)
                  lrs <- lrs/mlr
                  odd <- 1 == (1:length(lrs))%%2
                  lrs[odd] <- 1 - lrs[odd]
                  tab$Freq <- lrs
                  revid[[m]] <- mlr
                }
                else {
                  Fna <- paste("n_", ind, "_", ca2, sep = "")
                  add.edge(d, "Rlikd", Fna)
                  tab <- get.table(d, "Rlikd")
                  lrs <- 0.5 * tab[, Fna]/q[ca2]
                  mlr <- max(lrs)
                  lrs <- lrs/mlr
                  odd <- 1 == (1:length(lrs))%%2
                  lrs[odd] <- 1 - lrs[odd]
                  tab$Freq <- lrs
                  revid[[m]] <- mlr
                }
            }
            set.table(d, "Rlikd", tab, type = "cpt")
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
    Revid <<- revid
}
add.relative.likd.node <-
function (mixture, aca, ind = 1) 
{
    revid <- list()
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            add.node(d, "Rlikd", subtype = "boolean")
            if (max(aca[[m]]) == 2) {
                a <- which(aca[[m]] == max(aca[[m]]))
                Fna <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, "Rlikd", Fna)
                tab <- get.table(d, "Rlikd")
                lrs <- 0.5 * tab[, Fna]/q[a]
                mlr <- max(lrs)
                lrs <- lrs/mlr
                odd <- 1 == (1:length(lrs))%%2
                lrs[odd] <- 1 - lrs[odd]
                tab$Freq <- lrs
                revid[[m]] <- mlr
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                a <- ab[1]
                b <- ab[2]
                Fna <- paste("n_", ind, "_", a, sep = "")
                Fnb <- paste("n_", ind, "_", b, sep = "")
                add.edge(d, "Rlikd", Fna)
                add.edge(d, "Rlikd", Fnb)
                tab <- get.table(d, "Rlikd")
                lrs <- 0.25 * (tab[, Fna]/q[a] + tab[, Fnb]/q[b])
                mlr <- max(lrs)
                lrs <- lrs/mlr
                odd <- 1 == (1:length(lrs))%%2
                lrs[odd] <- 1 - lrs[odd]
                tab$Freq <- lrs
                revid[[m]] <- mlr
            }
            set.table(d, "Rlikd", tab, type = "cpt")
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
    Revid <<- revid
}
addU <-
function (IBD, nU) 
{
    extra <- matrix(0, nrow(IBD$patt), 2 * nU)
    m <- apply(IBD$patt, 1, max)
    for (r in 1:nrow(IBD$patt)) extra[r, ] <- m[r] + (1:(2 * 
        nU))
    IBD$patt <- cbind(IBD$patt, extra)
    IBD
}
as.gt <-
function (res00, s) 
{
    z <- subset(res00, Sim == s)
    allele <- z$Allele
    allele[allele == "X"] <- 0
    allele[allele == "Y"] <- 1
    allele <- matrix(as.numeric(allele), ncol = 2, byrow = TRUE)
    data.frame(marker = unique(z$Marker), allele1 = allele[, 
        1], allele2 = allele[, 2])
}
binary <-
function (n, d = NA) 
{
    res <- NULL
    while (ifelse(is.na(d), n > 0, d != 0)) {
        d <- abs(d) - 1
        res <- c(n%%2, res)
        n <- n%/%2
    }
    res
}
checkpeaks <-
function (epg, db, fix = 0.003) 
{
    first <- TRUE
    if ("height" %in% names(epg)) {
        for (m in unique(epg$marker)) {
            a <- epg[epg$marker == m, "allele"]
            adb <- db[db$marker == m, "allele"]
            w <- !a %in% adb
            if (any(w)) {
                if (first) {
                  cat("peaks at alleles not in database\n")
                  first <- FALSE
                }
                cat(m, ":", a[w], ":", epg[epg$marker == m, "height"][w], 
                  "\n")
                if (fix > 0) {
                  db <- merge(db, data.frame(marker = m, allele = a[w], 
                    frequency = fix), all = T)
                  f <- db[db$marker == m, ]$frequency
                  db[db$marker == m, ]$frequency <- f/sum(f)
                }
            }
        }
    }
    else {
        for (m in unique(epg$marker)) {
            a1 <- epg[epg$marker == m, "allele1"]
            a2 <- epg[epg$marker == m, "allele2"]
            adb <- db[db$marker == m, "allele"]
            for (a in list(a1, a2)) {
                w <- !a %in% adb
                if (any(w)) {
                  if (first) {
                    cat("observed alleles not in database\n")
                    first <- FALSE
                  }
                  cat(m, ":", a[w], "\n")
                  if (fix > 0) 
                    db <- merge(db, data.frame(marker = m, allele = a[w], 
                      frequency = fix), all = T)
                }
            }
            f <- db[db$marker == m, ]$frequency
            db[db$marker == m, ]$frequency <- f/sum(f)
        }
    }
    invisible(db)
}
convertIBD <-
function (x) 
{
    if (missing(x)) {
        x <- list(pr = 1, patt = matrix(c(1, 2, 1, 3), 1, 4))
    }
    if (is.character(x)) 
        x <- switch(pmatch(x, c("sibs", "parent-child", "half-sibs", 
            "cousins", "half-cousins", "second-cousins", "double-first-cousins", 
            "quadruple-half-first-cousins", "3cousins-cyclic", 
            "3cousins-star", "trio")), c(0.25, 0.5, 0.25), c(0, 
            1, 0), c(0.5, 0.5, 0), c(0.75, 0.25, 0), c(0.875, 
            0.125, 0), c(0.9375, 0.0625, 0), c(0.5625, 0.375, 
            0.0625), c(17, 14, 1)/32, list(pr = c(0.015625, 0.046875, 
            0.046875, 0.140625, 0.046875, 0.140625, 0.140625, 
            0.421875), patt = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 
            2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3, 3, 
            3, 3, 3, 4, 4, 4, 4, 2, 2, 3, 4, 1, 1, 3, 5, 3, 4, 
            4, 5, 3, 5, 5, 6), 8, 6)), list(pr = c(0.0625, 0.1875, 
            0.1875, 0.1875, 0.375), patt = matrix(c(1, 1, 1, 
            1, 1, 2, 2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 
            1, 4, 1, 3, 5, 4, 5, 5, 5, 6), 5, 6)), list(patt = c(1, 
            2, 3, 4, 1, 3)))
    if (is.list(x)) {
        if (length(x) == 2 && class(x[[1]]) %in% c("ped", "pedList")) 
            IBD <- pedigreeIBD(x[[1]], x[[2]])
        else IBD <- x
    }
    else {
        if (is.null(dim(x)) && length(x) == 3) {
            kappa <- x
            patt <- matrix(c(1, 2, 3, 4, 1, 2, 3, 1, 1, 2, 2, 
                1), 3, 4, byrow = TRUE)
            wk <- which(kappa > 0)
            IBD <- list(pr = kappa[wk], patt = patt[wk, , drop = FALSE])
        }
        else if (is.null(dim(x)) && length(x) == 9) {
            Delta <- x
            patt <- matrix(c(1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 
                2, 1, 1, 2, 3, 1, 2, 1, 1, 1, 2, 3, 3, 1, 2, 
                1, 2, 1, 2, 1, 3, 1, 2, 3, 4), 9, 4, byrow = TRUE)
            wk <- which(Delta > 0)
            IBD <- list(pr = Delta[wk], patt = patt[wk, , drop = FALSE])
        }
        else stop("invalid non-list x")
    }
    if (is.null(dim(IBD$patt))) 
        IBD$patt <- matrix(IBD$patt, nrow = 1)
    if (is.null(IBD$pr)) 
        IBD$pr <- rep(1/nrow(IBD$patt), nrow(IBD$patt))
    structure(IBD,class='IBD')
}
dbetabinom <-
function (x, n, alpha, beta) 
{
    lx <- max(length(x), length(n), length(alpha), length(beta))
    x <- rep_len(x, len = lx)
    n <- rep_len(n, len = lx)
    alpha <- rep_len(alpha, len = lx)
    beta <- rep_len(beta, len = lx)
    options(warn = -1)
    res <- choose(n, x) * exp(lgamma(alpha + x) + lgamma(beta + 
        n - x) + lgamma(alpha + beta) - (lgamma(alpha + beta + 
        n) + lgamma(alpha) + lgamma(beta)))
    options(warn = 0)
    res[is.na(res)] <- (n == x)[is.na(res)]
    res
}
delete.DQnodes <-
function (mixture, which = "DQ") 
{
    for (m in mixture$markers) {
        d <- mixture$domains[[m]]
        z <- get.nodes(d)
        if (!is.na(pmatch("D", unlist(strsplit(which, ""))))) 
            for (n in z[substring(z, 1, 1) == "D"]) delete.node(d, 
                n)
        if (!is.na(pmatch("Q", unlist(strsplit(which, ""))))) 
            for (n in z[substring(z, 1, 1) == "Q"]) delete.node(d, 
                n)
    }
}
extract <-
function (db, epg, gts, C, markers = sort(unique(db$marker))) 
{
    z <- list()
    for (m in markers) if (m != "AMEL") {
        e <- subset(epg, marker == m)
        z[[m]]$R <- sort(e$allele[e$height > C])
        for (gt in gts) {
            q <- get(gt)
            q <- sort(unlist(q[q$marker == m, 2:3]))
            names(q) <- NULL
            z[[m]][[gt]] <- q
        }
        dbm <- subset(db, marker == m)
        z[[m]]$alleles <- dbm$allele
        z[[m]]$frequency <- dbm$freq
    }
    z
}
fuzz <-
function (mixture, eps = 1e-08) 
{
    for (m in mixture$markers) {
        d <- mixture$domains[[m]]
        z <- get.nodes(d)
        for (n in z[substring(z, 1, 1) == "n"]) {
            tab <- get.table(d, n)
            tab$Freq <- pmax(eps, tab$Freq)
            set.table(d, n, tab, type = "cpt")
        }
    }
}
gt2aca <-
function (mixture, gt, eps = 0) 
{
    logGt <<- 0
    aca <- list()
    namesaca <- NULL
    for (m in 1:length(mixture$markers)) {
        md <- match(mixture$markers[m], as.character(gt$marker))
        if (!is.na(md)) {
            a1 <- match(gt$allele1[md], mixture$data[[m]]$allele)
            a2 <- match(gt$allele2[md], mixture$data[[m]]$allele)
            if (is.na(a1) || is.na(a2)) 
                warning("invalid allele in genotype")
            else {
                ac <- rep(0, length(mixture$data[[m]]$allele))
                ac[a1] <- ac[a1] + 1
                ac[a2] <- ac[a2] + 1
                if (eps > 0) {
                  AC <- matrix(eps, length(ac), 3)
                  AC[cbind(1:length(ac), ac + 1)] <- 1
                  aca <- c(aca, list(AC))
                }
                else {
                  aca <- c(aca, list(ac))
                }
                namesaca <- c(namesaca, mixture$markers[m])
                fq <- prod(mixture$data[[m]]$freq[c(a1, a2)])
                if (a1 != a2) 
                  fq <- 2 * fq
                logGt <<- logGt + log(fq)
            }
        }
    }
    names(aca) <- namesaca
    attr(aca, "logGt") <- logGt
    aca
}
intoMix <-
function (res00) 
{
    res0 <- subset(res00, !is.na(DNA))
    res0$Sim <- match(res0$Sim, unique(res0$Sim))
    res0
}
logL.UKX <-
function (mixture, initialize = TRUE) 
{
    if (initialize) 
        lapply(mixture$domains, initialize.domain)
    C <- mixture$C
    n.unknown <- mixture$n.unknown
    U <- mixture$U
    K <- mixture$K
    n_K <- lapply(mixture$data, function(d) subset(d, select = K))
    function(pararray) {
        logL.m <- function(m) {
            domain <- mixture$domains[[m]]
            d <- mixture$data[[m]]
            rho <- pararray[, "rho"]
            xi <- pararray[, "xi"]
            eta <- pararray[, "eta"]
            phi <- pararray[, "phi"]
            for (r in seq_len(mixture$ntraces)) {
                if (!all(is.na(d[, r + 1]))) {
                  evidence <- setCPT.O(domain, rho[[r]], xi[[r]], 
                    eta[[r]], phi[[r]][c(U, K)], d[, r + 1], 
                    C[[r]], n.unknown, n_K[[m]], attr(domain, 
                      "O")[[r]], d$gets_stutter, d$can_stutter, 
                    d$stutter.from)
                  lapply(seq_along(attr(domain, "O")[[r]]), function(i) set.finding(domain, 
                    attr(domain, "O")[[r]][i], evidence[, i]))
                }
            }
            extra.findings(domain, names(mixture$domains)[m])
            propagate(domain)
            get.normalization.constant(domain, log = TRUE)
        }
        sum(sapply(seq_along(mixture$domains), logL.m))
    }
}
logLX <-
function (mixture, presence.only = FALSE, initialize = TRUE) 
{
    if (mixture$n.unknown > 0) {
        if (presence.only) 
            logLpres.UK(mixture, initialize = initialize)
        else logL.UKX(mixture, initialize = initialize)
    }
    else {
        if (presence.only) 
            logLpres.K(mixture)
        else logL.K(mixture)
    }
}
loop.rpt.IBD <-
function (listdata, pars, IBD, typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs, quiet = FALSE, verbose = FALSE, 
    presence.only = FALSE, ...) 
{
    if (!is.null(targets)) {
        inds <- match(targets, contribs)
        inds <- inds[!is.na(inds)]
        if (!quiet) 
            cat("inds", inds, "\n")
        jcontr <- match(contribs, targets)
        jcontr <- jcontr[!is.na(jcontr)]
        if (!quiet) 
            cat("jcontr", jcontr, "\n")
        jtyped <- match(names(typed.gts), targets)
        if (!quiet) 
            cat("jtyped", jtyped, "\n")
    }
    IBD <- convertIBD(IBD)
    sumlogLR <- 0
    logLR<-NULL
    for (t in 1:length(listdata)) {
        data <- listdata[[t]]
        for (m in setdiff(unique(data$marker), "AMEL")) if (all(unlist(lapply(typed.gts, 
            function(x) {
                m %in% x$marker
            })))) {
            listdatam <- lapply(listdata, function(d) subset(d, 
                d$marker == m))
            mixDm <- DNAmixture(listdatam, ...)
            baseline <- protected(logL(mixDm, presence.only)(pars))
            LR <- 0
            wtsum <- 0
            for (r in 1:nrow(IBD$patt)) {
                mixDmr <- DNAmixture(listdatam, triangulate = F, 
                  compile = F, ...)
                wt <- IBD$pr[r] * rpt.IBD(mixDmr, list(patt = IBD$patt[r, 
                  ]), typed.gts, inds, jtyped, jcontr)
                if (length(wt) > 0 && wt > 0) {
                  log10LR <- (protected(logL(mixDmr, presence.only)(pars)) - 
                    baseline)/log(10)
                  LR <- LR + wt * 10^log10LR
                  wtsum <- wtsum + wt
                }
            }
            LR <- LR/wtsum
            if (verbose) 
                cat(m, "log10 LR", log10(LR), "; LR", 
                  LR, "\n")
		logLR<-c(logLR,log10(LR))
            sumlogLR <- sumlogLR + log10(LR)
        }
    }
    if (verbose) 
        cat("all trace all marker overall log10 LR", sumlogLR, 
            "LR", 10^sumlogLR, "\n")
    names(sumlogLR)<-'sumlog10LR'
    invisible(structure(sumlogLR,log10LR=logLR))
}
make.findings <-
function (z) 
{
    tf <- tempfile()
    cat("f<-function(domain,name){\n", file = tf)
    if (length(z) > 0) 
        for (j in 1:length(z)) {
            if (z[[j]][1] == "Male") {
                ind <- z[[j]]$ind
                cat("if(name==\"AMEL\"){\n", file = tf, append = T)
                cat("\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_1\",sep=\"\"), c(0,1,0))\n", sep = "", 
                  file = tf, append = T)
                cat("\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_2\",sep=\"\"), c(0,1,0))\n", sep = "", 
                  file = tf, append = T)
                cat("}\n", file = tf, append = T)
            }
            else if (z[[j]][1] == "Female") {
                ind <- z[[j]]$ind
                cat("if(name==\"AMEL\"){\n", file = tf, append = T)
                cat("\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_1\",sep=\"\"), c(0,0,2))\n", sep = "", 
                  file = tf, append = T)
                cat("\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_2\",sep=\"\"), c(2,0,0))\n", sep = "", 
                  file = tf, append = T)
                cat("}\n", file = tf, append = T)
            }
            else if (z[[j]][1] == "Aca") {
                ind <- z[[j]]$ind
                aca <- z[[j]]$aca
                cat("w<-", aca, "[[name]]\n", "if(!is.null(w))\n", 
                  "\t{\n", "\tif(is.null(dim(w))) for(k in 1:length(w))\n", 
                  "\tset.finding(domain,paste(\"n_\",", ind, 
                  ",\"_\",k,sep=\"\"),w[k])\n", "\telse for(k in 1:nrow(w))\n", 
                  "\tset.finding(domain,paste(\"n_\",", ind, 
                  ",\"_\",k,sep=\"\"),w[k,])\n", "\t}\n", sep = "", 
                  file = tf, append = T)
            }
            else if (z[[j]][1] == "Caca") {
                aca <- z[[j]]$aca
                cat("if(!is.null(", aca, "[[name]]))", " for(k in 1:length(", 
                  aca, "[[name]]))\n", "\tset.finding(domain,paste(\"Cn\",k,sep=\"\"),", 
                  aca, "[[name]][k])\n", sep = "", file = tf, 
                  append = T)
            }
            else if (z[[j]][1] == "Denom") {
                cat("\tset.finding(domain, \"denom\", 1)\n", 
                  sep = "", file = tf, append = T)
            }
            else if (z[[j]][1] == "Rlikd") {
                aca <- z[[j]]$aca
                cgt <- z[[j]]$cgt
                evid <- z[[j]]$evid
                cat("if(!is.null(", aca, "[[name]]))\n", "\tset.finding(domain,\"Rlikd\",c(0,", 
                  evid, "[[name]]))\n", sep = "", file = tf, 
                  append = T)
            }
        }
    cat("}\n", file = tf, append = T)
    source(tf, local = TRUE)
    f
}
make.profile <-
function (gt, name = "K") 
{
    marker <- NULL
    allele <- NULL
    count <- NULL
    for (i in 1:nrow(gt)) {
        a1 <- gt$allele1[i]
        a2 <- gt$allele2[i]
        if (a1 == a2) {
            marker <- c(marker, gt$marker[i])
            allele <- c(allele, a1)
            count <- c(count, 2)
        }
        else {
            marker <- c(marker, gt$marker[i], gt$marker[i])
            allele <- c(allele, a1, a2)
            count <- c(count, 1, 1)
        }
    }
    out <- data.frame(marker = marker, allele = allele, count = count)
    names(out)[3] <- name
    out
}
mixMLX <-
function (mixture, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, ...) 
{
    R <- mixture$ntraces
    k <- mixture$k
    U <- mixture$U
    K <- mixture$K
    contr <- c(U, K)
    n.unknown <- mixture$n.unknown
    x2phi <- function(x) {
        phi <- if (phi.eq) {
            rep(list(tail(x, -3 * R)), R)
        }
        else {
            split(tail(x, -3 * R), rep(1:R, each = k))
        }
        mapply(function(th, contributor) {
            names(th) <- contributor
            th
        }, phi, rep(list(contr), times = R), SIMPLIFY = FALSE)
    }
    x2phiU <- function(x) {
        lapply(x2phi(x), function(th) {
            head(th, n.unknown)
        })
    }
    x2arr <- function(x) {
        arr <- array(list(NULL), dim = c(R, 4), dimnames = list(NULL, 
            c("rho", "eta", "xi", "phi")))
        arr[1:(3 * R)] <- as.list(head(x, 3 * R))
        arr[, "phi"] <- x2phi(x)
        arr
    }
    parlist2x <- function(parlist) {
        rex <- unlist(parlist[, 1:3], use.names = FALSE)
        if (phi.eq) {
            phi <- parlist[[1, 4]][contr]
        }
        else {
            phi <- unlist(lapply(parlist[, 4], function(x) x[contr]), 
                use.names = FALSE)
        }
        c(rex, phi)
    }
    logl <- logLX(mixture)
    funvals <- numeric(0)
    minus.loglikelihood <- function(x) {
        xs <- x2arr(x)
        if (trace) 
            print.mixpar(xs)
        val <- -logl(xs)
        if (trace) 
            print(-val)
        val
    }
    lb <- rep(0, times = 3 * R + ifelse(phi.eq, k, k * R))
    ub <- rep(c(Inf, 1), times = c(2 * R, R + ifelse(phi.eq, 
        k, k * R)))
    if (phi.eq) {
        phi.sum.constraint <- function(x) {
            sum(tail(x, -3 * R))
        }
        eqB <- 1
    }
    else {
        phi.sum.constraint <- function(x) {
            sapply(x2phi(x), sum)
        }
        eqB <- rep(1, R)
    }
    if (!missing(constraints)) {
        eqfun <- function(x) {
            c(phi.sum.constraint(x), do.call(constraints, list(x2arr(x))))
        }
        eqB <- c(eqB, val)
    }
    else {
        eqfun <- phi.sum.constraint
    }
    if (n.unknown > 1 & order.unknowns) {
        phi.symmetry <- function(x) {
            phiU <- x2phiU(x)[[1]]
            diff(phiU, lag = 1)
        }
        n.diffs <- (n.unknown - 1)
        ineqLB <- rep(-1, n.diffs)
        ineqUB <- rep(0, n.diffs)
    }
    else phi.symmetry <- ineqLB <- ineqUB <- NULL
    x0 <- parlist2x(pars)
    soln <- solnp(x0, fun = minus.loglikelihood, LB = lb, UB = ub, 
        eqfun = eqfun, eqB = eqB, ineqfun = phi.symmetry, ineqLB = ineqLB, 
        ineqUB = ineqUB, control = list(trace = 0), ...)
    est <- x2arr(soln$pars)
    class(est) <- "mixpar"
    val <- -tail(soln$value, 1)
    out <- list(mle = est, lik = val, funvals = funvals, starting.point = x0, 
        minimization.output = soln)
}
ped2fam <-
function (H) 
{
    id <- H$ID
    fidx <- H$FIDX
    fidx[fidx == 0] <- NA
    midx <- H$MIDX
    midx[midx == 0] <- NA
    sex <- H$SEX
    if (any(sex != 1 & sex != 2)) {
        warning("all sexes should be male or female: others have been coerced to male")
        sex[sex != 2] <- 1
    }
    Familias::FamiliasPedigree(id = id, dadid = id[fidx], momid = id[midx], 
        sex = c("male", "female")[sex])
}
pedigreeIBD <-
function (x, target = c("c1", "c2", "c3"), 
    cond = TRUE, quiet = TRUE, verbose = FALSE) 
{
    if (is.ped(x)) {
        pdf <- as.data.frame(x)
        id <- pdf$id
        trim <- (id %in% target) | (id %in% ancestors(x, target))
        id <- id[trim]
        fid <- pdf$fid[trim]
        mid <- pdf$mid[trim]
        sex <- pdf$sex[trim]
    }
    else {
        id <- fid <- mid <- sex <- NULL
        for (ic in 1:length(x)) {
            pdf <- as.data.frame(x[[ic]])
            idi <- pdf$id
            fidi <- pdf$fid
            midi <- pdf$mid
            sexi <- pdf$sex
            trim <- (idi %in% target) | (idi %in% ancestors(x[[ic]], 
                target[target %in% idi]))
            id <- c(id, idi[trim])
            fid <- c(fid, fidi[trim])
            mid <- c(mid, midi[trim])
            sex <- c(sex, sexi[trim])
        }
    }
    x <- list(ID = id, FIDX = match(fid, id, nomatch = 0), MIDX = match(mid, 
        id, nomatch = 0), SEX = sex)
    ida <- id[1]
    pattern <- matrix(1:2, 1, 2)
    mlp <- 2
    count <- 1
    for (i in 2:length(id)) {
        idnext <- id[i]
        iid <- match(idnext, id)
        ida <- c(ida, idnext)
        wf <- match(fid[iid], ida)
        if (is.na(wf)) {
            pattern <- cbind(pattern, max(pattern) + 1)
        }
        else {
            pattern <- rbind(cbind(pattern, pattern[, 2 * wf - 
                1]), cbind(pattern, pattern[, 2 * wf]))
            count <- c(count, count)
        }
        wm <- match(mid[iid], ida)
        if (is.na(wm)) {
            pattern <- cbind(pattern, max(pattern) + 1)
        }
        else {
            pattern <- rbind(cbind(pattern, pattern[, 2 * wm - 
                1]), cbind(pattern, pattern[, 2 * wm]))
            count <- c(count, count)
        }
        keep <- prune(x, target, ida, verbose) | (ida %in% target)
        ida <- ida[keep]
        pattern <- pattern[, rep(keep, rep(2, length(keep))), 
            drop = FALSE]
        if (cond) 
            pattern <- t(apply(pattern, 1, function(r) minimalPattern(r)))
        else pattern <- t(apply(pattern, 1, function(r) match(r, 
            unique.default(r))))
        o <- statnet.common::order(pattern)
        pattern <- pattern[o, , drop = FALSE]
        nrp <- nrow(pattern)
        w <- rep(0, nrp)
        if (nrp > 1) 
            for (j in 1:(nrp - 1)) w[j] <- as.integer(any(pattern[j, 
                ] != pattern[j + 1, ]))
        w[nrp] <- 1
        cumct <- cumsum(count[o])
        count <- diff(c(0, cumct[which(w == 1)]))
        pattern <- pattern[w == 1, , drop = FALSE]
        mlp <- max(mlp, length(pattern))
        if (verbose) {
            cat(idnext, dim(pattern), "\n")
            dimnames(pattern) <- list(rep(" ", nrow(pattern)), 
                as.vector(rbind(ida, " ")))
            print(pattern)
        }
    }
    dimnames(pattern) <- list(rep(" ", nrow(pattern)), 
        as.vector(rbind(ida, " ")))
    mti <- match(target, ida)
    pattern <- pattern[, as.vector(rbind(2 * mti - 1, 2 * mti))]
    if (!quiet) 
        print(cbind(pr = count/sum(count), pattern = pattern))
    IBD <- structure(list(pr = count/sum(count), patt = pattern), 
        ped = ped(id, fid, mid, sex), mlp = mlp, target=target, class='IBD')
    IBD
}
plot.IBD <-
function (pattern, labels = NULL, probs = NULL, order = NULL, 
    colrs = c("black", "red", "blue"), digits = 3, nr = ceiling(sqrt(np))) 
{
    if (is.list(pattern)) {
        probs <- pattern$pr
        pattern <- pattern$patt
    }
    np <- nrow(pattern)
    if (!is.null(labels) && is.na(labels)) 
        labels <- apply(pattern, 1, function(x) paste(x, collapse = " "))
    chcol <- -1
    if (!is.null(order)) {
        if (length(order) > 1) {
            o <- order(order)
            chcol <- which(labels[o][-1] != labels[o][-np])
        }
        else {
            switch(pmatch(order, c("pattern", "probs", "labels")), 
                {
                  nout <- apply(pattern, 1, function(p) sum(0 != 
                    apply(matrix(p, nrow = 2), 2, diff)))
                  ndist <- apply(pattern, 1, max)
                  maxf <- apply(pattern, 1, function(p) max(tabulate(p, 
                    6)))
                  z <- apply(pattern, 1, function(p) sum(apply(matrix(apply((outer(p, 
                    p, "-") + diag(6)) == 0, 1, sum), 2, 3), 
                    2, min)))
                  o <- statnet.common::order(cbind(nout, ndist, 
                    maxf, z))
                  zz <- cbind(nout, ndist, maxf, z)[o, ]
                  chcol <- which(apply(apply(zz, 2, diff) != 
                    0, 1, any))
                }, {
                  o <- order(probs, decreasing = TRUE)
                  chcol <- which(diff(probs[o]) != 0)
                }, {
                  if (is.null(labels)) o <- 1:np else {
                    o <- order(labels)
                    chcol <- which(labels[o][-1] != labels[o][-np])
                  }
                })
            pattern <- pattern[o, ]
            if (!is.null(labels)) 
                labels <- labels[o]
            if (!is.null(probs)) 
                probs <- probs[o]
        }
    }
    nc <- ceiling(np/nr)
    opars <- par(mfrow = c(nr, nc), mar = c(0.5, 0.5, 0.5, 0.5))
    nind <- ncol(pattern)/2
    y <- rep((nind - 1):0, each = 2)
    x <- rep(c(0, 1), nind)
    curve <- spline(seq(0, 1, len = 5), c(0, 2/3, 1, 2/3, 0))
    sc <- c(0, 0.15, 0.4)
    icol <- length(colrs)
    for (f in 1:np) {
        patt <- pattern[f, ]
        plot(x, y, xlim = c(-0.5, 1.5), ylim = c(-0.9, nind - 
            0.5), type = "n", xaxt = "n", yaxt = "n", bty = "n", 
            asp = 1)
        xo <- yo <- 0
        for (i in 1:(2 * nind - 1)) for (j in (i + 1):(2 * nind)) if (patt[i] == 
            patt[j]) {
            if (i%%2 == 1 && j%%2 == 1 & abs(i - j) > 2) 
                lines(xo + x[j] - sc[abs(i - j)/2] * curve$y, 
                  yo + y[j] + (y[i] - y[j]) * curve$x, col = colrs[icol])
            else if (i%%2 == 0 && j%%2 == 0 & abs(i - j) > 2) 
                lines(xo + x[j] + sc[abs(i - j)/2] * curve$y, 
                  yo + y[j] + (y[i] - y[j]) * curve$x, col = colrs[icol])
            else lines(xo + x[c(i, j)], yo + y[c(i, j)], col = colrs[icol])
        }
        points(xo + x, yo + y, pch = 16)
        yt <- 0.2
        if (!is.null(labels)) {
            yt <- yt - 0.5
            text(0.5, yt, labels[f], cex = 1)
        }
        if (!is.null(probs)) {
            yt <- yt - 0.5
            text(0.5, yt, signif(probs[f], digits), cex = 1)
        }
        if (f %in% chcol) 
            icol <- (icol%%length(colrs)) + 1
    }
    par(opars)
}
print.IBD <-
function(IBD,targets=letters[1:(nc/2)])
{
IBD<-convertIBD(IBD)
nr<-length(IBD$pr)
nc<-length(IBD$patt)/nr
res<-cbind(IBD$pr,matrix(IBD$patt,nr,nc))
if(!is.null(attr(IBD,'target'))) targets<-attr(IBD,'target')
prmatrix(cbind(IBD$pr,IBD$patt),rowlab=rep('',nr),collab=c('pr',as.vector(rbind(targets,rep('',nc/2)))))
}
print.tablesize <-
function (s) 
{
    cat("total table size", s, "\n")
    invisible(s)
}
process.patterns <-
function (kinship, igt, jtyped, jcontr, q, keeplabels) 
{
    ntyped <- length(jtyped)
    nperms <- 2^ntyped
    apa <- NULL
    ape <- NULL
    pdenom <- NULL
    type <- NULL
    val <- NULL
    npp <- 0
    rlen <- 2 * length(jcontr)
    for (ipatt in 1:length(kinship$pr)) {
        kt <- kinship$patt[ipatt, as.vector(outer(c(-1, 0), 2 * 
            jtyped, "+"))]
        ukt <- unique(kt)
        ap <- rep(0, max(kinship$patt))
        for (iperm in 1:nperms) {
            ipp <- nperms * (ipatt - 1) + iperm
            possible <- TRUE
            if (ntyped != 0) {
                w <- binary(iperm - 1, ntyped)
                pkt <- kt[as.vector(rbind(w, 1 - w)) + rep(2 * 
                  (1:ntyped) - 1, rep(2, ntyped))]
                ap[pkt] <- igt
                possible <- all(ap[pkt] == igt)
            }
            if (possible) {
                apa <- c(apa, ipatt)
                ape <- c(ape, iperm)
                type <- rbind(type, rep(" ", rlen))
                val <- rbind(val, rep(0, rlen))
                if (is.null(dim(q))) 
                  pqa <- prod(q[ap[ukt]])
                else pqa <- prod(q[cbind(ukt, ap[ukt])])
                pdenom <- c(pdenom, kinship$pr[ipatt] * pqa)
                npp <- npp + 1
                for (ij in 1:(2 * length(jcontr))) {
                  j <- as.vector(outer(c(-1, 0), 2 * jcontr, 
                    "+"))[ij]
                  ku <- kinship$patt[ipatt, j]
                  for (ik in 1:length(ku)) {
                    ia <- match(ku[ik], ukt)
                    if (is.na(ia)) {
                      type[npp, ij] <- "d"
                      val[npp, ij] <- ku[ik]
                    }
                    else {
                      type[npp, ij] <- "c"
                      val[npp, ij] <- ap[ukt[ia]]
                    }
                  }
                }
            }
        }
    }
    res <- NULL
    if (!is.null(val)) {
        oval <- val
        if (keeplabels) {
            todraw <- NULL
            for (ir in 1:nrow(val)) {
                w <- which(type[ir, ] == "d")
                vw <- val[ir, w]
                todraw <- c(todraw, vw)
            }
            todraw <- unique(todraw)
        }
        else {
            ntodraw <- 0
            for (ir in 1:nrow(val)) {
                w <- which(type[ir, ] == "d")
                vw <- val[ir, w]
                uvw <- unique(vw)
                oval[ir, w] <- match(vw, uvw)
                ntodraw <- max(ntodraw, match(vw, uvw))
            }
            if (ntodraw > 0) 
                todraw <- 1:ntodraw
            else todraw <- NULL
        }
        res <- data.frame(apa, ape, pdenom, type = type, oval = oval, 
            stringsAsFactors = FALSE)
        attr(res, "todraw") <- todraw
        attr(res, "npp") <- npp
    }
    res
}
protected <-
function (x, default = -Inf) 
{
    if (is.numeric(tryCatch.W.E(x)$val)) 
        x
    else default
}
protected.mixML <-
function (mixture, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, default = -999999, 
    ...) 
{
    R <- mixture$ntraces
    k <- mixture$k
    U <- mixture$U
    K <- mixture$K
    contr <- c(U, K)
    n.unknown <- mixture$n.unknown
    x2phi <- function(x) {
        phi <- if (phi.eq) {
            rep(list(tail(x, -3 * R)), R)
        }
        else {
            split(tail(x, -3 * R), rep(1:R, each = k))
        }
        mapply(function(th, contributor) {
            names(th) <- contributor
            th
        }, phi, rep(list(contr), times = R), SIMPLIFY = FALSE)
    }
    x2phiU <- function(x) {
        lapply(x2phi(x), function(th) {
            head(th, n.unknown)
        })
    }
    x2arr <- function(x) {
        arr <- array(list(NULL), dim = c(R, 4), dimnames = list(NULL, 
            c("rho", "eta", "xi", "phi")))
        arr[1:(3 * R)] <- as.list(head(x, 3 * R))
        arr[, "phi"] <- x2phi(x)
        arr
    }
    parlist2x <- function(parlist) {
        rex <- unlist(parlist[, 1:3], use.names = FALSE)
        if (phi.eq) {
            phi <- parlist[[1, 4]][contr]
        }
        else {
            phi <- unlist(lapply(parlist[, 4], function(x) x[contr]), 
                use.names = FALSE)
        }
        c(rex, phi)
    }
    logl <- logL(mixture)
    funvals <- numeric(0)
    minus.loglikelihood <- function(x) {
        xs <- x2arr(x)
        if (trace) 
            print.mixpar(xs)
        val <- -protected(logl(xs), default)
        if (trace) 
            print(-val)
        val
    }
    lb <- rep(0, times = 3 * R + ifelse(phi.eq, k, k * R))
    ub <- rep(c(Inf, 1), times = c(2 * R, R + ifelse(phi.eq, 
        k, k * R)))
    if (phi.eq) {
        phi.sum.constraint <- function(x) {
            sum(tail(x, -3 * R))
        }
        eqB <- 1
    }
    else {
        phi.sum.constraint <- function(x) {
            sapply(x2phi(x), sum)
        }
        eqB <- rep(1, R)
    }
    if (!missing(constraints)) {
        eqfun <- function(x) {
            c(phi.sum.constraint(x), do.call(constraints, list(x2arr(x))))
        }
        eqB <- c(eqB, val)
    }
    else {
        eqfun <- phi.sum.constraint
    }
    if (n.unknown > 1 & order.unknowns) {
        phi.symmetry <- function(x) {
            phiU <- x2phiU(x)[[1]]
            diff(phiU, lag = 1)
        }
        n.diffs <- (n.unknown - 1)
        ineqLB <- rep(-1, n.diffs)
        ineqUB <- rep(0, n.diffs)
    }
    else phi.symmetry <- ineqLB <- ineqUB <- NULL
    x0 <- parlist2x(pars)
    soln <- solnp(x0, fun = minus.loglikelihood, LB = lb, UB = ub, 
        eqfun = eqfun, eqB = eqB, ineqfun = phi.symmetry, ineqLB = ineqLB, 
        ineqUB = ineqUB, control = list(trace = 0, delta = 1e-07, 
            tol = 1e-08), ...)
    est <- x2arr(soln$pars)
    class(est) <- "mixpar"
    val <- -tail(soln$value, 1)
    out <- list(mle = est, lik = val, funvals = funvals, starting.point = x0, 
        minimization.output = soln)
}
prune <-
function (x, target, S, verbose = FALSE) 
{
    id <- x$ID
    ifid <- x$FIDX
    imid <- x$MIDX
    T <- NULL
    stack <- target
    repeat {
        now <- stack[1]
        stack <- stack[-1]
        if (!now %in% S) {
            inow <- match(now, id)
            f <- ifid[inow]
            if (f != 0) {
                cf <- id[f]
                if (cf %in% S) 
                  T <- unique(c(T, cf))
                else if (!cf %in% stack) 
                  stack <- c(stack, cf)
            }
            m <- imid[inow]
            if (m != 0) {
                cm <- id[m]
                if (cm %in% S) 
                  T <- unique(c(T, cm))
                else if (!cm %in% stack) 
                  stack <- c(stack, cm)
            }
        }
        if (length(stack) == 0) 
            break
    }
    if (verbose) 
        cat("can delete", S[!S %in% T], "\nmust keep ", S[S %in% 
            T], "\n")
    S %in% T
}
read.db <-
function (file = "db.csv", sep = ifelse(Sys.getenv("USERNAME") == 
    "Julia", ";", ",")) 
{
    z <- read.csv(file, stringsAsFactors = FALSE, sep = sep)
    data.frame(marker = as.factor(z$marker), allele = as.numeric(z$allele), 
        height = z$frequency)
}
read.epg <-
function (file = "epg.csv", sep = ifelse(Sys.getenv("USERNAME") == 
    "Julia", ";", ",")) 
{
    z <- read.csv(file, stringsAsFactors = FALSE, sep = sep)
    data.frame(marker = as.factor(z$marker), allele = as.numeric(z$allele), 
        height = z$height)
}
replace.tables.for.UAF <-
function (mixture, M, compile = TRUE) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        na <- nrow(mixture$data[[m]])
        q <- mixture$data[[m]]$freq
        q <- q/sum(q)
        alpha <- M * q
        beta <- c(rev(cumsum(rev(alpha)))[-1], 0)
        z <- get.nodes(d)
        n.unknown <- length(grep("^n_.*_1$", z))
        for (n in z[substring(z, 1, 1) == "n"]) for (par in get.parents(d, 
            n)) delete.edge(d, n, par)
        n <- outer(1:n.unknown, 1:na, function(i, a) paste("n", 
            i, a, sep = "_"))
        S <- outer(1:n.unknown, 1:na, function(i, a) paste("S", 
            i, a, sep = "_"))
        T <- outer(1:n.unknown, 1:na, function(i, a) paste("T", 
            i, a, sep = "_"))
        U <- outer(1:n.unknown, 1:na, function(i, a) paste("U", 
            i, a, sep = "_"))
        for (a in 1:na) {
            if (n.unknown > 1) 
                for (i in 1:(n.unknown - 1)) {
                  add.node(d, T[i, a], states = 0:(2 * i), subtype = "numbered")
                  add.node(d, U[i, a], states = 0:(2 * i), subtype = "numbered")
                }
            for (i in 1:n.unknown) {
                if (a == 1) {
                  if (i == 1) {
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), rep(alpha[1], lnval), rep(beta[1], 
                      lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), alpha[a] + tab[, T[i - 1, a]], 
                      (beta[a] + pmax(0, tab[, U[i - 1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                else {
                  if (i == 1) {
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], rep(alpha[a], lnval), 
                      rep(beta[a], lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], alpha[a] + tab[, T[i - 
                      1, a]], (beta[a] + pmax(0, tab[, U[i - 
                      1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, T[i, a], n[i, a])
                    tab <- get.table(d, T[i, a])
                    tab$Freq <- ifelse(tab[, T[i, a]] == tab[, 
                      n[i, a]], 1, 0)
                    set.table(d, T[i, a], tab)
                  }
                  else {
                    add.edge(d, T[i, a], n[i, a])
                    add.edge(d, T[i, a], T[i - 1, a])
                    tab <- get.table(d, T[i, a])
                    sumvals <- tab[, n[i, a]] + tab[, T[i - 1, 
                      a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      T[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, T[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, U[i, a], S[i, a])
                    tab <- get.table(d, U[i, a])
                    tab$Freq <- ifelse(tab[, U[i, a]] == 2 - 
                      tab[, S[i, a]], 1, 0)
                    set.table(d, U[i, a], tab)
                  }
                  else {
                    add.edge(d, U[i, a], S[i, a])
                    add.edge(d, U[i, a], U[i - 1, a])
                    tab <- get.table(d, U[i, a])
                    sumvals <- 2 - tab[, S[i, a]] + tab[, U[i - 
                      1, a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      U[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, U[i, a], tab)
                  }
                }
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
replace.Ui.tables <-
function (mixture, aca, ind = 1) 
{
    if (is.data.frame(aca)) 
        aca <- gt2aca(mixture, aca)
    for (m in names(mixture$domains)) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            nC <- aca[[m]]
            na <- length(nC)
            if (max(aca[[m]]) == 2) {
                am <- which(aca[[m]] == max(aca[[m]]))
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a] - (am == 1), 
                      1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a] - (am == a), 
                      pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                        0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, tab[, S2am1]), 1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                am <- ab[1]
                bm <- ab[2]
                add.node(d, "ibdyet", subtype = "boolean")
                tab <- get.table(d, "ibdyet")
                tab$Freq <- c(0.5, 0.5)
                set.table(d, "ibdyet", tab, type = "cpt")
                add.edge(d, paste("n_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("n_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == am) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                          0, s, 1))
                    }
                  }
                  else if (a == bm) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], 1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], pmax(0, 1 - tab[, S2am1]), 
                        q[a]/ifelse(s > 0, s, 1))
                    }
                  }
                  else if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a], 1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a], pmax(0, 1 - 
                      tab[, S2am1]), q[a]/ifelse(s > 0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                  if (a == bm) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
        }
        compile(d)
    }
}
require.compiled <-
function (mixture) 
{
    for (m in mixture$markers) {
        d <- mixture$dom[[m]]
        if (!unclass(summary(d))$domain$compiled) {
            compile(d)
            cat(m, "compiled\n")
        }
        else cat(m, "already compiled\n")
    }
}
require.edge <-
function (d, n, p) 
{
    if (!p %in% get.parents(d, n)) 
        add.edge(d, n, p)
}
rGTs <-
function (nreps, IBD, db, DNA, sex = rep(0, ncontr), nU = 0) 
{
    IBD <- convertIBD(IBD)
    if (nU > 0) 
        IBD <- addU(IBD, nU)
    ncontr <- ncol(IBD$patt)/2
    res <- list()
    markers <- as.character(unique(db$marker))
    nmark <- length(markers)
    z <- matrix(0, nmark, ncol(IBD$patt))
    dimnames(z) <- list(markers, NULL)
    for (ir in 1:nreps) {
        for (m in markers) {
            dbm <- db[db$marker == m, ]
            zi <- IBD$patt[sample(length(IBD$pr), 1, repl = TRUE, 
                pr = IBD$pr), ]
            g <- sample(dbm$allele, max(IBD$patt), repl = TRUE, 
                pr = dbm$frequency)
            z[m, ] <- g[zi]
        }
        Sim <- rep(1:ncontr, each = 2 * nmark)
        Marker <- rep(markers, each = 2, times = ncontr)
        Allele <- as.character(as.vector(aperm(array(z, c(nmark, 
            2, ncontr)), c(2, 1, 3))))
        w <- Marker == "AMEL"
        if (!is.null(w)) {
            amel <- rep("X", 2 * ncontr)
            for (ic in 1:ncontr) amel[2 * ic] <- switch(sex[ic] + 
                1, sample(c("X", "Y"), 1), "Y", "X")
            Allele[w] <- amel
        }
        sn <- rep("rGTs", 2 * ncontr * nmark)
        res[[ir]] <- data.frame(Sim, Sample.Name = sn, Marker, 
            Allele, DNA = DNA[Sim], stringsAsFactors = FALSE)
    }
    res
}
rni <-
function (seed = 0) 
{
    if (seed == 0) {
        op <- options(digits.secs = 6)
        tt <- Sys.time()
        options(op)
        st <- strftime(tt, "%s")
        t <- as.numeric(paste(substring(st, nchar(st) - 5), substring(as.character(tt), 
            21, 23), sep = ""))
        set.seed(t)
        seed <- sample(999999999, 1)
        cat("to rerun type rni(", seed, ")\n", sep = "")
    }
    set.seed(seed)
    invisible(seed)
}
rpt.AMEL <-
function (mixture, sex, compile = TRUE) 
{
    if ("AMEL" %in% mixture$markers) {
        d <- mixture$domains[["AMEL"]]
        for (i in seq_along(sex)) {
            ni1 <- paste("n", i, 1, sep = "_")
            if (ni1 %in% get.nodes(d)) {
                tab <- get.table(d, ni1)
                freq <- switch(sex[i] + 1, c(0, 0.5, 0.5), c(0, 
                  1, 0), c(0, 0, 1))
                tab$Freq <- freq
                set.table(d, ni1, tab, type = "cpt")
            }
        }
    }
    if (compile) 
        for (d in mixture$domains) compile(d)
}
rpt.IBD <-
function (mixture, IBD, typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs, quiet = FALSE, all.freq = NULL, 
    compile = TRUE) 
{
    if (!is.null(targets)) {
        jcontr <- match(contribs, targets)
	  inds<-which(!is.na(jcontr))
        if (!quiet) 
            cat("inds", inds, "\n")
        jcontr <- jcontr[!is.na(jcontr)]
        if (!quiet) 
            cat("jcontr", jcontr, "\n")
        jtyped <- match(names(typed.gts), targets)
        if (!quiet) 
            cat("jtyped", jtyped, "\n")
    }
    IBD <- convertIBD(IBD)
    if (length(jtyped) != length(typed.gts)) 
        stop("jtyped and typed.gts incompatible")
    if (length(jcontr) != length(inds)) 
        stop("jcontr and inds incompatible")
    if (!all(jtyped %in% (1:ncol(IBD$patt)/2)) || !all(jcontr %in% 
        (1:ncol(IBD$patt)/2))) 
        stop("jtyped, jcontr and IBD incompatible")
    ntyped <- length(jtyped)
    ap <- rep(0, max(IBD$patt))
    ptypedgts <- NULL
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (all(unlist(lapply(typed.gts, function(x) {
            m %in% x$marker
        })))) {
            keeplabels <- FALSE
            if (missing(all.freq) || is.null(all.freq)) {
                q <- mixture$data[[m]]$freq
                q <- q/sum(q)
            }
            else {
                if (is.data.frame(all.freq)) {
                  dbm <- subset(all.freq, marker == m)
                  q <- dbm$frequency[match(mixture$data[[m]]$allele, 
                    dbm$allele)]
                  q <- q/sum(q)
                }
                else {
                  keeplabels <- TRUE
                  q <- matrix(0, max(IBD$patt), nrow(mixture$data[[m]]))
                  for (k in 1:max(IBD$patt)) {
                    dbm <- subset(all.freq[[k]], marker == m)
                    q[k, ] <- dbm$frequency[match(mixture$data[[m]]$allele, 
                      dbm$allele)]
                    q[k, ] <- q[k, ]/sum(q[k, ])
                  }
                }
            }
            igt <- NULL
            if (ntyped > 0) 
                for (l in 1:ntyped) {
                  gt <- typed.gts[[l]]
                  a1 <- match(gt$allele1[gt$marker == m], mixture$data[[m]]$allele)
                  a2 <- match(gt$allele2[gt$marker == m], mixture$data[[m]]$allele)
                  igt <- c(igt, a1, a2)
                }
            z <- process.patterns(IBD, igt, jtyped, jcontr, q, 
                keeplabels)
            if (is.null(z)) {
                ptypedgts <- c(ptypedgts, 0)
                names(ptypedgts)[length(ptypedgts)] <- m
            }
            else {
                nodes <- get.nodes(d)
                for (i in inds) {
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("n", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                  }
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("S", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                    delete.node(d, n)
                  }
                }
                add.node(d, "pattperm", states = 1:attr(z, 
                  "npp"))
                tab <- get.table(d, "pattperm")
                f <- as.vector(z$pdenom)
                ptypedgts <- c(ptypedgts, sum(f))
                names(ptypedgts)[length(ptypedgts)] <- m
                tab$Freq <- f/sum(f)
                set.table(d, "pattperm", tab)
                na <- nrow(mixture$data[[m]])
                for (k in attr(z, "todraw")) {
                  if (is.null(dim(q))) {
                    qstar <- q/rev(cumsum(rev(q)))
                  }
                  else {
                    qstar <- q[k, ]/rev(cumsum(rev(q[k, ])))
                  }
                  for (a in 1:na) {
                    npka <- paste("np", k, "_", a, 
                      sep = "")
                    add.node(d, npka, states = 0:1)
                    if (a == 1) {
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[1], qstar[1])
                    }
                    else {
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.edge(d, npka, Spka1)
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[a], qstar[a], 1, 
                        0)
                    }
                    set.table(d, npka, tab)
                    if (a != 1 & a != na) {
                      Spka <- paste("Sp", k, "_", 
                        a, sep = "")
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.node(d, Spka, states = 0:1)
                      add.edge(d, Spka, npka)
                      add.edge(d, Spka, Spka1)
                      tab <- get.table(d, Spka)
                      tab$Freq <- as.numeric(tab[, Spka] == pmax(tab[, 
                        npka], tab[, Spka1]))
                      set.table(d, Spka, tab)
                    }
                  }
                }
                for (l in 1:length(inds)) {
                  i <- inds[l]
                  j <- l
                  t <- cbind(z[[paste("type.", 2 * j - 
                    1, sep = "")]], z[[paste("type.", 
                    2 * j, sep = "")]])
                  npj <- cbind(paste("np", z[[paste("oval.", 
                    2 * j - 1, sep = "")]], sep = ""), 
                    paste("np", z[[paste("oval.", 
                      2 * j, sep = "")]], sep = ""))
                  ov <- cbind(z[[paste("oval.", 2 * j - 
                    1, sep = "")]], z[[paste("oval.", 
                    2 * j, sep = "")]])
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    add.edge(d, nia, "pattperm")
                    for (ipp in 1:attr(z, "npp")) {
                      if (t[ipp, 1] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 1], 
                          a, sep = "_"))
                      if (t[ipp, 2] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 2], 
                          a, sep = "_"))
                    }
                  }
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    tab <- get.table(d, nia)
                    for (ipp in 1:attr(z, "npp")) {
                      w <- tab[, "pattperm"] == ipp
                      value <- rep(0, sum(w))
                      if (t[ipp, 1] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          1], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 1])
                      if (t[ipp, 2] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          2], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 2])
                      tab$Freq[w] <- as.numeric(tab[w, nia] == 
                        value)
                    }
                    set.table(d, nia, tab)
                  }
                }
            }
        }
        if (compile) {
            compile(d)
        }
    }
    if ("AMEL" %in% mixture$markers) {
        ptypedgts <- c(ptypedgts, 1)
        names(ptypedgts)[length(ptypedgts)] <- "AMEL"
        if (compile) {
            compile(mixture$domains[["AMEL"]])
        }
    }
    invisible(ptypedgts)
}
rpt.typed.child <-
function (mixture, aca, ind = 1) 
{
    if (is.data.frame(aca)) 
        aca <- gt2aca(mixture, aca)
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            nC <- aca[[m]]
            na <- length(nC)
            if (max(aca[[m]]) == 2) {
                am <- which(aca[[m]] == max(aca[[m]]))
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a] - (am == 1), 
                      1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a] - (am == a), 
                      pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                        0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, tab[, S2am1]), 1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                am <- ab[1]
                bm <- ab[2]
                add.node(d, "ibdyet", subtype = "boolean")
                tab <- get.table(d, "ibdyet")
                tab$Freq <- c(0.5, 0.5)
                set.table(d, "ibdyet", tab, type = "cpt")
                add.edge(d, paste("n_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("n_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == am) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                          0, s, 1))
                    }
                  }
                  else if (a == bm) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], 1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], pmax(0, 1 - tab[, S2am1]), 
                        q[a]/ifelse(s > 0, s, 1))
                    }
                  }
                  else if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a], 1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a], pmax(0, 1 - 
                      tab[, S2am1]), q[a]/ifelse(s > 0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                  if (a == bm) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.parents <-
function (mixture, Mgt, Fgt, ind = 1, compile = TRUE) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Mgt$marker && m %in% Fgt$marker) {
            z <- get.nodes(d)
            for (n in z[substring(z, 1, 3 + floor(log10(ind))) == 
                paste("S", ind, sep = "_")]) delete.node(d, n)
            na <- nrow(mixture$data[[m]])
            ma1 <- match(Mgt$allele1[Mgt$marker == m], mixture$data[[m]]$allele)
            ma2 <- match(Mgt$allele2[Mgt$marker == m], mixture$data[[m]]$allele)
            fa1 <- match(Fgt$allele1[Fgt$marker == m], mixture$data[[m]]$allele)
            fa2 <- match(Fgt$allele2[Fgt$marker == m], mixture$data[[m]]$allele)
            add.node(d, "patt", states = 1:4)
            tab <- get.table(d, "patt")
            tab$Freq <- rep(0.25, 4)
            set.table(d, "patt", tab, type = "cpt")
            for (a in 1:na) {
                nia <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, nia, "patt")
                tab <- get.table(d, nia)
                sum <- ifelse(tab[, "patt"] == 1, (a == ma1) + 
                  (a == fa1), 0) + ifelse(tab[, "patt"] == 2, 
                  (a == ma2) + (a == fa1), 0) + ifelse(tab[, 
                  "patt"] == 3, (a == ma1) + (a == fa2), 0) + 
                  ifelse(tab[, "patt"] == 4, (a == ma2) + (a == 
                    fa2), 0)
                tab$Freq <- ifelse((tab[, nia] == sum), 1, 0)
                set.table(d, nia, tab, type = "cpt")
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.relative <-
function (mixture, Rgt, kinship = c(0.25, 0.5, 0.25), ind = 1, 
    compile = TRUE) 
{
    if (is.list(kinship)) {
        if (ncol(kinship$patt) != 4) 
            stop("can only specify 2-person pattern")
        if (any(apply(matrix(t(kinship$patt), nrow = 2), 2, diff) == 
            0)) 
            stop("this function does not handle autozygosity")
        kappas <- rep(0, 3)
        for (j in 1:nrow(kinship$patt)) {
            patt <- kinship$patt[j, ]
            k <- 1 + sum(outer(patt[1:2], patt[3:4], "=="))
            kappas[k] <- kappas[k] + kinship$pr[j]
        }
        kinship <- kappas
    }
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Rgt$marker) {
            z <- get.nodes(d)
            for (n in z[substring(z, 1, 3 + floor(log10(ind))) == 
                paste("S", ind, sep = "_")]) delete.node(d, n)
            ra1 <- match(Rgt$allele1[Rgt$marker == m], mixture$data[[m]]$allele)
            ra2 <- match(Rgt$allele2[Rgt$marker == m], mixture$data[[m]]$allele)
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            qstar <- q/rev(cumsum(rev(q)))
            na <- length(q)
            add.node(d, "patt", states = 1:4)
            tab <- get.table(d, "patt")
            tab$Freq <- c(kinship[1], kinship[2]/2, kinship[2]/2, 
                kinship[3])
            set.table(d, "patt", tab, type = "cpt")
            for (a in 1:na) {
                nsia <- paste("ns_", ind, "_", a, sep = "")
                add.node(d, nsia, states = 0:2)
                add.edge(d, nsia, "patt")
                if (a != 1) {
                  Siam1 <- ifelse(a == 2, paste("ns_", ind, "_1", 
                    sep = ""), paste("S_", ind, "_", a - 1, sep = ""))
                  add.edge(d, nsia, Siam1)
                  tab <- get.table(d, nsia)
                  tsofar <- tab[, Siam1]
                }
                else {
                  tab <- get.table(d, nsia)
                  tsofar <- rep(0, nrow(tab))
                }
                trials <- pmax(0, c(2, 1, 1, 0)[tab[, "patt"]] - 
                  tsofar)
                tab$Freq <- dbinom(tab[, nsia], trials, qstar[a])
                set.table(d, nsia, tab, type = "cpt")
                if (a != 1 & a != na) {
                  Siam1 <- ifelse(a == 2, paste("ns_", ind, "_1", 
                    sep = ""), paste("S_", ind, "_", a - 1, sep = ""))
                  Sia <- paste("S_", ind, "_", a, sep = "")
                  add.node(d, Sia, states = 0:2)
                  add.edge(d, Sia, nsia)
                  add.edge(d, Sia, Siam1)
                  tab <- get.table(d, Sia)
                  sum <- tab[, nsia] + tab[, Siam1]
                  tab$Freq <- ifelse((tab[, Sia] == sum) | (sum > 
                    2), 1, 0)
                  set.table(d, Sia, tab, type = "cpt")
                }
                nia <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, nia, nsia)
                add.edge(d, nia, "patt")
                tab <- get.table(d, nia)
                sum <- ifelse(tab[, "patt"] == 1, tab[, nsia], 
                  0) + ifelse(tab[, "patt"] == 2, tab[, nsia] + 
                  (a == ra1), 0) + ifelse(tab[, "patt"] == 3, 
                  tab[, nsia] + (a == ra2), 0) + ifelse(tab[, 
                  "patt"] == 4, tab[, nsia] + (a == ra1) + (a == 
                  ra2), 0)
                tab$Freq <- ifelse((tab[, nia] == sum) | (sum > 
                  2), 1, 0)
                set.table(d, nia, tab, type = "cpt")
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.relatives <-
function (mixture, IBD, typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs, quiet = FALSE, all.freq = NULL, 
    compile = TRUE) 
{
    if (!is.null(targets)) {
        jcontr <- match(contribs, targets)
	  inds<-which(!is.na(jcontr))
        if (!quiet) 
            cat("inds", inds, "\n")
        jcontr <- jcontr[!is.na(jcontr)]
        if (!quiet) 
            cat("jcontr", jcontr, "\n")
        jtyped <- match(names(typed.gts), targets)
        if (!quiet) 
            cat("jtyped", jtyped, "\n")
    }
    IBD <- convertIBD(IBD)
    if (length(jtyped) != length(typed.gts)) 
        stop("jtyped and typed.gts incompatible")
    if (length(jcontr) != length(inds)) 
        stop("jcontr and inds incompatible")
    if (!all(jtyped %in% (1:ncol(IBD$patt)/2)) || !all(jcontr %in% 
        (1:ncol(IBD$patt)/2))) 
        stop("jtyped, jcontr and IBD incompatible")
    ntyped <- length(jtyped)
    ap <- rep(0, max(IBD$patt))
    ptypedgts <- NULL
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (all(unlist(lapply(typed.gts, function(x) {
            m %in% x$marker
        })))) {
            keeplabels <- FALSE
            if (missing(all.freq) || is.null(all.freq)) {
                q <- mixture$data[[m]]$freq
                q <- q/sum(q)
            }
            else {
                if (is.data.frame(all.freq)) {
                  dbm <- subset(all.freq, marker == m)
                  q <- dbm$frequency[match(mixture$data[[m]]$allele, 
                    dbm$allele)]
                  q <- q/sum(q)
                }
                else {
                  keeplabels <- TRUE
                  q <- matrix(0, max(IBD$patt), nrow(mixture$data[[m]]))
                  for (k in 1:max(IBD$patt)) {
                    dbm <- subset(all.freq[[k]], marker == m)
                    q[k, ] <- dbm$frequency[match(mixture$data[[m]]$allele, 
                      dbm$allele)]
                    q[k, ] <- q[k, ]/sum(q[k, ])
                  }
                }
            }
            igt <- NULL
            if (ntyped > 0) 
                for (l in 1:ntyped) {
                  gt <- typed.gts[[l]]
                  a1 <- match(gt$allele1[gt$marker == m], mixture$data[[m]]$allele)
                  a2 <- match(gt$allele2[gt$marker == m], mixture$data[[m]]$allele)
                  igt <- c(igt, a1, a2)
                }
            z <- process.patterns(IBD, igt, jtyped, jcontr, q, 
                keeplabels)
            if (is.null(z)) {
                ptypedgts <- c(ptypedgts, 0)
                names(ptypedgts)[length(ptypedgts)] <- m
            }
            else {
                nodes <- get.nodes(d)
                for (i in inds) {
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("n", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                  }
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("S", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                    delete.node(d, n)
                  }
                }
                add.node(d, "pattperm", states = 1:attr(z, 
                  "npp"))
                tab <- get.table(d, "pattperm")
                f <- as.vector(z$pdenom)
                ptypedgts <- c(ptypedgts, sum(f))
                names(ptypedgts)[length(ptypedgts)] <- m
                tab$Freq <- f/sum(f)
                set.table(d, "pattperm", tab)
                na <- nrow(mixture$data[[m]])
                for (k in attr(z, "todraw")) {
                  if (is.null(dim(q))) {
                    qstar <- q/rev(cumsum(rev(q)))
                  }
                  else {
                    qstar <- q[k, ]/rev(cumsum(rev(q[k, ])))
                  }
                  for (a in 1:na) {
                    npka <- paste("np", k, "_", a, 
                      sep = "")
                    add.node(d, npka, states = 0:1)
                    if (a == 1) {
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[1], qstar[1])
                    }
                    else {
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.edge(d, npka, Spka1)
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[a], qstar[a], 1, 
                        0)
                    }
                    set.table(d, npka, tab)
                    if (a != 1 & a != na) {
                      Spka <- paste("Sp", k, "_", 
                        a, sep = "")
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.node(d, Spka, states = 0:1)
                      add.edge(d, Spka, npka)
                      add.edge(d, Spka, Spka1)
                      tab <- get.table(d, Spka)
                      tab$Freq <- as.numeric(tab[, Spka] == pmax(tab[, 
                        npka], tab[, Spka1]))
                      set.table(d, Spka, tab)
                    }
                  }
                }
                for (l in 1:length(inds)) {
                  i <- inds[l]
                  j <- l
                  t <- cbind(z[[paste("type.", 2 * j - 
                    1, sep = "")]], z[[paste("type.", 
                    2 * j, sep = "")]])
                  npj <- cbind(paste("np", z[[paste("oval.", 
                    2 * j - 1, sep = "")]], sep = ""), 
                    paste("np", z[[paste("oval.", 
                      2 * j, sep = "")]], sep = ""))
                  ov <- cbind(z[[paste("oval.", 2 * j - 
                    1, sep = "")]], z[[paste("oval.", 
                    2 * j, sep = "")]])
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    add.edge(d, nia, "pattperm")
                    for (ipp in 1:attr(z, "npp")) {
                      if (t[ipp, 1] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 1], 
                          a, sep = "_"))
                      if (t[ipp, 2] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 2], 
                          a, sep = "_"))
                    }
                  }
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    tab <- get.table(d, nia)
                    for (ipp in 1:attr(z, "npp")) {
                      w <- tab[, "pattperm"] == ipp
                      value <- rep(0, sum(w))
                      if (t[ipp, 1] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          1], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 1])
                      if (t[ipp, 2] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          2], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 2])
                      tab$Freq[w] <- as.numeric(tab[w, nia] == 
                        value)
                    }
                    set.table(d, nia, tab)
                  }
                }
            }
        }
        if (compile) {
            compile(d)
        }
    }
    if ("AMEL" %in% mixture$markers) {
        ptypedgts <- c(ptypedgts, 1)
        names(ptypedgts)[length(ptypedgts)] <- "AMEL"
        if (compile) {
            compile(mixture$domains[["AMEL"]])
        }
    }
    invisible(ptypedgts)
}
rpt.UAF <-
function (mixture, M, compile = TRUE) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        na <- nrow(mixture$data[[m]])
        q <- mixture$data[[m]]$freq
        q <- q/sum(q)
        alpha <- M * q
        beta <- c(rev(cumsum(rev(alpha)))[-1], 0)
        z <- get.nodes(d)
        n.unknown <- length(grep("^n_.*_1$", z))
        for (n in z[substring(z, 1, 1) == "n"]) for (par in get.parents(d, 
            n)) delete.edge(d, n, par)
        n <- outer(1:n.unknown, 1:na, function(i, a) paste("n", 
            i, a, sep = "_"))
        S <- outer(1:n.unknown, 1:na, function(i, a) paste("S", 
            i, a, sep = "_"))
        T <- outer(1:n.unknown, 1:na, function(i, a) paste("T", 
            i, a, sep = "_"))
        U <- outer(1:n.unknown, 1:na, function(i, a) paste("U", 
            i, a, sep = "_"))
        for (a in 1:na) {
            if (n.unknown > 1) 
                for (i in 1:(n.unknown - 1)) {
                  add.node(d, T[i, a], states = 0:(2 * i), subtype = "numbered")
                  add.node(d, U[i, a], states = 0:(2 * i), subtype = "numbered")
                }
            for (i in 1:n.unknown) {
                if (a == 1) {
                  if (i == 1) {
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), rep(alpha[1], lnval), rep(beta[1], 
                      lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), alpha[a] + tab[, T[i - 1, a]], 
                      (beta[a] + pmax(0, tab[, U[i - 1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                else {
                  if (i == 1) {
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], rep(alpha[a], lnval), 
                      rep(beta[a], lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], alpha[a] + tab[, T[i - 
                      1, a]], (beta[a] + pmax(0, tab[, U[i - 
                      1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, T[i, a], n[i, a])
                    tab <- get.table(d, T[i, a])
                    tab$Freq <- ifelse(tab[, T[i, a]] == tab[, 
                      n[i, a]], 1, 0)
                    set.table(d, T[i, a], tab)
                  }
                  else {
                    add.edge(d, T[i, a], n[i, a])
                    add.edge(d, T[i, a], T[i - 1, a])
                    tab <- get.table(d, T[i, a])
                    sumvals <- tab[, n[i, a]] + tab[, T[i - 1, 
                      a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      T[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, T[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, U[i, a], S[i, a])
                    tab <- get.table(d, U[i, a])
                    tab$Freq <- ifelse(tab[, U[i, a]] == 2 - 
                      tab[, S[i, a]], 1, 0)
                    set.table(d, U[i, a], tab)
                  }
                  else {
                    add.edge(d, U[i, a], S[i, a])
                    add.edge(d, U[i, a], U[i - 1, a])
                    tab <- get.table(d, U[i, a])
                    sumvals <- 2 - tab[, S[i, a]] + tab[, U[i - 
                      1, a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      U[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, U[i, a], tab)
                  }
                }
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
size <-
function (mixture) 
{
    res <- sum(unlist(lapply(mixture$dom, function(d) sum(unlist(lapply(summary(d, 
        jt = TRUE)$jt, function(j) j$size))))))
    class(res) <- "tablesize"
    res
}
tryCatch.W.E <-
function (expr) 
{
    W <- NULL
    w.handler <- function(w) {
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), 
        warning = w.handler), warning = W)
}
wlr <-
function (sep, Cgt, db, ind = 1, Mgt = NULL) 
{
    LR <- 1
    for (m in Cgt$marker) {
        sepm <- summary(sep)[[m]]
        pa1 <- sepm[[paste("U", ind, ".1", sep = "")]]
        pa2 <- sepm[[paste("U", ind, ".2", sep = "")]]
        Prob <- sepm$Prob
        if (is.null(Mgt)) {
            ca1 <- Cgt[Cgt$marker == m, ]$allele1
            ca2 <- Cgt[Cgt$marker == m, ]$allele2
            dbm <- db[db$marker == m, ]
            qca1 <- dbm$freq[dbm$allele == ca1]
            qca2 <- dbm$freq[dbm$allele == ca2]
            Den <- qca1 * qca2
            if (ca1 != ca2) 
                Den <- 2 * Den
            if (ca1 == ca2) {
                Num <- 0.5 * qca1 * ((pa1 == ca1) + (pa2 == ca1))
            }
            else {
                Num <- 0.5 * qca1 * ((pa1 == ca2) + (pa2 == ca2)) + 
                  0.5 * qca2 * ((pa1 == ca1) + (pa2 == ca1))
            }
            LRmfgt <- Num/Den
        }
        else {
            ca1 <- Cgt[Cgt$marker == m, ]$allele1
            ca2 <- Cgt[Cgt$marker == m, ]$allele2
            ma1 <- Mgt[Mgt$marker == m, ]$allele1
            ma2 <- Mgt[Mgt$marker == m, ]$allele2
            w <- match(ca1, c(ma1, ma2))
            if (is.na(w)) {
                w <- match(ca2, c(ma1, ma2))
                if (!is.na(w)) {
                  t <- ca1
                  ca1 <- ca2
                  ca2 <- t
                }
            }
            if (is.na(w)) {
                LRmfgt <- 0
            }
            else {
                if (w == 2) {
                  t <- ma1
                  ma1 <- ma2
                  ma2 <- t
                }
                dbm <- db[db$marker == m, ]
                qca1 <- dbm$freq[dbm$allele == ca1]
                qca2 <- dbm$freq[dbm$allele == ca2]
                Num <- 0.5 * ((pa1 == ca2) + (pa2 == ca2))
                Den <- qca2
                if (ca2 == ma2 && ca1 != ca2) {
                  Num <- Num + 0.5 * ((pa1 == ca1) + (pa2 == 
                    ca1))
                  Den <- Den + qca1
                }
                LRmfgt <- Num/Den
            }
        }
        LRm <- sum(LRmfgt * Prob)/sum(Prob)
        LR <- LR * LRm
    }
    LR
}
write.db <-
function (dbname, file = "db.csv") 
{
    DB <- getDb(dbname)
    marker <- allele <- frequency <- NULL
    for (m in names(DB)[-(1:3)]) {
        f <- DB[, m]
        a <- DB[!is.na(f), "Allele"]
        f <- f[!is.na(f)]
        if (length(f) > 0) {
            marker <- c(marker, rep(m, length(f)))
            allele <- c(allele, a)
            frequency <- c(frequency, f)
        }
    }
    marker <- as.factor(marker)
    allele[allele == "X"] <- 0
    allele[allele == "Y"] <- 1
    db <- data.frame(marker, allele = as.numeric(allele), frequency = as.numeric(frequency))
    if (file != "") 
        write.csv(db, "db.csv", row.names = FALSE)
    else db
}
write.epg <-
function (res, file = "epg.csv", C = 0) 
{
    epg <- subset(res[, c("Marker", "Allele", "Height")], Height > 
        C)
    names(epg) <- tolower(names(epg))
    epg$allele[epg$allele == "X"] <- 0
    epg$allele[epg$allele == "Y"] <- 1
    epg$allele <- as.numeric(epg$allele)
    if (file != "") 
        write.csv(epg, file, row.names = FALSE)
    else epg
}
