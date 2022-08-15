library(CellMix)
library(linseed)

library(limma)
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)

bwr <- colorRampPalette(c("blue", "white", "red"))(50)

addDeMarkerList <- function(simulation) {
    stopifnot(!is.null(simulation$pure_sample_names))
    pure_es <- ExpressionSet(
        assayData = simulation$data[, simulation$pure_sample_names],
        phenoData = AnnotatedDataFrame(as.data.frame(t(simulation$proportions[, simulation$pure_sample_names])))
    )
    simulation$marker_list <- getDeMarkerList(pure_es)
    simulation
}

barplotMetric <- function(predicts, truth, metric.f, name, colors) {
    metrics <- predicts %>%
        purrr::map(
            ~ sapply(seq_len(nrow(.)), function(i) metric.f(.[i, ], truth[i, ])) %>% `names<-`(rownames(truth))
        ) %>%
        as.data.frame() %>%
        as.matrix()
    par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    barplot(metrics, main = name, beside = TRUE, col = colors)
    legend("topright", inset = c(-0.2, 0), legend = rownames(metrics), fill = colors)
}

guessOrder <- function(predicted, actual) {
    ctn <- nrow(predicted)
    allPerms <- combinat::permn(ctn)

    vals <- sapply(allPerms, function(perm) {
        sum(abs(predicted[perm, ] - actual))
    })
    perm <- allPerms[[which.min(vals)]]
    perm
}

meltPredictedAndActualBasis <- function(predicted, actual) {
    dp <- as.data.frame(predicted)
    da <- as.data.frame(actual)
    dp$gene <- rownames(dp)
    da$gene <- rownames(da)
    dpm <- melt(dp, id.vars = "gene", variable.name = "tissue", value.name = "predicted")
    dam <- melt(da, id.vars = "gene", variable.name = "tissue", value.name = "actual")
    merge(dpm, dam, by = c("gene", "tissue"))
}

getMeltedBasis <- function(mix, fits, lo) {
    basis.true <- basis(mix)
    if (quantile(basis.true, 0.99) > 1000) {
        basis.true <- log2(basis.true + 1)
    }
    colnames(basis.true) <- tolower(colnames(basis.true))
    basis.preds <- fits %>% purrr::map(~ basis(.))
    basis.preds$linseed <- lo$signatures
    basis.preds <- basis.preds %>% purrr::map(~ log2(. + 1))
    basis.melted <- list()
    for (name in names(basis.preds)) {
        bp <- basis.preds[[name]]
        colnames(bp) <- colnames(basis.true)[
            guessOrder(
                t(basis.true[rownames(bp), ]),
                t(bp)
            )
        ]
        basis.melted[[name]] <- meltPredictedAndActualBasis(bp, basis.true)
        basis.melted[[name]]$model <- name
    }
    basis.all <- do.call("rbind", basis.melted)
    rownames(basis.all) <- NULL
    basis.all
}

plotMeltedBasis <- function(meltedBasis) {
    ggplot(meltedBasis, aes(actual, predicted)) +
        geom_point(shape = ".") +
        stat_cor(method = "pearson", aes(label = ..rr.label..), label.y.npc = 0.95, cor.coef.name = "R") +
        stat_cor(method = "spearman", aes(label = ..rr.label..), label.y.npc = 0.85, cor.coef.name = "rho") +
        facet_grid(cols = vars(tissue), rows = vars(model)) +
        theme(text = element_text(size = 20))
}

fitDeconvolutionModel <- function(edata, kind, ntypes, markers = NULL, verbose = FALSE) {
    if (kind == "dsa") {
        CellMix::ged(edata, x = markers, "DSA", verbose = TRUE)
    } else if (kind == "lee") {
        CellMix::ged(edata, x = ntypes, method = "ssFrobenius", data = markers, verbose = TRUE, sscale = FALSE)
    } else if (kind == "brunet") {
        CellMix::ged(edata, x = ntypes, method = "ssKL", verbose = TRUE, sscale = FALSE)
    } else if (kind == "deconf") {
        CellMix::ged(edata, x = ntypes, method = "deconf", verbose = TRUE)
    } else if (kind == "linseed") {
        lo <- linseed::LinseedObject$new(edata)
        lo$calculatePairwiseLinearity()
        lo$calculateSpearmanCorrelation()
        lo$calculateSignificanceLevel(100)
        lo$filterDatasetByPval(0.01)
        lo$setCellTypeNumber(ntypes)
        lo$project("full")
        lo$smartSearchCorners(dataset = "filtered", error = "norm")
        lo$deconvolveByEndpoints(dataset = "full", error = "raw")
        lo
    } else if (kind == "linseed2") {
        lo2 <- init_lo2(edata, min(nrow(edata), 10000))
        samples_flt <- do_mad_cutoff(lo2, samples = T)
        genes_flt <- do_mad_cutoff(lo2, samples = F)
        edata_flt <- edata[rownames(edata) %in% genes_flt, colnames(edata) %in% samples_flt]
        lo2 <- init_lo2(edata_flt, min(nrow(edata_flt), 10000))
        run_block(lo2, iterations = 10000)
        lo2
    }
}
