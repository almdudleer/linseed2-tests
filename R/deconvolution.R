library(CellMix)
library(linseed)

library(limma)
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)

source("R/docker_linseed_pipeline/LinseedMetadata.R")
source("R/docker_linseed_pipeline/SinkhornNNLSLinseedC.R")
Rcpp::sourceCpp("R/docker_linseed_pipeline/pipeline.cpp")

bwr <- colorRampPalette(c("blue", "white", "red"))(50)

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
        top_genes = 493
        global_iterations = 2000
        coef_hinge_H <- 10
        coef_hinge_W <- 10
        coef_pos_D_h <- 0.01
        coef_pos_D_w <- 0.01
        coef_der_X <- 0.00001
        coef_der_Omega <- 0.0000001

        lo2 <- SinkhornNNLSLinseed$new(
            dataset = edata,
            path = "linseed2_results",
            analysis_name = "deconv_comparison_2",
            cell_types = ntypes,
            global_iterations = global_iterations,
            coef_der_X = coef_der_X,
            coef_der_Omega = coef_der_Omega,
            coef_pos_D_h = coef_pos_D_h,
            coef_pos_D_w = coef_pos_D_w,
            coef_hinge_H = coef_hinge_H,
            coef_hinge_W = coef_hinge_W
        )
        lo2$selectTopGenes(min(top_genes, nrow(lo2$filtered_data)))
        lo2$scaleDataset(iterations = 20)
        lo2$getSvdProjectionsNew(k = ntypes)
        lo2$selectInitOmega()
        lo2$runOptimization()
        lo2
    }
}
