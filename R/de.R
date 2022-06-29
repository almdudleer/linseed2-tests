getDeMarkers <- function(tissue,
                         ps_data,
                         n_markers_each = 30,
                         adj_pval_threshold = 0.01,
                         stat = "logFC",
                         plot_de = T) {
    design <- model.matrix(as.formula(stringr::str_interp("~1 + ${tissue}")), data = Biobase::pData(ps_data))
    colnames(design) <- c("intercept", tissue)
    fit <- limma::lmFit(ps_data, design)
    fit2 <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = tissue, levels = design))
    fit2 <- limma::eBayes(fit2, trend = T)
    de <- limma::topTable(fit2, adjust.method = "BH", number = Inf, sort.by = stat)
    de <- de[de$logFC > 0, ]
    de_markers <- rownames(head(de, n_markers_each))
    de_markers <- rowMeans(Biobase::exprs(ps_data)[de_markers, Biobase::pData(ps_data)[tissue] != 0])
    if (plot_de) {
        pheatmap::pheatmap(
            ps_data[rownames(de), ],
            cluster_rows = F, cluster_cols = F,
            show_rownames = F, color = bwr
        )
    }
    de_markers
}

getDeMarkerList <- function(pure_es) {
    exprs(pure_es) <- limma::normalizeBetweenArrays(log2(exprs(pure_es) + 1), method = "quantile")
    ml_raw <- colnames(pData(pure_es)) %>%
        purrr::set_names() %>%
        purrr::map(~getDeMarkers(., pure_es, plot_de = F))
    marker_names <- ml_raw %>% purrr::map(names) %>% purrr::flatten()
    marker_names <- marker_names[!(duplicated(marker_names) | duplicated(marker_names, fromLast = T))]
    ml_raw <- ml_raw %>% purrr::map(~.[names(.) %in% marker_names])
    CellMix::MarkerList(ml_raw)
}
