source("R/de.R")
source("R/metrics.R")

# Logic
generateBasis <- function(n_genes, n_cell_types, sd_ = 0.2) {
    data_1 <- c(
        rnorm(n_genes * 2, mean = 4, sd = 0.75),
        rnorm(n_genes * 3, mean = 10, sd = 1.5)
    )
    basis <- matrix(0, nrow = n_genes, ncol = n_cell_types)
    sds_ <- rnorm(n_cell_types, mean = sd_, sd = 0.02)

    for (i in 1:n_genes) {
        basis[i, 1] <- data_1[sample(seq_len(length(data_1)), 1)] * rnorm(1, mean = 1, sd = sds_[1])
        for (j in 2:n_cell_types) {
            # DISCUSS: the only difference between the cell types is their deviation from the first type
            # which is seen on differential expression
            basis[i, j] <- basis[i, 1] * rnorm(1, mean = 1, sd = sds_[j])
        }
    }

    basis <- limma::normalizeBetweenArrays(basis, method = "quantile")
    basis <- 2^basis
    rownames(basis) <- paste0("gene_", 1:n_genes)
    colnames(basis) <- paste0("cell_type_", 1:n_cell_types)

    basis
}

generateProportions <- function(n_samples, n_cell_types) {
    proportions <- linseed:::sampleFromSimplexUniformly(n_samples, n_cell_types, 100000)
    colnames(proportions) <- paste0("sample_", 1:n_samples)
    rownames(proportions) <- paste0("cell_type_", 1:n_cell_types)

    proportions
}

addNoiseToData <- function(data, noise_deviation) {
    # DISCUSS: additive or multiplicative noise
    noise_mask <- matrix(rnorm(length(data), sd = noise_deviation), nrow = nrow(data), ncol = ncol(data))
    noisy_data <- data + 2^noise_mask
    noisy_data[noisy_data < 0] <- 0
    noisy_data
}

generatePureSamples <- function(basis, samples_per_cell_type, noise_deviation = 0.05) {
    data <- basis[, rep(seq_len(ncol(basis)), each = samples_per_cell_type)] * matrix(
        rnorm(nrow(basis) * ncol(basis) * samples_per_cell_type, mean = 1, sd = noise_deviation),
        nrow = nrow(basis),
        ncol = ncol(basis) * samples_per_cell_type
    )
    proportions <- diag(1, ncol(basis), ncol(basis))[rep(seq_len(ncol(basis)), each = samples_per_cell_type), ]
    colnames(proportions) <- colnames(basis)
    rownames(proportions) <- paste0("pure_sample_", seq_len(nrow(proportions)))
    proportions <- t(proportions)
    colnames(data) <- colnames(proportions)
    list(proportions = proportions, data = data)
}

generateBasisSamples <- function(basis) {
    proportions <- diag(1, ncol(basis), ncol(basis))
    rownames(proportions) <- colnames(basis)
    colnames(proportions) <- paste0("basis_sample_", seq_len(ncol(proportions)))
    data <- basis
    colnames(data) <- colnames(proportions)
    list(proportions = proportions, data = data)
}


# Builder
createSimulation <- function(n_genes, n_samples, n_cell_types) {
    basis <- generateBasis(n_genes, n_cell_types)
    proportions <- generateProportions(n_samples, n_cell_types)
    data <- basis %*% proportions

    data[data < 0] <- 0

    list(
        basis = basis,
        proportions = proportions,
        data = data,
        mixed_sample_names = colnames(data)
    )
}

withPureSamples <- function(simulation, samples_per_cell_type, noise_deviation = 0.05) {
    pure_samples <- generatePureSamples(
        simulation$basis,
        samples_per_cell_type,
        noise_deviation = noise_deviation
    )
    simulation$pure_sample_names <- colnames(pure_samples$data)
    simulation$data <- cbind(simulation$data, pure_samples$data)
    simulation$proportions <- cbind(simulation$proportions, pure_samples$proportions)
    simulation$pure_samples_noise_deviation <- noise_deviation
    simulation
}

withBasisSamples <- function(simulation) {
    basis_samples <- generateBasisSamples(simulation$basis)
    simulation$basis_sample_names <- colnames(basis_samples$data)
    simulation$data <- cbind(simulation$data, basis_samples$data)
    simulation$proportions <- cbind(simulation$proportions, basis_samples$proportions)
    simulation
}

withNoise <- function(simulation, noise_deviation) {
    simulation$data <- addNoiseToData(simulation$data, noise_deviation)
    simulation$data_noise_deviation <- noise_deviation
    simulation
}

withDeMarkerList <- function(simulation) {
    stopifnot(!is.null(simulation$pure_sample_names))
    pure_es <- ExpressionSet(
        assayData = simulation$data[, simulation$pure_sample_names],
        phenoData = AnnotatedDataFrame(as.data.frame(t(simulation$proportions[, simulation$pure_sample_names])))
    )
    simulation$marker_list <- getDeMarkerList(pure_es)
    simulation
}


# All together
generateData <- function(n_genes, n_samples, n_cell_types, noise_deviation = 0, pure_samples_per_cell_type = 0) {
    simulation <- createSimulation(n_genes, n_samples, n_cell_types)
    if (noise_deviation > 0) {
        simulation <- withNoise(simulation, noise_deviation)
    }
    if (pure_samples_per_cell_type > 0) {
        simulation <- withPureSamples(pure_samples_per_cell_type)
    }
    simulation
}
