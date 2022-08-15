LINSEED2_PATH <- "../../docker_linseed_pipeline/"
R_SRC_PATH <- "../R/"

source(paste0(R_SRC_PATH, "de.R"))
source(paste0(R_SRC_PATH, "deconvolution.R"))
source(paste0(R_SRC_PATH, "linseed_2_utils.R"))
source(paste0(R_SRC_PATH, "metrics.R"))
source(paste0(LINSEED2_PATH, "simulation.R"))
source(paste0(LINSEED2_PATH, "LinseedMetadata.R"))
source(paste0(LINSEED2_PATH, "SinkhornNNLSLinseedC.R"))
Rcpp::sourceCpp(paste0(LINSEED2_PATH, "pipeline.cpp"))