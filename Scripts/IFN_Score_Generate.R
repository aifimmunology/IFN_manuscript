library(dplyr)

ifn_score_generate <- function(donor_df, 
                                cell_type_level, 
                                cell_type,
                                stim_mat_dir) {
  #' Generate cell type specific IFNa and IFNg response scores from subject-level DEG results
  #'
  #' @param donor_df Data frame containing differential expression results. 
  #'   First two columns must be gene names and log2FC values, respectively.
  #' @param cell_type_level Character string: "L1", "L2", or "PBMC".
  #'   - "L1": Use broad level 1 cell types
  #'   - "L2": Use more granular level 2 cell types
  #'   - "PBMC": If input differentials from bulk RNAseq, use PBMC-level 
  #' @param cell_type Character string indicating the specific cell type.
  #'   Must match valid choices for the selected `cell_type_level`.
  #'   Ignored if `cell_type_level` is "PBMC".
  #' @param stim_mat_dir Path to directory containing the IFNa/IFNg stimulation matrix `IFNa_IFNg_Scoring_ISGs_Mat.csv`.
  #'
  #' @return Data frame with columns:
  #'   - `IFNa_score`: numeric IFNa response score.
  #'   - `IFNg_score`: numeric IFNg response score.
  
  # Valid cell types for L1 and L2
  celltypes_l1 <- c("Tcell", "Bcell", "NK", "Monocyte")
  celltypes_l2 <- c("CD4 Naive", "CD8 Naive", "CD4 Memory", "CD8 Memory", 
                    "NK.CD56hi", "NK.CD56dim", "Treg", "gdT", "MAIT", 
                    "B Naive", "B Memory", "Plasma")
  
  # Cell types that do not respond to IFNg
  IFNg_nonresponding <- c("NK.CD56hi", "NK.CD56dim", "gdT", "MAIT", 
                          "Memory CD8", "Plasma")
  
  if (cell_type_level == "PBMC") {
    cell_type <- "PBMC"
  } else if (cell_type_level == "L1") {
    if (!cell_type %in% celltypes_l1) {
      stop("For L1, cell_type must be one of: ", paste(celltypes_l1, collapse = ", "))
    }
  } else if (cell_type_level == "L2") {
    if (!cell_type %in% celltypes_l2) {
      stop("For L2, cell_type must be one of: ", paste(celltypes_l2, collapse = ", "))
    }
  } else {
    stop('cell_type_level must be "L1", "L2", or "PBMC"')
  }

  # Load IFN response matrix and extract cell type specific signature
  isg_df <- data.table::fread(file.path(stim_mat_dir, "IFNa_IFNg_Scoring_ISGs_Mat.csv")) %>%
                dplyr::filter(celltype_level == cell_type_level &
                              subset == cell_type) %>%
                dplyr::select(Gene, IFNa, IFNg) %>%
                dplyr::mutate(Gene = as.character(Gene))
    
  # Ensure donor_df has correct column names and data type
  colnames(donor_df)[1:2] <- c("Gene", "Log2FC")


  if (!is.character(donor_df$Gene) && !is.factor(donor_df$Gene)) {
    stop("First column (Gene) must be of type character.")
  }
  if (!is.numeric(donor_df$Log2FC)) {
    stop("Second column (Log2FC) must be numeric.")
  }
    
  donor_df <- donor_df %>% 
                dplyr::select(Gene, Log2FC) %>%
                dplyr::mutate(Gene = as.character(Gene))
    
  # Join with ISG matrix, remove genes not in user results, and scale 
  nnls_mat <- dplyr::left_join(donor_df, isg_df, by = "Gene") %>%
    tidyr::drop_na() %>%
    tibble::column_to_rownames("Gene") %>%
    as.matrix() %>%
    scale(center = FALSE)
  
  # Run NNLS
  nnls_res <- nnls::nnls(A = nnls_mat[, 2:3], b = nnls_mat[, 1])
  
  res_df <- data.frame(
    IFNa_score = nnls_res$x[1],
    IFNg_score = nnls_res$x[2]
  )
  
  # Set IFNg score to zero for non-responding cell types
  if (cell_type %in% IFNg_nonresponding) {
    res_df$IFNg_score <- 0
  }
  
  return(res_df)
}