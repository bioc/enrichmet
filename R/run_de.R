#' Run Differential Expression Analysis using LIMMA
#'
#' Performs differential expression analysis between two sample groups using LIMMA,
#' generates volcano plots, and prepares results for downstream enrichment analysis.
#'
#' @param df A data frame or matrix with metabolites as rows and samples as columns. 
#'           Rownames should be metabolite identifiers.
#' @param group1_pattern Character string pattern to identify samples in group 1
#' @param group2_pattern Character string pattern to identify samples in group 2
#' @param plot_volcano Logical indicating whether to generate volcano plot (default: TRUE)
#' @param top_n_labels Number of top features to label in volcano plot (default: 20)
#' @param fc_threshold Fold change threshold for significance (default: 1)
#' @param pval_threshold P-value threshold for significance (default = 0.05)
#' @param min_samples Minimum number of samples required per group (default = 3)
#' @param split_complex_ids Logical indicating whether to split complex KEGG IDs like 
#'        "C00025|C00979" into individual metabolites (default = TRUE)
#'
#' @return A list containing:
#'   \item{full_results}{Full results with all metabolites}
#'   \item{kegg_ready}{Results filtered for KEGG metabolites only}
#'   \item{volcano_plot}{ggplot object of volcano plot (if plot_volcano = TRUE)}
#'   \item{summary_stats}{Summary statistics of the analysis}
#'
#' @examples
#' 
#' # Example showing clear differential expression patterns
#' set.seed(123)
#' 
#' # Create data with 30 metabolites and 20 samples
#' n_metabolites <- 30
#' n_samples <- 20
#' simulated_data <- matrix(0, nrow = n_metabolites, ncol = n_samples)
#' 
#' # Control samples (columns 1-10)
#' simulated_data[, 1:10] <- rnorm(n_metabolites * 10, mean = 10, sd = 2)
#' 
#' # Treatment samples (columns 11-20) with clear patterns
#' simulated_data[1:5, 11:20] <- rnorm(5 * 10, mean = 15, sd = 2)   # Upregulated
#' simulated_data[6:10, 11:20] <- rnorm(5 * 10, mean = 5, sd = 2)    # Downregulated
#' simulated_data[11:30, 11:20] <- rnorm(20 * 10, mean = 10, sd = 2) # Not significant
#' 
#' # Set row and column names
#' rownames(simulated_data) <- paste0("Met_", 1:30)
#' colnames(simulated_data) <- c(paste0("Control_", 1:10), paste0("Treatment_", 1:10))
#' 
#' # Run differential expression analysis
#' results <- run_de(
#'   df = simulated_data,
#'   group1_pattern = "Control",
#'   group2_pattern = "Treatment",
#'   plot_volcano = TRUE,
#'   fc_threshold = 1,
#'   pval_threshold = 0.05
#' )
#' 
#' # View significant metabolites
#' sig_metabolites <- results$full_results[results$full_results$Significant != "Not significant", ]
#' head(sig_metabolites)
#' 
#' # View summary statistics
#' results$summary_stats
#' 
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline theme_bw xlab ylab theme ggtitle
#' @importFrom ggrepel geom_text_repel
#' @export
run_de <- function(df, group1_pattern, group2_pattern, plot_volcano = TRUE, 
                   top_n_labels = 20, fc_threshold = 1, pval_threshold = 0.05,
                   min_samples = 3, split_complex_ids = FALSE) {   # default changed to FALSE
    
    # Input validation
    if (!is.data.frame(df) && !is.matrix(df)) {
        stop("Input must be a data frame or matrix")
    }
    
    # Convert to matrix if data frame
    if (is.data.frame(df)) {
        expr <- as.matrix(df)
    } else {
        expr <- df
    }
    
    # Detect sample sets
    g1_cols <- grep(group1_pattern, colnames(expr), value = TRUE)
    g2_cols <- grep(group2_pattern, colnames(expr), value = TRUE)
    
    if (length(g1_cols) == 0 || length(g2_cols) == 0) {
        stop("Pattern did not match any samples in column names")
    }
    
    if (length(g1_cols) < min_samples || length(g2_cols) < min_samples) {
        warning("Fewer than ", min_samples, 
                " samples in one or both groups. Results may be unreliable.")
    }
    
    # Subset expression matrix
    expr_subset <- expr[, c(g1_cols, g2_cols), drop = FALSE]
    
    # Group factor
    group <- factor(
        ifelse(colnames(expr_subset) %in% g1_cols, "G1", "G2"), 
        levels = c("G1", "G2")
    )
    
    # Design and fit
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    # Check for constant rows
    constant_rows <- apply(expr_subset, 1, function(x) {
        sd_val <- stats::sd(x, na.rm = TRUE)
        !is.na(sd_val) && sd_val == 0
    })
    
    if (any(constant_rows)) {
        warning("Removing ", sum(constant_rows), 
                " metabolites with zero variance")
        expr_subset <- expr_subset[!constant_rows, , drop = FALSE]
    }
    
    fit <- limma::lmFit(expr_subset, design)
    contrast.matrix <- limma::makeContrasts(G1vsG2 = G1 - G2, levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    
    # Results
    res <- limma::topTable(
        fit2, 
        coef = "G1vsG2", 
        number = nrow(expr_subset), 
        sort.by = "none"
    )
    
    res_df <- res %>%
        dplyr::mutate(
            met_id = rownames(res),
            pval = .data$P.Value,
            padj = .data$adj.P.Val,
            log2fc = .data$logFC,
            ave_expr = .data$AveExpr
        ) %>%
        dplyr::select(.data$met_id, .data$log2fc, .data$pval, 
                      .data$padj, .data$ave_expr)
    
    # Add significance column
    res_df$Significant <- "Not significant"
    res_df$Significant[res_df$padj < pval_threshold & 
                           res_df$log2fc > fc_threshold] <- "Up"
    res_df$Significant[res_df$padj < pval_threshold & 
                           res_df$log2fc < -fc_threshold] <- "Down"
    
    # Create volcano plot
    volcano_obj <- NULL
    if (plot_volcano) {
        plot_data <- res_df
        plot_data$label <- ifelse(
            plot_data$padj < pval_threshold & 
                abs(plot_data$log2fc) > fc_threshold,
            plot_data$met_id, 
            ""
        )
        
        significant_features <- plot_data[
            plot_data$Significant != "Not significant", 
        ]
        non_significant_features <- plot_data[
            plot_data$Significant == "Not significant", 
        ]
        
        if (nrow(significant_features) > 0) {
            significant_features <- significant_features[
                order(significant_features$padj), 
            ]
        }
        if (nrow(non_significant_features) > 0) {
            non_significant_features <- non_significant_features[
                order(non_significant_features$padj), 
            ]
        }
        
        if (nrow(significant_features) >= top_n_labels) {
            top_features <- significant_features[
                seq_len(min(top_n_labels, nrow(significant_features))), 
            ]
        } else {
            top_features <- rbind(
                significant_features,
                non_significant_features[
                    seq_len(min(
                        top_n_labels - nrow(significant_features), 
                        nrow(non_significant_features)
                    )), 
                ]
            )
        }
        
        top_features <- top_features[
            stats::complete.cases(top_features), 
        ]
        
        up_count <- sum(plot_data$Significant == "Up")
        down_count <- sum(plot_data$Significant == "Down")
        total_significant <- up_count + down_count
        
        volcano_obj <- ggplot2::ggplot(
            plot_data, 
            ggplot2::aes(x = .data$log2fc, y = -log10(.data$padj))
        ) +
            ggplot2::geom_point(
                ggplot2::aes(color = .data$Significant, 
                             alpha = .data$Significant), 
                size = 2
            ) +
            ggplot2::geom_vline(
                xintercept = c(-fc_threshold, fc_threshold), 
                linetype = "dashed", 
                color = "black", 
                alpha = 0.7, 
                linewidth = 0.8
            ) +
            ggplot2::geom_hline(
                yintercept = -log10(pval_threshold), 
                linetype = "dashed", 
                color = "black", 
                alpha = 0.7, 
                linewidth = 0.8
            ) +
            ggrepel::geom_text_repel(
                data = top_features, 
                ggplot2::aes(label = .data$met_id),
                size = 3, 
                max.overlaps = Inf,
                box.padding = 0.5, 
                point.padding = 0.2,
                min.segment.length = 0.2,
                segment.color = "grey50",
                segment.alpha = 0.5
            ) +
            ggplot2::scale_color_manual(
                values = c(
                    "Up" = "#E41A1C", 
                    "Down" = "#377EB8", 
                    "Not significant" = "grey70"
                )
            ) +
            ggplot2::scale_alpha_manual(
                values = c(
                    "Up" = 0.8, 
                    "Down" = 0.8, 
                    "Not significant" = 0.4
                )
            ) +
            ggplot2::theme_bw(base_size = 14) +
            ggplot2::xlab("log2 Fold Change") +
            ggplot2::ylab("-log10 FDR") +
            ggplot2::theme(
                legend.title = ggplot2::element_blank(),
                legend.position = "right",
                panel.grid.minor = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5)
            ) +
            ggplot2::ggtitle(
                paste0("Volcano Plot: ", group1_pattern, " vs ", group2_pattern),
                subtitle = paste0(
                    total_significant, " significant metabolites (", 
                    up_count, " up, ", down_count, " down) | ",
                    "FC threshold: +/-", fc_threshold, 
                    ", FDR threshold: ", pval_threshold
                )
            )
        
        message("Volcano plot created successfully")
        message("Total significant metabolites: ", total_significant)
        message("  Up-regulated: ", up_count)
        message("  Down-regulated: ", down_count)
    }
    
    # ------------------------------------------------------------------
    # FIXED: KEGG ID extraction – always keep ONLY the first valid ID
    # (consistent with perform_gsea_analysis). Never duplicate statistics.
    # ------------------------------------------------------------------
    res_kegg <- res_df
    
    # Helper that extracts the first valid KEGG ID (same logic as GSEA)
    extract_first_kegg <- function(met_id) {
        if (is.na(met_id) || met_id == "") return(NA_character_)
        
        # Already clean
        if (grepl("^C\\d{5}$", met_id)) return(met_id)
        
        # Complex formats: take everything after last underscore, then first ID
        if (grepl("_", met_id)) {
            last_part <- sub(".*_", "", met_id)
            candidates <- unlist(strsplit(last_part, "\\|"))
        } else {
            candidates <- unlist(strsplit(met_id, "\\|"))
        }
        
        candidates <- trimws(candidates)
        valid <- candidates[grepl("^C\\d+", candidates)]
        
        if (length(valid) > 0) {
            return(valid[1])          # ALWAYS first only
        }
        return(NA_character_)
    }
    
    res_kegg$kegg_id <- vapply(res_kegg$met_id, extract_first_kegg, 
                               FUN.VALUE = character(1))
    
    # Warn if complex IDs were present
    complex_present <- any(grepl("\\|", res_df$met_id), na.rm = TRUE)
    if (complex_present) {
        warning(
            "Complex KEGG IDs (containing '|') detected. ",
            "Only the FIRST valid KEGG ID is retained for each metabolite ",
            "to preserve independence assumptions of Fisher's exact test and GSEA. ",
            "Secondary IDs are discarded."
        )
    }
    
    # Remove rows without a valid KEGG ID
    res_kegg <- res_kegg[!is.na(res_kegg$kegg_id) & res_kegg$kegg_id != "", ]
    rownames(res_kegg) <- NULL
    
    # Optional: if user still sets split_complex_ids = TRUE we now ignore it
    # and issue a message (parameter kept for backward compatibility)
    if (isTRUE(split_complex_ids)) {
        message(
            "Note: split_complex_ids = TRUE is ignored. ",
            "Duplication of statistics is no longer performed ",
            "(it violated independence assumptions). ",
            "Only the first KEGG ID is kept."
        )
    }
    
    # Summary statistics
    summary_stats <- list(
        total_metabolites = nrow(res_df),
        significant_metabolites = sum(res_df$Significant != "Not significant"),
        upregulated = sum(res_df$Significant == "Up"),
        downregulated = sum(res_df$Significant == "Down"),
        group1_samples = length(g1_cols),
        group2_samples = length(g2_cols),
        fc_threshold = fc_threshold,
        pval_threshold = pval_threshold,
        split_complex_ids = FALSE,          # always FALSE now
        unique_kegg_ids = length(unique(res_kegg$kegg_id))
    )
    
    # Return comprehensive results
    results_list <- list(
        full_results = res_df,
        kegg_ready = res_kegg,
        summary_stats = summary_stats
    )
    
    if (plot_volcano) {
        results_list$volcano_plot <- volcano_obj
    }
    
    message("Differential analysis completed successfully")
    message("Total metabolites analyzed: ", nrow(res_df))
    message("Significant metabolites (FDR < ", pval_threshold, 
            ", |FC| > ", fc_threshold, "): ", 
            sum(res_df$Significant != "Not significant"))
    
    return(results_list)
}