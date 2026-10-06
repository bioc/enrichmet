#' Convert a SummarizedExperiment to a da_results Object
#'
#' Reconstructs the \code{da_results} structure expected by \code{enrichmet()}
#' from a \code{SummarizedExperiment}. Metabolite identifiers, fold changes,
#' and significance statistics are read from \code{rowData()}, so no external
#' \code{da_out} object is required.
#'
#' @param se A \code{SummarizedExperiment} (or an object that inherits from it)
#'   whose \code{rowData()} contains at least the columns \code{met_id},
#'   \code{log2fc}, \code{pval}, and \code{padj}.
#' @param fc_threshold Numeric. Absolute log2 fold change threshold used to
#'   label metabolites as \code{"Up"} or \code{"Down"} when the
#'   \code{Significant} column is not already present in \code{rowData()}.
#'   Default is \code{1}.
#' @param pval_threshold Numeric. Adjusted p value threshold used together with
#'   \code{fc_threshold} to assign significance labels when \code{Significant}
#'   is absent. Default is \code{0.05}.
#' @param split_complex_ids Logical. If \code{TRUE}, metabolite identifiers
#'   that embed multiple KEGG IDs (for example
#'   \code{neg_00096_C00025|C00979}) are split into one row per KEGG ID in the
#'   reconstructed \code{kegg_ready} table. If \code{FALSE}, only the first
#'   KEGG ID per metabolite is retained. Default is \code{TRUE}.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{full_results}{Data frame with columns \code{met_id},
#'     \code{log2fc}, \code{pval}, \code{padj}, optionally \code{ave_expr},
#'     and \code{Significant}.}
#'   \item{kegg_ready}{Data frame with the same columns as
#'     \code{full_results} plus a \code{kegg_id} column, one row per KEGG ID
#'     when \code{split_complex_ids = TRUE}.}
#'   \item{summary_stats}{List of summary counts and thresholds, mirroring
#'     the output of \code{run_de()}.}
#' }
#'
#' @details
#' \code{da_results_from_se()} performs the following steps:
#' \enumerate{
#'   \item Extracts \code{rowData(se)} as a data frame and verifies that
#'     \code{met_id}, \code{log2fc}, \code{pval}, and \code{padj} are present.
#'   \item Recomputes the \code{Significant} column if it is absent, using the
#'     supplied \code{fc_threshold} and \code{pval_threshold}, matching the
#'     conventions used by \code{run_de()}.
#'   \item Rebuilds the \code{kegg_ready} table by extracting KEGG IDs from
#'     \code{met_id}.
#'   \item Rebuilds the \code{summary_stats} list from the assay dimensions and
#'     the significance labels.
#' }
#' The helper does not require the original \code{da_out} object. Once the
#' differential statistics live in \code{rowData()}, the
#' \code{SummarizedExperiment} is self contained for the purposes of
#' \code{enrichmet()}.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' # Example abundance matrix
#' set.seed(123)
#' mat <- matrix(
#'     rnorm(100, mean = 10, sd = 2),
#'     nrow = 10,
#'     dimnames = list(
#'         paste0("Met_", 1:10),
#'         c(paste0("Control_", 1:5), paste0("Treat_", 1:5))
#'     )
#' )
#'
#' # Differential statistics stored in rowData()
#' rd <- DataFrame(
#'     met_id   = rownames(mat),
#'     log2fc   = rnorm(10),
#'     pval     = runif(10, 0, 0.1),
#'     padj     = runif(10, 0, 0.2),
#'     ave_expr = rowMeans(mat)
#' )
#'
#' se <- SummarizedExperiment(
#'     assays  = list(counts = mat),
#'     rowData = rd
#' )
#'
#' # Rebuild da_results from the SummarizedExperiment alone
#' da_out_from_se <- da_results_from_se(se)
#' str(da_out_from_se, max.level = 1)
#'
#' @seealso \code{\link{run_de}}, \code{\link{enrichmet}}
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom methods is
#' @export
da_results_from_se <- function(se,
                               fc_threshold      = 1,
                               pval_threshold    = 0.05,
                               split_complex_ids = TRUE) {
    
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object")
    }
    
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    
    required <- c("met_id", "log2fc", "pval", "padj")
    missing  <- setdiff(required, colnames(rd))
    if (length(missing) > 0) {
        stop("rowData(se) is missing required columns: ",
             paste(missing, collapse = ", "))
    }
    
    # Recompute significance labels if not already present
    if (!"Significant" %in% colnames(rd)) {
        rd$Significant <- "Not significant"
        rd$Significant[rd$padj < pval_threshold &
                           rd$log2fc >  fc_threshold] <- "Up"
        rd$Significant[rd$padj < pval_threshold &
                           rd$log2fc < -fc_threshold] <- "Down"
    }
    
    # Rebuild kegg_ready
    if (split_complex_ids) {
        kegg_list <- list()
        for (i in seq_len(nrow(rd))) {
            met_id     <- rd$met_id[i]
            candidates <- character(0)
            
            if (grepl("_", met_id)) {
                last_part  <- sub(".*_", "", met_id)
                candidates <- unlist(strsplit(last_part, "\\|"))
            } else if (grepl("^C\\d+", met_id)) {
                candidates <- unlist(strsplit(met_id, "\\|"))
            }
            
            candidates <- trimws(candidates)
            candidates <- unique(candidates[
                candidates != "" & grepl("^C\\d+", candidates)
            ])
            
            for (kegg_id in candidates) {
                new_row         <- rd[i, ]
                new_row$kegg_id <- kegg_id
                kegg_list[[length(kegg_list) + 1]] <- new_row
            }
        }
        
        if (length(kegg_list) > 0) {
            kegg_ready <- do.call(rbind, kegg_list)
            rownames(kegg_ready) <- NULL
        } else {
            kegg_ready <- rd
            kegg_ready$kegg_id <- NA_character_
            kegg_ready <- kegg_ready[!is.na(kegg_ready$kegg_id), ]
        }
    } else {
        kegg_ready <- rd
        kegg_ready$kegg_id <- vapply(
            strsplit(rd$met_id, "_"),
            function(x) {
                hit <- grep("^C\\d+", x, value = TRUE)
                if (length(hit) > 0) trimws(hit[1]) else NA_character_
            },
            FUN.VALUE = character(1)
        )
        kegg_ready <- kegg_ready[!is.na(kegg_ready$kegg_id), ]
        rownames(kegg_ready) <- NULL
    }
    
    # Trim whitespace from all kegg_id values
    kegg_ready$kegg_id <- trimws(kegg_ready$kegg_id)
    
    # Rebuild summary statistics
    summary_stats <- list(
        total_metabolites       = nrow(rd),
        significant_metabolites = sum(rd$Significant != "Not significant"),
        upregulated             = sum(rd$Significant == "Up"),
        downregulated           = sum(rd$Significant == "Down"),
        fc_threshold            = fc_threshold,
        pval_threshold          = pval_threshold,
        split_complex_ids       = split_complex_ids,
        unique_kegg_ids         = length(unique(kegg_ready$kegg_id))
    )
    
    list(
        full_results  = rd,
        kegg_ready    = kegg_ready,
        summary_stats = summary_stats
    )
}