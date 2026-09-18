
# ==============================================================================
# TCGA-LGG DNMT3A Main Figure 4
# Tumor Microenvironment
# FINAL V3 - Immune/TME concordance; xCell excluded
#
# Final figure panels:
#   A. Hedges' g heatmap: ESTIMATE + MCP-counter + Bindea cell-type ssGSEA
#   B. Positive-positive expression/methylation immune concordance
#   C. Method-wise concordance summary
#   xCell is excluded from the final analysis/figure
#
# Primary group comparisons:
#   1. DNMT3A expression: High vs Low
#   2. cg21629895 methylation: Low vs High
#
# IMPORTANT:
# - All categorical High/Low group assignments are read directly from the
#   completed DNMT3A master file.
# - No median group is recalculated in this script.
# - RNA expression is expected to be log2(norm_count + 1).
# - The RNA input should be the annotated legacy TCGA-LGG HiSeqV2 matrix with
#   Gene_Symbol plus TCGA sample-core columns.
# - Effect sizes are Hedges' g.
# - Positive Hedges' g means higher values in the adverse-direction group:
#       DNMT3A High vs Low
#       cg21629895 Low vs High
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Package setup
# ------------------------------------------------------------------------------

cran_packages <- c(
    "data.table",
    "dplyr",
    "tidyr",
    "tibble",
    "stringr",
    "ggplot2",
    "cowplot",
    "patchwork",
    "scales",
    "grid",
    "ggrepel",
    "msigdbr",
    "tidyestimate",
    "IOBR"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
    "GSVA",
    "MCPcounter"
)

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )
    }
}

required_packages <- c(
    cran_packages,
    "GSVA",
    "MCPcounter"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {
    stop(
        paste0(
            "The following required packages could not be loaded: ",
            paste(missing_packages, collapse = ", ")
        )
    )
}

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scales)
library(grid)


# ------------------------------------------------------------------------------
# 1. Publication settings
# ------------------------------------------------------------------------------

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 225
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"
FONT_PANEL_TAG <- 6.8
FONT_TITLE <- 6.4
FONT_AXIS_TITLE <- 5.6
FONT_AXIS_TEXT <- 5.2
FONT_LEGEND <- 5.0
FONT_STAT <- 4.9

INK <- "#202124"
MUTED <- "#6B6F76"

COMPARISON_COLORS <- c(
    "DNMT3A High vs Low" = "#6B5B95",
    "cg21629895 Low vs High" = "#3F7F6D"
)

INTEGRATED_COMPARISON_COLORS <- c(
    "DNMT3A High vs Low" = "#6B5B95",
    "cg21629895 Low vs High" = "#3F7F6D",
    "11-CpG promoter Low vs High" = "#C47A3A"
)

HEATMAP_LOW <- "#3B6FB6"
HEATMAP_MID <- "#F7F7F7"
HEATMAP_HIGH <- "#C94C4C"

pt_to_mm <- function(pt) {
    pt * 0.3527778
}

safe_numeric <- function(x) {
    suppressWarnings(
        as.numeric(
            as.character(x)
        )
    )
}


# ------------------------------------------------------------------------------
# 2. Select input files
# ------------------------------------------------------------------------------

cat(
    "\n1. Select the completed DNMT3A master CSV file.\n",
    sep = ""
)

master_file <- file.choose()

cat(
    "\n2. Select the annotated TCGA-LGG HiSeqV2 RNA expression CSV file.\n",
    sep = ""
)

rna_file <- file.choose()


file_1 <- data.table::fread(
    master_file,
    check.names = FALSE
)

file_2 <- data.table::fread(
    rna_file,
    check.names = FALSE
)

# Identify the two selected files by their actual column structure.
# This prevents an accidental master/RNA file swap from generating a
# misleading "all required columns are missing" error.

master_signature <- c(
    "Patient_ID",
    "Sample_Core",
    "DNMT3A_expression",
    "DNMT3A_expression_Median_Group",
    "cg21629895",
    "cg21629895_Median_Group"
)

is_master_like <- function(x) {
    all(
        master_signature %in%
            names(x)
    )
}

is_rna_like <- function(x) {
    "Gene_Symbol" %in%
        names(x)
}

if (
    is_master_like(file_1) &&
    is_rna_like(file_2)
) {

    master <- file_1
    rna_raw <- file_2

} else if (
    is_rna_like(file_1) &&
    is_master_like(file_2)
) {

    cat(
        "\nThe two selected files were reversed. They were automatically reassigned correctly.\n",
        sep = ""
    )

    master <- file_2
    rna_raw <- file_1

} else {

    cat("\nSelected file 1 first 20 columns:\n")
    print(
        head(
            names(file_1),
            20
        )
    )

    cat("\nSelected file 2 first 20 columns:\n")
    print(
        head(
            names(file_2),
            20
        )
    )

    stop(
        paste0(
            "The selected files could not be identified as the required pair.\n",
            "Please select:\n",
            "1. The completed DNMT3A master CSV containing Patient_ID, Sample_Core, ",
            "DNMT3A_expression, DNMT3A_expression_Median_Group, cg21629895, ",
            "and cg21629895_Median_Group\n",
            "2. The annotated TCGA-LGG HiSeqV2 RNA CSV containing Gene_Symbol."
        )
    )
}

cat("\nMaster file identified successfully:\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")

cat("\nRNA matrix identified successfully:\n")
cat("Rows:", nrow(rna_raw), "\n")
cat("Columns:", ncol(rna_raw), "\n")


# ------------------------------------------------------------------------------
# 3. Select output directory
# ------------------------------------------------------------------------------

choose_output_dir <- function() {

    selected <- NULL

    if (.Platform$OS.type == "windows") {

        selected <- tryCatch(
            utils::choose.dir(
                default = Sys.getenv("USERPROFILE"),
                caption = "Select a SHORT output folder for DNMT3A Figure 4"
            ),
            error = function(e) NULL
        )
    }

    if (
        is.null(selected) ||
        length(selected) == 0 ||
        is.na(selected) ||
        !nzchar(selected)
    ) {

        desktop <- file.path(
            Sys.getenv("USERPROFILE"),
            "Desktop"
        )

        if (!dir.exists(desktop)) {
            desktop <- getwd()
        }

        selected <- file.path(
            desktop,
            "DNMT3A_Fig4_TME"
        )
    }

    if (!dir.exists(selected)) {
        dir.create(
            selected,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    normalizePath(
        selected,
        winslash = "/",
        mustWork = TRUE
    )
}

OUTPUT_DIR <- choose_output_dir()

cat("\nAll Figure 4 outputs will be saved to:\n")
cat(OUTPUT_DIR, "\n")


# ------------------------------------------------------------------------------
# 4. Validate master-file columns
# ------------------------------------------------------------------------------

required_master_columns <- c(
    "Patient_ID",
    "Sample_Core",
    "DNMT3A_expression",
    "DNMT3A_expression_Median_Group",
    "cg21629895",
    "cg21629895_Median_Group"
)

missing_master_columns <- setdiff(
    required_master_columns,
    names(master)
)

if (length(missing_master_columns) > 0) {
    stop(
        paste0(
            "Required master-file columns are missing:\n",
            paste(
                missing_master_columns,
                collapse = ", "
            )
        )
    )
}

if (!"Gene_Symbol" %in% names(rna_raw)) {
    stop(
        "The RNA file must contain a Gene_Symbol column."
    )
}


# ------------------------------------------------------------------------------
# 5. Prepare master-file grouping variables
# ------------------------------------------------------------------------------

# These High/Low assignments are read directly from the completed master file.
# No median cutoff is recalculated in this script.

master[
    ,
    DNMT3A_Group := factor(
        trimws(
            as.character(
                DNMT3A_expression_Median_Group
            )
        ),
        levels = c(
            "Low",
            "High"
        )
    )
]

master[
    ,
    cg21629895_Group := factor(
        trimws(
            as.character(
                cg21629895_Median_Group
            )
        ),
        levels = c(
            "High",
            "Low"
        )
    )
]

master_group_qc <- data.table::data.table(
    Group = c(
        "DNMT3A mRNA expression",
        "cg21629895 methylation"
    ),
    Source_Column = c(
        "DNMT3A_expression_Median_Group",
        "cg21629895_Median_Group"
    ),
    Low_N = c(
        sum(
            master$DNMT3A_Group == "Low",
            na.rm = TRUE
        ),
        sum(
            master$cg21629895_Group == "Low",
            na.rm = TRUE
        )
    ),
    High_N = c(
        sum(
            master$DNMT3A_Group == "High",
            na.rm = TRUE
        ),
        sum(
            master$cg21629895_Group == "High",
            na.rm = TRUE
        )
    ),
    Missing_N = c(
        sum(
            is.na(
                master$DNMT3A_Group
            )
        ),
        sum(
            is.na(
                master$cg21629895_Group
            )
        )
    )
)

cat("\nMaster-file grouping QC:\n")
print(master_group_qc)

if (
    any(
        master_group_qc$Missing_N > 0
    )
) {
    warning(
        "Missing median-group assignments were detected in the master file."
    )
}


# ------------------------------------------------------------------------------
# 6. Prepare RNA expression matrix
# ------------------------------------------------------------------------------

sample_columns <- names(rna_raw)[
    grepl(
        "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
        names(rna_raw)
    )
]

if (length(sample_columns) == 0) {
    stop(
        paste0(
            "No TCGA sample-core columns were detected in the RNA file. ",
            "Expected names similar to TCGA-XX-XXXX-01."
        )
    )
}

# Keep primary tumor samples only.
sample_columns <- sample_columns[
    substr(
        sample_columns,
        14,
        15
    ) == "01"
]

overlap_samples <- intersect(
    master$Sample_Core,
    sample_columns
)

if (length(overlap_samples) < 2) {
    stop(
        "Too few overlapping primary tumor samples between master and RNA files."
    )
}

cat("\nRNA/master sample overlap:\n")
cat("Master samples:", length(unique(master$Sample_Core)), "\n")
cat("RNA primary tumor samples:", length(sample_columns), "\n")
cat("Overlapping samples:", length(overlap_samples), "\n")


rna_dt <- rna_raw[
    !is.na(Gene_Symbol) &
        Gene_Symbol != "",
    c(
        "Gene_Symbol",
        overlap_samples
    ),
    with = FALSE
]

# Collapse duplicate gene symbols by mean if any are present.
rna_long <- data.table::melt(
    rna_dt,
    id.vars = "Gene_Symbol",
    variable.name = "Sample_Core",
    value.name = "Expression"
)

rna_long[
    ,
    Expression := safe_numeric(Expression)
]

rna_collapsed <- rna_long[
    is.finite(Expression),
    .(
        Expression = mean(
            Expression,
            na.rm = TRUE
        )
    ),
    by = .(
        Gene_Symbol,
        Sample_Core
    )
]

rna_wide <- data.table::dcast(
    rna_collapsed,
    Gene_Symbol ~ Sample_Core,
    value.var = "Expression"
)

gene_symbols <- rna_wide$Gene_Symbol

expr_matrix <- as.matrix(
    rna_wide[
        ,
        -1,
        with = FALSE
    ]
)

rownames(expr_matrix) <- gene_symbols

storage.mode(expr_matrix) <- "numeric"

# Reorder expression columns to the master-file sample order.
sample_order <- master$Sample_Core[
    master$Sample_Core %in% colnames(expr_matrix)
]

sample_order <- unique(sample_order)

expr_matrix <- expr_matrix[
    ,
    sample_order,
    drop = FALSE
]

analysis_master <- master[
    match(
        sample_order,
        Sample_Core
    )
]

if (!all(
    analysis_master$Sample_Core == colnames(expr_matrix)
)) {
    stop(
        "Master/RNA sample order mismatch after alignment."
    )
}

cat("\nAligned expression matrix:\n")
cat("Genes:", nrow(expr_matrix), "\n")
cat("Samples:", ncol(expr_matrix), "\n")


# ------------------------------------------------------------------------------
# 7. General effect-size/statistics helper
# ------------------------------------------------------------------------------

hedges_g <- function(
        x_group1,
        x_group0
) {

    x1 <- x_group1[
        is.finite(x_group1)
    ]

    x0 <- x_group0[
        is.finite(x_group0)
    ]

    n1 <- length(x1)
    n0 <- length(x0)

    if (
        n1 < 3 ||
        n0 < 3
    ) {
        return(
            c(
                g = NA_real_,
                lower = NA_real_,
                upper = NA_real_
            )
        )
    }

    s1 <- stats::sd(x1)
    s0 <- stats::sd(x0)

    pooled_sd <- sqrt(
        (
            (n1 - 1) * s1^2 +
                (n0 - 1) * s0^2
        ) /
            (
                n1 + n0 - 2
            )
    )

    if (
        !is.finite(pooled_sd) ||
        pooled_sd == 0
    ) {
        return(
            c(
                g = 0,
                lower = 0,
                upper = 0
            )
        )
    }

    d <- (
        mean(
            x1,
            na.rm = TRUE
        ) -
            mean(
                x0,
                na.rm = TRUE
            )
    ) / pooled_sd

    dfree <- n1 + n0 - 2

    correction <- 1 - 3 / (
        4 * dfree - 1
    )

    g <- correction * d

    se_g <- sqrt(
        (
            n1 + n0
        ) /
            (
                n1 * n0
            ) +
            (
                g^2
            ) /
                (
                    2 * dfree
                )
    )

    c(
        g = g,
        lower = g - 1.96 * se_g,
        upper = g + 1.96 * se_g
    )
}


compare_score_matrix <- function(
        score_matrix,
        feature_type
) {

    score_matrix <- score_matrix[
        ,
        analysis_master$Sample_Core,
        drop = FALSE
    ]

    comparison_specs <- list(
        list(
            name = "DNMT3A High vs Low",
            group = analysis_master$DNMT3A_Group,
            group1 = "High",
            group0 = "Low"
        ),
        list(
            name = "cg21629895 Low vs High",
            group = analysis_master$cg21629895_Group,
            group1 = "Low",
            group0 = "High"
        )
    )

    output <- list()

    for (spec in comparison_specs) {

        for (feature_name in rownames(score_matrix)) {

            values <- safe_numeric(
                score_matrix[
                    feature_name,
                    ,
                    drop = TRUE
                ]
            )

            group_vector <- spec$group

            x1 <- values[
                group_vector == spec$group1
            ]

            x0 <- values[
                group_vector == spec$group0
            ]

            effect <- hedges_g(
                x1,
                x0
            )

            p_value <- tryCatch(
                stats::wilcox.test(
                    x1,
                    x0,
                    exact = FALSE
                )$p.value,
                error = function(e) NA_real_
            )

            output[[
                length(output) + 1
            ]] <- data.table::data.table(
                Feature_Type = feature_type,
                Feature = feature_name,
                Comparison = spec$name,
                Group1 = spec$group1,
                Group0 = spec$group0,
                N_Group1 = sum(
                    is.finite(x1)
                ),
                N_Group0 = sum(
                    is.finite(x0)
                ),
                Median_Group1 = stats::median(
                    x1,
                    na.rm = TRUE
                ),
                Median_Group0 = stats::median(
                    x0,
                    na.rm = TRUE
                ),
                Hedges_g = effect["g"],
                CI_Lower = effect["lower"],
                CI_Upper = effect["upper"],
                P_Value = p_value
            )
        }
    }

    result <- data.table::rbindlist(
        output,
        fill = TRUE
    )

    result[
        ,
        FDR := stats::p.adjust(
            P_Value,
            method = "BH"
        ),
        by = Comparison
    ]

    result
}


format_q <- function(q) {

    if (is.na(q)) {
        return("")
    }

    if (q < 0.001) {
        return("q<0.001")
    }

    paste0(
        "q=",
        sprintf(
            "%.3f",
            q
        )
    )
}


# ------------------------------------------------------------------------------
# 8. ESTIMATE
# ------------------------------------------------------------------------------

# The legacy R-Forge estimate package can fail in newer R installations because
# package-data objects such as common_genes or SI_geneset are not available at
# runtime. This script therefore uses tidyestimate, a modern implementation of
# the ESTIMATE algorithm that contains its required common-gene set and stromal/
# immune signatures as package data.
#
# RNA-seq data are analyzed with is_affymetrix = FALSE. Therefore, only the
# relative StromalScore, ImmuneScore, and ESTIMATEScore are used; no absolute
# tumor-purity estimate is calculated.

estimate_input <- data.frame(
    Gene_Symbol = rownames(
        expr_matrix
    ),
    expr_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Restrict the matrix to the ESTIMATE cross-platform common-gene universe.
# Gene identifiers are HGNC symbols. Alias reassignment is disabled to avoid
# heuristic gene-name substitutions.
estimate_filtered <- tidyestimate::filter_common_genes(
    estimate_input,
    id = "hgnc_symbol",
    tidy = TRUE,
    tell_missing = FALSE,
    find_alias = FALSE
)

if (
    is.null(
        estimate_filtered
    ) ||
    nrow(
        estimate_filtered
    ) < 100
) {
    stop(
        paste0(
            "tidyestimate common-gene filtering returned too few genes: ",
            ifelse(
                is.null(estimate_filtered),
                0,
                nrow(estimate_filtered)
            ),
            "."
        )
    )
}

cat("\nESTIMATE common-gene filtering completed.\n")
cat("Genes before filtering:", nrow(expr_matrix), "\n")
cat("Genes after filtering:", nrow(estimate_filtered), "\n")
cat("Samples:", ncol(expr_matrix), "\n")

estimate_scores_df <- tidyestimate::estimate_score(
    estimate_filtered,
    is_affymetrix = FALSE
)

required_estimate_columns <- c(
    "sample",
    "stromal",
    "immune",
    "estimate"
)

missing_estimate_columns <- setdiff(
    required_estimate_columns,
    names(
        estimate_scores_df
    )
)

if (
    length(
        missing_estimate_columns
    ) > 0
) {
    stop(
        paste0(
            "tidyestimate output is missing expected columns: ",
            paste(
                missing_estimate_columns,
                collapse = ", "
            )
        )
    )
}

estimate_scores_dt <- data.table::as.data.table(
    estimate_scores_df
)

estimate_scores_dt[
    ,
    sample := as.character(
        sample
    )
]

missing_estimate_samples <- setdiff(
    analysis_master$Sample_Core,
    estimate_scores_dt$sample
)

if (
    length(
        missing_estimate_samples
    ) > 0
) {
    stop(
        paste0(
            "tidyestimate output is missing ",
            length(
                missing_estimate_samples
            ),
            " aligned TCGA samples."
        )
    )
}

estimate_scores_dt <- estimate_scores_dt[
    match(
        analysis_master$Sample_Core,
        sample
    )
]

if (
    !all(
        estimate_scores_dt$sample ==
            analysis_master$Sample_Core
    )
) {
    stop(
        "ESTIMATE sample order does not match the aligned master file."
    )
}

estimate_matrix <- rbind(
    ImmuneScore = safe_numeric(
        estimate_scores_dt$immune
    ),
    StromalScore = safe_numeric(
        estimate_scores_dt$stromal
    ),
    ESTIMATEScore = safe_numeric(
        estimate_scores_dt$estimate
    )
)

colnames(
    estimate_matrix
) <- estimate_scores_dt$sample

if (
    any(
        !is.finite(
            estimate_matrix
        )
    )
) {
    stop(
        "Non-finite values were detected in the ESTIMATE score matrix."
    )
}

cat("\nESTIMATE scoring completed successfully with tidyestimate.\n")
cat("Score matrix:", nrow(estimate_matrix), "scores x",
    ncol(estimate_matrix), "samples\n")

estimate_stats <- compare_score_matrix(
    estimate_matrix,
    "ESTIMATE"
)


# ------------------------------------------------------------------------------
# 9. MCP-counter
# ------------------------------------------------------------------------------

mcp_matrix <- MCPcounter::MCPcounter.estimate(
    expr_matrix,
    featuresType = "HUGO_symbols"
)

mcp_matrix <- as.matrix(
    mcp_matrix
)

mcp_matrix <- mcp_matrix[
    ,
    analysis_master$Sample_Core,
    drop = FALSE
]

mcp_stats <- compare_score_matrix(
    mcp_matrix,
    "MCP-counter"
)


# ------------------------------------------------------------------------------
# 11. Cell-type ssGSEA using published Bindea immune-cell signatures
# ------------------------------------------------------------------------------

# Cell-type ssGSEA is based on the published immune-cell signatures from:
# Bindea et al., Immunity 2013; PMID: 24138885.
#
# IOBR distributes these curated signatures in signature_collection.
# We select broad immune-cell populations so that Panel D can use the same
# forest-plot presentation as MCP-counter and xCell.

signature_env <- new.env(
    parent = emptyenv()
)

suppressWarnings(
    utils::data(
        list = "signature_collection",
        package = "IOBR",
        envir = signature_env
    )
)

if (
    !exists(
        "signature_collection",
        envir = signature_env,
        inherits = FALSE
    )
) {
    stop(
        paste0(
            "IOBR::signature_collection could not be loaded. ",
            "Please reinstall/update IOBR and rerun the script."
        )
    )
}

signature_collection_iobr <- get(
    "signature_collection",
    envir = signature_env,
    inherits = FALSE
)

bindea_all <- signature_collection_iobr[
    grep(
        "Bindea_et_al",
        names(
            signature_collection_iobr
        ),
        ignore.case = TRUE
    )
]

if (length(bindea_all) < 10) {
    stop(
        paste0(
            "Fewer than 10 Bindea immune-cell signatures were found in IOBR. ",
            "Detected signatures: ",
            paste(
                names(bindea_all),
                collapse = ", "
            )
        )
    )
}

cat("\nBindea immune-cell signatures detected in IOBR:\n")
print(
    names(
        bindea_all
    )
)


# Match one published Bindea signature to each broad cell population.
# Matching is done against the installed IOBR signature names and therefore
# does not hard-code unpublished or study-specific marker genes.

pick_bindea_signature <- function(
        display_label,
        include_patterns,
        exclude_patterns = NULL
) {

    available_names <- names(
        bindea_all
    )

    keep <- rep(
        TRUE,
        length(
            available_names
        )
    )

    if (length(include_patterns) > 0) {

        include_hit <- rep(
            FALSE,
            length(
                available_names
            )
        )

        for (pat in include_patterns) {
            include_hit <- include_hit |
                grepl(
                    pat,
                    available_names,
                    ignore.case = TRUE
                )
        }

        keep <- keep &
            include_hit
    }

    if (!is.null(exclude_patterns)) {

        exclude_hit <- rep(
            FALSE,
            length(
                available_names
            )
        )

        for (pat in exclude_patterns) {
            exclude_hit <- exclude_hit |
                grepl(
                    pat,
                    available_names,
                    ignore.case = TRUE
                )
        }

        keep <- keep &
            !exclude_hit
    }

    hits <- available_names[
        keep
    ]

    if (length(hits) == 0) {
        return(
            NULL
        )
    }

    # Prefer the shortest matching signature name when more than one exists.
    hits <- hits[
        order(
            nchar(
                hits
            )
        )
    ]

    selected_name <- hits[1]

    out <- bindea_all[
        selected_name
    ]

    names(out) <- display_label

    out
}


bindea_target_specs <- list(
    "T cells" = list(
        include = c(
            "(^|_)T[_ ]?cells",
            "T_cells"
        ),
        exclude = c(
            "CD8",
            "CD4",
            "Th1",
            "Th2",
            "Treg",
            "regulatory",
            "gamma",
            "cytotoxic"
        )
    ),
    "CD8 T cells" = list(
        include = c(
            "CD8"
        ),
        exclude = NULL
    ),
    "B cells" = list(
        include = c(
            "(^|_)B[_ ]?cells",
            "B_cells"
        ),
        exclude = c(
            "memory",
            "plasma"
        )
    ),
    "NK cells" = list(
        include = c(
            "NK",
            "natural_killer"
        ),
        exclude = NULL
    ),
    "Cytotoxic lymphocytes" = list(
        include = c(
            "cytotoxic"
        ),
        exclude = NULL
    ),
    "Macrophages" = list(
        include = c(
            "macroph"
        ),
        exclude = NULL
    ),
    "Dendritic cells" = list(
        include = c(
            "dendritic",
            "(^|_)DC"
        ),
        exclude = NULL
    ),
    "Neutrophils" = list(
        include = c(
            "neutroph"
        ),
        exclude = NULL
    ),
    "Th1 cells" = list(
        include = c(
            "Th1"
        ),
        exclude = NULL
    ),
    "Treg cells" = list(
        include = c(
            "Treg",
            "regulatory"
        ),
        exclude = NULL
    )
)


bindea_selected <- list()

for (display_label in names(
    bindea_target_specs
)) {

    spec <- bindea_target_specs[[
        display_label
    ]]

    selected <- pick_bindea_signature(
        display_label = display_label,
        include_patterns = spec$include,
        exclude_patterns = spec$exclude
    )

    if (!is.null(selected)) {
        bindea_selected <- c(
            bindea_selected,
            selected
        )
    }
}


# Require enough broad cell populations for a publication-quality cell-type panel.
if (length(bindea_selected) < 6) {
    stop(
        paste0(
            "Only ",
            length(bindea_selected),
            " broad Bindea cell-type signatures could be matched. ",
            "Available Bindea signatures in the installed IOBR version:\n",
            paste(
                names(bindea_all),
                collapse = "\n"
            )
        )
    )
}


# Remove duplicate gene sets if two labels unexpectedly resolve to the same genes.
signature_key <- vapply(
    bindea_selected,
    function(x) {
        paste(
            sort(
                unique(
                    as.character(
                        x
                    )
                )
            ),
            collapse = "|"
        )
    },
    character(1)
)

bindea_selected <- bindea_selected[
    !duplicated(
        signature_key
    )
]


cat("\nSelected Bindea cell-type signatures for ssGSEA:\n")
print(
    names(
        bindea_selected
    )
)


run_ssgsea <- function(
        expression_matrix,
        gene_sets
) {

    if (
        exists(
            "ssgseaParam",
            where = asNamespace(
                "GSVA"
            ),
            inherits = FALSE
        )
    ) {

        param <- GSVA::ssgseaParam(
            expression_matrix,
            gene_sets,
            normalize = TRUE,
            minSize = 5
        )

        return(
            GSVA::gsva(
                param,
                verbose = FALSE
            )
        )
    }

    GSVA::gsva(
        expression_matrix,
        gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = TRUE,
        min.sz = 5,
        max.sz = 500,
        verbose = FALSE
    )
}


ssgsea_matrix <- run_ssgsea(
    expr_matrix,
    bindea_selected
)

ssgsea_matrix <- as.matrix(
    ssgsea_matrix
)

storage.mode(
    ssgsea_matrix
) <- "numeric"

missing_ssgsea_samples <- setdiff(
    analysis_master$Sample_Core,
    colnames(
        ssgsea_matrix
    )
)

if (length(
    missing_ssgsea_samples
) > 0) {
    stop(
        paste0(
            "Cell-type ssGSEA output is missing ",
            length(
                missing_ssgsea_samples
            ),
            " aligned TCGA samples."
        )
    )
}

ssgsea_matrix <- ssgsea_matrix[
    ,
    analysis_master$Sample_Core,
    drop = FALSE
]

ssgsea_stats <- compare_score_matrix(
    ssgsea_matrix,
    "ssGSEA"
)

cat("\nCell-type ssGSEA completed successfully.\n")
cat("Cell signatures:", nrow(ssgsea_matrix), "\n")
cat("Samples:", ncol(ssgsea_matrix), "\n")


# ------------------------------------------------------------------------------
# 12. Save raw score matrices
# ------------------------------------------------------------------------------

save_score_matrix <- function(
        score_matrix,
        filename
) {

    out <- data.table::as.data.table(
        score_matrix,
        keep.rownames = "Feature"
    )

    data.table::fwrite(
        out,
        file.path(
            OUTPUT_DIR,
            filename
        ),
        na = ""
    )
}

save_score_matrix(
    estimate_matrix,
    "Fig4_ESTIMATE_Scores.csv"
)

save_score_matrix(
    mcp_matrix,
    "Fig4_MCPcounter_Scores.csv"
)



save_score_matrix(
    ssgsea_matrix,
    "Fig4_Immune_ssGSEA_Scores.csv"
)


# ------------------------------------------------------------------------------
# 13. Combined statistics
# ------------------------------------------------------------------------------

all_stats <- data.table::rbindlist(
    list(
        estimate_stats,
        mcp_stats,
        ssgsea_stats
    ),
    fill = TRUE
)

data.table::fwrite(
    all_stats,
    file.path(
        OUTPUT_DIR,
        "Fig4_All_TME_EffectSize_Statistics.csv"
    ),
    na = ""
)

data.table::fwrite(
    master_group_qc,
    file.path(
        OUTPUT_DIR,
        "Fig4_Master_Group_QC.csv"
    ),
    na = ""
)





# ------------------------------------------------------------------------------
# 14. FINAL MAIN FIGURE 4
# Immune/TME concordance between DNMT3A expression and cg21629895 methylation
# FINAL V2 - manuscript layout
# ------------------------------------------------------------------------------

# Final layout:
#
#   A. Large Hedges' g heatmap on the left
#
#                           B. Positive-positive concordance scatter
#
#                           C. Method-wise concordance summary
#
# Included methods:
#   - ESTIMATE
#   - MCP-counter
#   - Cell-type ssGSEA using Bindea immune-cell signatures
#
# xCell is intentionally excluded from the final main figure.
#
# Comparison directions:
#   DNMT3A mRNA: High vs Low
#   cg21629895 methylation: Low vs High
#
# Positive Hedges' g therefore indicates:
#   - higher score in DNMT3A-high tumors
#   - higher score in cg21629895-low tumors
#
# Panel B displays the full four-quadrant effect-size landscape.
# The Spearman correlation is calculated using ALL immune/TME features.
#
# Panel C summarizes ALL features and therefore preserves concordant-negative,
# discordant-significant, and partial/not-significant categories.

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 158
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"
FONT_PANEL <- 7.0
FONT_TITLE <- 6.5
FONT_AXIS_TITLE <- 5.6
FONT_AXIS_TEXT <- 5.0
FONT_LEGEND <- 5.0
FONT_HEATMAP <- 5.0
FONT_CELL <- 5.0
FONT_STAT <- 5.0

INK <- "#222222"
MUTED <- "#666666"

METHOD_FILL <- c(
    "ESTIMATE" = "#DDEAF7",
    "MCP-counter" = "#DCEEDC",
    "Cell-type ssGSEA" = "#FCE4CF"
)

METHOD_TEXT <- c(
    "ESTIMATE" = "#1F5A99",
    "MCP-counter" = "#2E7D4F",
    "Cell-type ssGSEA" = "#C54E12"
)

METHOD_POINT_COLORS <- c(
    "ESTIMATE" = "#1F5A99",
    "MCP-counter" = "#2E7D4F",
    "Cell-type ssGSEA" = "#E6A400"
)

METHOD_SHAPES <- c(
    "ESTIMATE" = 15,
    "MCP-counter" = 16,
    "Cell-type ssGSEA" = 17
)

CONCORDANCE_COLORS <- c(
    "Concordant positive" = "#C9474D",
    "Concordant negative" = "#2C6DB2",
    "Discordant significant" = "#7A7A7A",
    "Partial / not significant" = "#D9D9D9"
)

pt_to_mm <- function(pt) {
    pt * 0.3527778
}


# ------------------------------------------------------------------------------
# 15. Prespecified immune/TME features
# ------------------------------------------------------------------------------

estimate_order <- c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
)

mcp_order <- c(
    "T cells",
    "CD8 T cells",
    "Cytotoxic lymphocytes",
    "B lineage",
    "NK cells",
    "Monocytic lineage",
    "Myeloid dendritic cells",
    "Neutrophils",
    "Endothelial cells",
    "Fibroblasts"
)

ssgsea_order <- names(
    bindea_selected
)


# ------------------------------------------------------------------------------
# 16. Build shared immune/TME effect-size table
# ------------------------------------------------------------------------------

all_stats_final <- data.table::rbindlist(
    list(
        estimate_stats,
        mcp_stats,
        ssgsea_stats
    ),
    fill = TRUE
)

all_stats_final <- all_stats_final[
    Comparison %in% c(
        "DNMT3A High vs Low",
        "cg21629895 Low vs High"
    )
]

all_stats_final[
    ,
    Method := data.table::fcase(
        Feature_Type == "ESTIMATE", "ESTIMATE",
        grepl("MCP", Feature_Type, ignore.case = TRUE), "MCP-counter",
        grepl("ssGSEA", Feature_Type, ignore.case = TRUE), "Cell-type ssGSEA",
        default = as.character(Feature_Type)
    )
]

all_stats_final[
    ,
    Comparison_Display := data.table::fcase(
        Comparison == "DNMT3A High vs Low",
        "DNMT3A mRNA\nHigh vs Low",
        Comparison == "cg21629895 Low vs High",
        "cg21629895 methylation\nLow vs High",
        default = NA_character_
    )
]

comparison_levels <- c(
    "DNMT3A mRNA\nHigh vs Low",
    "cg21629895 methylation\nLow vs High"
)

all_stats_final[
    ,
    Comparison_Display := factor(
        Comparison_Display,
        levels = comparison_levels
    )
]

all_stats_final[
    ,
    Sig_Label := data.table::fcase(
        is.finite(FDR) & FDR < 0.001, "***",
        is.finite(FDR) & FDR < 0.01, "**",
        is.finite(FDR) & FDR < 0.05, "*",
        default = "ns"
    )
]


# ------------------------------------------------------------------------------
# 17. Stable feature order
# ------------------------------------------------------------------------------

method_order <- c(
    "ESTIMATE",
    "MCP-counter",
    "Cell-type ssGSEA"
)

feature_order <- list(
    "ESTIMATE" = estimate_order,
    "MCP-counter" = mcp_order,
    "Cell-type ssGSEA" = ssgsea_order
)

row_map <- data.table::rbindlist(
    lapply(
        method_order,
        function(method_name) {

            available_features <- all_stats_final[
                Method == method_name,
                unique(Feature)
            ]

            present <- feature_order[[method_name]][
                feature_order[[method_name]] %in% available_features
            ]

            data.table::data.table(
                Method = method_name,
                Feature = present
            )
        }
    ),
    fill = TRUE
)

row_map[
    ,
    Row_ID := seq_len(.N)
]

all_stats_final <- merge(
    all_stats_final,
    row_map,
    by = c(
        "Method",
        "Feature"
    ),
    all.x = TRUE,
    sort = FALSE
)

all_stats_final <- all_stats_final[
    !is.na(Row_ID) &
        !is.na(Comparison_Display)
]


# ------------------------------------------------------------------------------
# 18. Panel A: Hedges' g heatmap
# ------------------------------------------------------------------------------

heat_limit <- max(
    abs(all_stats_final$Hedges_g),
    na.rm = TRUE
)

if (
    !is.finite(heat_limit) ||
    heat_limit <= 0
) {
    heat_limit <- 1
}

heat_limit <- ceiling(
    heat_limit * 10
) / 10


method_blocks <- row_map[
    ,
    .(
        ymin = min(Row_ID) - 0.5,
        ymax = max(Row_ID) + 0.5,
        ymid = mean(c(min(Row_ID), max(Row_ID)))
    ),
    by = Method
]

method_blocks[
    ,
    Method := factor(
        Method,
        levels = method_order
    )
]


method_strip <- ggplot(
    method_blocks,
    aes(
        ymin = ymin,
        ymax = ymax,
        fill = Method
    )
) +
    geom_rect(
        aes(
            xmin = 0,
            xmax = 1
        ),
        color = "white",
        linewidth = 0.75
    ) +
    geom_text(
        aes(
            x = 0.5,
            y = ymid,
            label = Method,
            color = Method
        ),
        size = pt_to_mm(5.0),
        fontface = "bold",
        lineheight = 0.90
    ) +
    scale_fill_manual(
        values = METHOD_FILL,
        guide = "none"
    ) +
    scale_color_manual(
        values = METHOD_TEXT,
        guide = "none"
    ) +
    scale_y_reverse(
        limits = c(
            max(row_map$Row_ID) + 0.5,
            0.5
        ),
        expand = c(0, 0)
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        expand = c(0, 0)
    ) +
    coord_cartesian(
        clip = "off"
    ) +
    theme_void() +
    theme(
        plot.margin = margin(
            0,
            -2,
            0,
            0
        )
    )


heatmap_core_base <- ggplot(
    all_stats_final,
    aes(
        x = Comparison_Display,
        y = Row_ID,
        fill = Hedges_g
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.72
    ) +
    geom_text(
        aes(
            label = Sig_Label
        ),
        size = pt_to_mm(FONT_CELL),
        fontface = "bold",
        color = INK
    ) +
    scale_fill_gradient2(
        low = "#2C6DB2",
        mid = "#F7F7F7",
        high = "#C9474D",
        midpoint = 0,
        limits = c(
            -heat_limit,
            heat_limit
        ),
        oob = scales::squish,
        name = "Hedges' g"
    ) +
    scale_y_reverse(
        limits = c(
            max(row_map$Row_ID) + 0.5,
            0.5
        ),
        breaks = row_map$Row_ID,
        labels = row_map$Feature,
        expand = c(0, 0)
    ) +
    scale_x_discrete(
        drop = FALSE
    ) +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_classic(
        base_size = 6,
        base_family = FONT_FAMILY
    ) +
    theme(
        axis.text.x = element_text(
            size = FONT_AXIS_TEXT,
            face = "bold",
            color = INK,
            lineheight = 0.90,
            margin = margin(t = 3)
        ),
        axis.text.y = element_text(
            size = FONT_HEATMAP,
            color = INK,
            hjust = 1,
            margin = margin(r = 1.5),
            lineheight = 0.90
        ),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = margin(
            0,
            2,
            0,
            -1
        )
    )


# Compact horizontal Hedges' g scale bar for Panel A.
heatmap_scale_plot <- ggplot(
    data.frame(
        x = c(
            -heat_limit,
            heat_limit
        ),
        y = 1
    ),
    aes(
        x = x,
        y = y,
        color = x
    )
) +
    geom_point(
        alpha = 0
    ) +
    scale_color_gradient2(
        low = "#2C6DB2",
        mid = "#F7F7F7",
        high = "#C9474D",
        midpoint = 0,
        limits = c(
            -heat_limit,
            heat_limit
        ),
        name = "Hedges' g"
    ) +
    guides(
        color = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(28, "mm"),
            barheight = unit(2.2, "mm"),
            ticks = TRUE,
            frame.colour = "#A0A0A0",
            frame.linewidth = 0.25
        )
    ) +
    theme_void(
        base_family = FONT_FAMILY
    ) +
    theme(
        legend.position = "bottom",
        legend.title = element_text(
            size = FONT_LEGEND,
            face = "bold"
        ),
        legend.text = element_text(
            size = FONT_LEGEND
        ),
        legend.margin = margin(
            0,
            0,
            0,
            0
        )
    )

heatmap_scale <- cowplot::get_legend(
    heatmap_scale_plot
)


panel_A_body <- cowplot::plot_grid(
    method_strip,
    heatmap_core_base,
    nrow = 1,
    rel_widths = c(
        0.20,
        0.80
    ),
    align = "h",
    axis = "tb"
)

panel_A_title <- cowplot::ggdraw() +
    cowplot::draw_label(
        "A",
        x = 0.002,
        y = 0.98,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    ) +
    cowplot::draw_label(
        "Immune/TME effect-size profile",
        x = 0.055,
        y = 0.98,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_TITLE,
        color = INK
    )

panel_A_footer <- cowplot::ggdraw() +
    cowplot::draw_label(
        "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001; ns = not significant",
        x = 0.00,
        y = 0.75,
        hjust = 0,
        vjust = 1,
        size = 5.0,
        color = MUTED,
        fontfamily = FONT_FAMILY
    )

panel_A <- cowplot::plot_grid(
    panel_A_title,
    panel_A_body,
    heatmap_scale,
    panel_A_footer,
    ncol = 1,
    rel_heights = c(
        0.07,
        0.80,
        0.07,
        0.06
    )
)


# ------------------------------------------------------------------------------
# 19. One-row-per-feature concordance table
# ------------------------------------------------------------------------------

expr_stats <- all_stats_final[
    Comparison == "DNMT3A High vs Low",
    .(
        Method,
        Feature,
        Hedges_g_mRNA = Hedges_g,
        FDR_mRNA = FDR
    )
]

meth_stats <- all_stats_final[
    Comparison == "cg21629895 Low vs High",
    .(
        Method,
        Feature,
        Hedges_g_Methylation = Hedges_g,
        FDR_Methylation = FDR
    )
]

concordance_dt <- merge(
    expr_stats,
    meth_stats,
    by = c(
        "Method",
        "Feature"
    ),
    all = FALSE,
    sort = FALSE
)

concordance_dt[
    ,
    Both_Significant :=
        is.finite(FDR_mRNA) &
        is.finite(FDR_Methylation) &
        FDR_mRNA < 0.05 &
        FDR_Methylation < 0.05
]

concordance_dt[
    ,
    Concordance_Class := data.table::fcase(
        Both_Significant &
            Hedges_g_mRNA > 0 &
            Hedges_g_Methylation > 0,
        "Concordant positive",

        Both_Significant &
            Hedges_g_mRNA < 0 &
            Hedges_g_Methylation < 0,
        "Concordant negative",

        Both_Significant &
            sign(Hedges_g_mRNA) !=
            sign(Hedges_g_Methylation),
        "Discordant significant",

        default = "Partial / not significant"
    )
]

concordance_levels <- c(
    "Concordant positive",
    "Concordant negative",
    "Discordant significant",
    "Partial / not significant"
)

concordance_dt[
    ,
    Concordance_Class := factor(
        Concordance_Class,
        levels = concordance_levels
    )
]

concordance_dt[
    ,
    Shared_Effect_Strength := pmin(
        abs(Hedges_g_mRNA),
        abs(Hedges_g_Methylation)
    )
]




# ------------------------------------------------------------------------------
# 20. Panel B: optimized full-quadrant effect-size concordance
# ------------------------------------------------------------------------------

# The complete four-quadrant landscape is retained.
# Every feature significant in BOTH comparisons (BH FDR < 0.05 in each) is labeled.
# Labels are shortened only for display and duplicate feature names are given
# method suffixes so that the plotted labels remain unambiguous.

cor_test_all <- suppressWarnings(
    stats::cor.test(
        concordance_dt$Hedges_g_mRNA,
        concordance_dt$Hedges_g_Methylation,
        method = "spearman",
        exact = FALSE
    )
)

rho_value_all <- unname(
    cor_test_all$estimate
)

rho_p_all <- cor_test_all$p.value

rho_text <- paste0(
    "All features: Spearman \u03c1 = ",
    sprintf("%.2f", rho_value_all),
    "; P ",
    ifelse(
        rho_p_all < 0.001,
        "< 0.001",
        paste0(
            "= ",
            sprintf("%.3f", rho_p_all)
        )
    ),
    "; n = ",
    nrow(concordance_dt)
)

# Label ALL features significant in both comparisons.
label_dt <- concordance_dt[
    Both_Significant == TRUE &
        is.finite(Hedges_g_mRNA) &
        is.finite(Hedges_g_Methylation)
]

# Create concise, unambiguous display labels.
label_dt[
    ,
    Label_Display := Feature
]

label_dt[
    Feature == "ESTIMATEScore",
    Label_Display := "ESTIMATE"
]

label_dt[
    Feature == "ImmuneScore",
    Label_Display := "Immune"
]

label_dt[
    Feature == "StromalScore",
    Label_Display := "Stromal"
]

label_dt[
    Feature == "Cytotoxic lymphocytes",
    Label_Display := "Cytotoxic\nlymphocytes"
]

label_dt[
    Feature == "Myeloid dendritic cells",
    Label_Display := "Myeloid DC"
]

label_dt[
    Feature == "Dendritic cells",
    Label_Display := "DC"
]

# Add method suffix only when the same feature label appears more than once.
duplicate_labels <- label_dt[
    ,
    .N,
    by = Label_Display
][
    N > 1,
    Label_Display
]

label_dt[
    Label_Display %in% duplicate_labels,
    Label_Display := paste0(
        Label_Display,
        ifelse(
            Method == "MCP-counter",
            " (MCP)",
            ifelse(
                Method == "Cell-type ssGSEA",
                " (ssGSEA)",
                " (EST)"
            )
        )
    )
]

# Use one symmetric axis range for clean quadrant comparison.
scatter_limit <- max(
    abs(
        c(
            concordance_dt$Hedges_g_mRNA,
            concordance_dt$Hedges_g_Methylation
        )
    ),
    na.rm = TRUE
)

if (
    !is.finite(scatter_limit) ||
    scatter_limit <= 0
) {
    scatter_limit <- 1
}

scatter_limit <- ceiling(
    scatter_limit * 10
) / 10

# Give labels a small amount of extra plotting room without changing the data scale.
scatter_plot_limit <- scatter_limit * 1.08

# Significant labels farther from the origin get slightly stronger priority.
label_dt[
    ,
    Label_Priority := sqrt(
        Hedges_g_mRNA^2 +
            Hedges_g_Methylation^2
    )
]

setorder(
    label_dt,
    -Label_Priority
)

panel_B_base <- ggplot(
    concordance_dt,
    aes(
        x = Hedges_g_mRNA,
        y = Hedges_g_Methylation
    )
) +
    # Very light emphasis of the prespecified concordant-positive direction.
    annotate(
        "rect",
        xmin = 0,
        xmax = Inf,
        ymin = 0,
        ymax = Inf,
        fill = "#FAF4F4",
        alpha = 0.35
    ) +
    geom_hline(
        yintercept = 0,
        linewidth = 0.30,
        color = "#969696",
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linewidth = 0.30,
        color = "#969696",
        linetype = "dashed"
    ) +
    geom_abline(
        intercept = 0,
        slope = 1,
        linewidth = 0.28,
        color = "#C0C0C0",
        linetype = "22"
    ) +
    # Plot nonsignificant/partial features first and lighter.
    geom_point(
        data = concordance_dt[
            Both_Significant == FALSE
        ],
        aes(
            color = Method,
            shape = Method
        ),
        size = 1.65,
        alpha = 0.38,
        stroke = 0.25
    ) +
    # Significant-in-both features are larger and fully opaque.
    geom_point(
        data = concordance_dt[
            Both_Significant == TRUE
        ],
        aes(
            color = Method,
            shape = Method
        ),
        size = 2.45,
        alpha = 1.00,
        stroke = 0.30
    ) +
    ggrepel::geom_label_repel(
        data = label_dt,
        aes(
            label = Label_Display
        ),
        size = pt_to_mm(4.9),
        family = FONT_FAMILY,
        color = INK,
        fill = scales::alpha(
            "white",
            0.82
        ),
        label.size = 0,
        label.padding = unit(
            0.07,
            "lines"
        ),
        box.padding = 0.28,
        point.padding = 0.16,
        min.segment.length = 0,
        segment.size = 0.18,
        segment.color = "#B5B5B5",
        max.overlaps = Inf,
        force = 2.6,
        force_pull = 0.12,
        max.time = 5,
        max.iter = 30000,
        seed = 20260913,
        show.legend = FALSE
    ) +
    scale_color_manual(
        values = METHOD_POINT_COLORS,
        breaks = method_order,
        name = NULL
    ) +
    scale_shape_manual(
        values = METHOD_SHAPES,
        breaks = method_order,
        name = NULL
    ) +
    scale_x_continuous(
        limits = c(
            -scatter_plot_limit,
            scatter_plot_limit
        ),
        breaks = scales::pretty_breaks(n = 5),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
        limits = c(
            -scatter_plot_limit,
            scatter_plot_limit
        ),
        breaks = scales::pretty_breaks(n = 5),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
        title = "Expression-methylation concordance",
        subtitle = rho_text,
        x = "Hedges' g: DNMT3A mRNA High vs Low",
        y = "Hedges' g: cg21629895 Low vs High"
    ) +
    coord_equal(
        clip = "off"
    ) +
    theme_classic(
        base_size = 6,
        base_family = FONT_FAMILY
    ) +
    theme(
        plot.title = element_text(
            size = FONT_TITLE,
            face = "bold",
            color = INK,
            hjust = 0,
            margin = margin(b = 1)
        ),
        plot.subtitle = element_text(
            size = 5.0,
            color = MUTED,
            hjust = 0,
            margin = margin(b = 3)
        ),
        axis.title = element_text(
            size = FONT_AXIS_TITLE,
            face = "bold",
            color = INK
        ),
        axis.text = element_text(
            size = FONT_AXIS_TEXT,
            color = INK
        ),
        axis.line = element_line(
            linewidth = 0.35,
            color = INK
        ),
        axis.ticks = element_line(
            linewidth = 0.30,
            color = INK
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(
            size = FONT_LEGEND
        ),
        legend.key.width = unit(
            3.2,
            "mm"
        ),
        legend.spacing.x = unit(
            1.5,
            "mm"
        ),
        plot.margin = margin(
            3,
            7,
            3,
            3
        )
    )

panel_B <- cowplot::ggdraw(
    panel_B_base
) +
    cowplot::draw_label(
        "B",
        x = 0.002,
        y = 0.998,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )


# ------------------------------------------------------------------------------
# 21. Panel C: method-wise concordance summary
# ------------------------------------------------------------------------------

summary_dt <- concordance_dt[
    ,
    .N,
    by = .(
        Method,
        Concordance_Class
    )
]

summary_grid <- data.table::CJ(
    Method = method_order,
    Concordance_Class = concordance_levels,
    unique = TRUE
)

summary_grid[
    ,
    Method := factor(
        Method,
        levels = method_order
    )
]

summary_grid[
    ,
    Concordance_Class := factor(
        Concordance_Class,
        levels = concordance_levels
    )
]

summary_dt[
    ,
    Method := factor(
        Method,
        levels = method_order
    )
]

summary_dt[
    ,
    Concordance_Class := factor(
        Concordance_Class,
        levels = concordance_levels
    )
]

summary_dt <- merge(
    summary_grid,
    summary_dt,
    by = c(
        "Method",
        "Concordance_Class"
    ),
    all.x = TRUE,
    sort = FALSE
)

summary_dt[
    is.na(N),
    N := 0L
]

summary_dt[
    ,
    Total := sum(N),
    by = Method
]

summary_dt[
    ,
    Percent := data.table::fifelse(
        Total > 0,
        100 * N / Total,
        0
    )
]

summary_dt[
    ,
    Percent_Label := data.table::fifelse(
        N > 0 &
            Percent >= 8,
        paste0(
            sprintf(
                "%.1f",
                Percent
            ),
            "%"
        ),
        ""
    )
]

method_labels <- c(
    "ESTIMATE" = "ESTIMATE\n(n = 3)",
    "MCP-counter" = "MCP-counter\n(n = 10)",
    "Cell-type ssGSEA" = paste0(
        "Cell-type ssGSEA\n(n = ",
        length(ssgsea_order),
        ")"
    )
)


panel_C_base <- ggplot(
    summary_dt,
    aes(
        x = Method,
        y = Percent,
        fill = Concordance_Class
    )
) +
    geom_col(
        width = 0.70,
        color = "white",
        linewidth = 0.35
    ) +
    geom_text(
        aes(
            label = Percent_Label
        ),
        position = position_stack(
            vjust = 0.5
        ),
        size = pt_to_mm(5.0),
        color = INK,
        fontface = "bold"
    ) +
    scale_fill_manual(
        values = CONCORDANCE_COLORS,
        breaks = concordance_levels,
        name = NULL
    ) +
    scale_x_discrete(
        labels = method_labels
    ) +
    scale_y_continuous(
        limits = c(
            0,
            100
        ),
        breaks = c(
            0,
            25,
            50,
            75,
            100
        ),
        expand = c(
            0,
            0
        )
    ) +
    labs(
        title = "Concordance summary by immune method",
        x = NULL,
        y = "Features (%)"
    ) +
    theme_classic(
        base_size = 6,
        base_family = FONT_FAMILY
    ) +
    theme(
        plot.title = element_text(
            size = FONT_TITLE,
            face = "bold",
            color = INK,
            hjust = 0,
            margin = margin(b = 2)
        ),
        axis.title.y = element_text(
            size = FONT_AXIS_TITLE,
            face = "bold",
            color = INK
        ),
        axis.text.y = element_text(
            size = FONT_AXIS_TEXT,
            color = INK
        ),
        axis.text.x = element_text(
            size = 5.0,
            color = INK,
            angle = 0,
            hjust = 0.5,
            vjust = 1,
            lineheight = 0.90
        ),
        axis.line = element_line(
            linewidth = 0.35,
            color = INK
        ),
        axis.ticks = element_line(
            linewidth = 0.30,
            color = INK
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.text = element_text(
            size = 5.0
        ),
        legend.key.height = unit(
            2.6,
            "mm"
        ),
        legend.key.width = unit(
            4.2,
            "mm"
        ),
        legend.margin = margin(
            t = 2,
            r = 0,
            b = 0,
            l = 0
        ),
        plot.margin = margin(
            2,
            4,
            7,
            4
        )
    ) +
    guides(
        fill = guide_legend(
            nrow = 2,
            byrow = TRUE
        )
    )

panel_C <- cowplot::ggdraw(
    panel_C_base
) +
    cowplot::draw_label(
        "C",
        x = 0.002,
        y = 0.998,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )


# ------------------------------------------------------------------------------
# 22. Final layout matching the requested manuscript composition
# ------------------------------------------------------------------------------

right_column <- cowplot::plot_grid(
    panel_B,
    panel_C,
    ncol = 1,
    rel_heights = c(
        0.58,
        0.42
    ),
    align = "v",
    axis = "lr"
)

main_body <- cowplot::plot_grid(
    panel_A,
    right_column,
    nrow = 1,
    rel_widths = c(
        0.53,
        0.47
    ),
    align = "hv",
    axis = "tblr"
)

main_title <- cowplot::ggdraw() +
    cowplot::draw_label(
        "DNMT3A expression and cg21629895 methylation are associated with coordinated immune/TME differences",
        x = 0.00,
        y = 0.82,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_TITLE,
        color = INK,
        fontfamily = FONT_FAMILY
    ) +
    cowplot::draw_label(
        "Positive effects indicate higher scores in DNMT3A-high or cg21629895-low tumors",
        x = 0.00,
        y = 0.24,
        hjust = 0,
        vjust = 1,
        size = 5.0,
        color = MUTED,
        fontfamily = FONT_FAMILY
    )

final_figure_4 <- cowplot::plot_grid(
    main_title,
    main_body,
    ncol = 1,
    rel_heights = c(
        0.075,
        0.925
    )
)


# ------------------------------------------------------------------------------
# 23. Save helpers
# ------------------------------------------------------------------------------

save_pdf <- function(
        plot_object,
        filename,
        width_mm,
        height_mm
) {

    grDevices::pdf(
        file = filename,
        width = width_mm / 25.4,
        height = height_mm / 25.4,
        family = "sans",
        useDingbats = FALSE
    )

    print(
        plot_object
    )

    grDevices::dev.off()
}

save_png <- function(
        plot_object,
        filename,
        width_mm,
        height_mm
) {

    ggplot2::ggsave(
        filename = filename,
        plot = plot_object,
        width = width_mm,
        height = height_mm,
        units = "mm",
        dpi = FIGURE_DPI,
        bg = "white",
        limitsize = FALSE
    )
}


# ------------------------------------------------------------------------------
# 24. Save final Figure 4
# ------------------------------------------------------------------------------

combined_pdf <- file.path(
    OUTPUT_DIR,
    "Fig4_Immune_TME_Concordance_FINAL_V6_SUBMISSION.pdf"
)

combined_png <- file.path(
    OUTPUT_DIR,
    "Fig4_Immune_TME_Concordance_FINAL_V6_SUBMISSION_600dpi.png"
)

save_pdf(
    final_figure_4,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

save_png(
    final_figure_4,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)


# ------------------------------------------------------------------------------
# 25. Save individual panels
# ------------------------------------------------------------------------------

save_pdf(
    panel_A,
    file.path(
        OUTPUT_DIR,
        "Fig4A_Immune_TME_HedgesG_Heatmap_FINAL_V2.pdf"
    ),
    100,
    136
)

save_pdf(
    panel_B,
    file.path(
        OUTPUT_DIR,
        "Fig4B_FullQuadrant_Concordance_FINAL_V6.pdf"
    ),
    80,
    66
)

save_pdf(
    panel_C,
    file.path(
        OUTPUT_DIR,
        "Fig4C_Concordance_Summary_FINAL_V2.pdf"
    ),
    80,
    66
)


# ------------------------------------------------------------------------------
# 26. Save statistics and QC
# ------------------------------------------------------------------------------

data.table::fwrite(
    all_stats_final,
    file.path(
        OUTPUT_DIR,
        "Fig4A_Immune_TME_EffectSize_Data_FINAL_V2.csv"
    ),
    na = ""
)

data.table::fwrite(
    concordance_dt,
    file.path(
        OUTPUT_DIR,
        "Fig4B_Immune_EffectSize_Concordance_Data_FINAL_V2.csv"
    ),
    na = ""
)

data.table::fwrite(
    summary_dt,
    file.path(
        OUTPUT_DIR,
        "Fig4C_Immune_Concordance_Summary_Data_FINAL_V2.csv"
    ),
    na = ""
)

concordance_qc <- data.table::data.table(
    Metric = c(
        "Total immune/TME features",
        "Features significant in both comparisons and labeled in Panel B",
        "All-feature Spearman rho",
        "All-feature Spearman P",
        "Concordant positive",
        "Concordant negative",
        "Discordant significant",
        "Partial or not significant",
        "xCell status",
        "Figure width mm",
        "Figure height mm"
    ),
    Value = c(
        nrow(concordance_dt),
        nrow(label_dt),
        rho_value_all,
        rho_p_all,
        sum(
            concordance_dt$Concordance_Class ==
                "Concordant positive"
        ),
        sum(
            concordance_dt$Concordance_Class ==
                "Concordant negative"
        ),
        sum(
            concordance_dt$Concordance_Class ==
                "Discordant significant"
        ),
        sum(
            concordance_dt$Concordance_Class ==
                "Partial / not significant"
        ),
        "Excluded",
        FIGURE_WIDTH_MM,
        FIGURE_HEIGHT_MM
    )
)

data.table::fwrite(
    concordance_qc,
    file.path(
        OUTPUT_DIR,
        "Fig4_Immune_TME_Concordance_QC_FINAL_V2.csv"
    ),
    na = ""
)


# ------------------------------------------------------------------------------
# 27. Console report
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FIGURE 4 V6: IMMUNE/TME CONCORDANCE\n")
cat("============================================================\n")

cat("\nMethods included:\n")
cat("ESTIMATE\n")
cat("MCP-counter\n")
cat("Cell-type ssGSEA using Bindea signatures\n")
cat("xCell: excluded\n")

cat("\nPanel A:\n")
cat("Hedges' g heatmap with compact horizontal scale bar\n")

cat("\nPanel B:\n")
cat("Positive-positive quadrant only\n")
cat("Displayed features =", nrow(label_dt), "\n")
cat("All-feature Spearman rho =", rho_value_all, "\n")
cat("All-feature P =", rho_p_all, "\n")

cat("\nPanel C:\n")
cat("Method-wise concordance summary for ALL features\n")

cat("\nFigure size:\n")
cat(FIGURE_WIDTH_MM, "mm x", FIGURE_HEIGHT_MM, "mm\n")

cat("\nMain PDF:\n")
cat(combined_pdf, "\n")

cat("\nMain PNG:\n")
cat(combined_png, "\n")

cat("\n============================================================\n")
