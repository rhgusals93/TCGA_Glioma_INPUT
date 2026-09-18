# =============================================================================
# TCGA-LGG DNMT3A Integrated Figure
# FINAL V9.3 - Reference-layout version
#
# Final panel structure:
#
# A. Concordant DEG heatmap
#
# B. Shared positive Hallmark enrichment programs
#    - DNMT3A mRNA High vs Low
#    - cg21629895 methylation Low vs High
#
# C. Shared positive VIPER transcription-factor programs
#    - DNMT3A mRNA High vs Low
#    - cg21629895 methylation Low vs High
#
# D. Continuous Hallmark associations
#    - Left subpanel: continuous DNMT3A mRNA expression
#    - Right subpanel: continuous cg21629895 methylation
#    - These are arranged like the TCGA / BeatAML pair in the reference.
#
# E. GO biological-process programs
#    - Upper GOChord: DNMT3A-high-associated positive DEGs
#    - Lower GOChord: cg21629895-low-associated positive DEGs
#    - Both GOChord plots are contained within one right-hand panel.
#
# Layout model:
#   A full width
#   B | C
#   D | E
#     |
#
# This follows the visual organization of the uploaded Figure 5 reference:
# two matched middle panels, one continuous-association panel at bottom left
# containing two side-by-side continuous plots, and one GO panel at bottom
# right containing two vertically stacked GOChord plots.
#
# Input files required: TWO
# 1. Completed TCGA-LGG DNMT3A master CSV
# 2. Annotated TCGA-LGG legacy HiSeqV2 RNA expression CSV
#
# RNA-seq unit:
#   log2(norm_count + 1)
#
# IMPORTANT:
#   - No additional log transformation is applied.
#   - limma is used for DEG analysis.
#   - Continuous Hallmark scores are calculated from the normalized continuous
#     expression matrix using GSVA with a Gaussian kernel.
#   - Spearman correlation is used for continuous associations.
#   - BH FDR is calculated separately for DNMT3A expression and cg21629895.
#
# Final publication format:
#   - 180 mm wide
#   - 248 mm high
#   - vector PDF + 600-dpi PNG
#   - standard grDevices::pdf()
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Package setup
# -----------------------------------------------------------------------------

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
    "pheatmap",
    "RColorBrewer",
    "msigdbr",
    "GOplot"
)

bioc_packages <- c(
    "limma",
    "clusterProfiler",
    "dorothea",
    "viper",
    "Biobase",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "GSVA"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(
            pkg,
            dependencies = TRUE
        )
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )
    }
}

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(stringr)
    library(ggplot2)
    library(cowplot)
    library(patchwork)
    library(scales)
    library(pheatmap)
    library(limma)
})


# -----------------------------------------------------------------------------
# 1. Global settings
# -----------------------------------------------------------------------------

set.seed(1234)

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 246
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"
FONT_PANEL <- 6.8
FONT_TITLE <- 6.2
FONT_AXIS_TITLE <- 5.6
FONT_AXIS_TEXT <- 5.2
FONT_LEGEND <- 5.0
FONT_LABEL <- 4.7
FONT_HEATMAP <- 4.6
FONT_CELL <- 4.2

INK <- "#222222"
MUTED <- "#666666"
GRID_LIGHT <- "#ECECEC"

MRNA_COLOR <- "#D1495B"
METH_COLOR <- "#3B82C4"

FDR_CUTOFF <- 0.05
DEG_ABS_LOG2FC_CUTOFF <- 0.50
HEATMAP_TOP_N <- 50L

TOP_SHARED_HALLMARK <- 13L
TOP_SHARED_VIPER <- 20L

DOROTHEA_LEVELS <- c("A", "B", "C")
VIPER_MIN_TARGETS <- 5L

# Match the uploaded Figure 5 V5 GOChord DEG threshold.
GO_DEG_FDR_CUTOFF <- 0.05
GO_DEG_LOG2FC_CUTOFF <- 0.50
GO_ORA_FDR_CUTOFF <- 0.05
GOCHORD_MAX_TERMS <- 6L
GOCHORD_MAX_GENES_PER_TERM <- 6L
GOCHORD_MAX_DISPLAY_GENES <- 20L

GOCHORD_GO_COLORS <- c(
    "#D83B32",
    "#E4E64F",
    "#69AE4B",
    "#66C4C8",
    "#3D6196",
    "#A15A9E"
)

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

standardize_gene <- function(x) {

    x <- toupper(
        trimws(
            as.character(x)
        )
    )

    x[x == ""] <- NA_character_

    x
}


# -----------------------------------------------------------------------------
# 2. Select input files
# -----------------------------------------------------------------------------

cat("\n1. Select the completed TCGA-LGG DNMT3A master CSV file.\n")
file_1_path <- file.choose()

cat("\n2. Select the annotated TCGA-LGG HiSeqV2 RNA expression CSV file.\n")
file_2_path <- file.choose()

file_1 <- data.table::fread(
    file_1_path,
    check.names = FALSE
)

file_2 <- data.table::fread(
    file_2_path,
    check.names = FALSE
)

names(file_1) <- trimws(
    sub(
        "^\ufeff",
        "",
        names(file_1)
    )
)

names(file_2) <- trimws(
    sub(
        "^\ufeff",
        "",
        names(file_2)
    )
)

master_signature <- c(
    "Patient_ID",
    "Sample_Core",
    "DNMT3A_expression",
    "DNMT3A_expression_Median_Group",
    "cg21629895",
    "cg21629895_Median_Group",
    "Molecular_3Group"
)

is_master_like <- function(x) {
    all(
        master_signature %in% names(x)
    )
}

is_rna_like <- function(x) {
    "Gene_Symbol" %in% names(x)
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

    cat("\nThe two selected files were reversed. They were reassigned automatically.\n")

    master <- file_2
    rna_raw <- file_1

} else {

    cat("\nFile 1 first 25 columns:\n")
    print(head(names(file_1), 25))

    cat("\nFile 2 first 25 columns:\n")
    print(head(names(file_2), 25))

    stop(
        paste0(
            "The selected files could not be identified as the required pair.\n",
            "Master must contain: ",
            paste(master_signature, collapse = ", "),
            "\nRNA file must contain Gene_Symbol."
        )
    )
}


# -----------------------------------------------------------------------------
# 3. Output directory
# -----------------------------------------------------------------------------

choose_output_dir <- function() {

    selected <- NULL

    if (.Platform$OS.type == "windows") {
        selected <- tryCatch(
            utils::choose.dir(
                default = Sys.getenv("USERPROFILE"),
                caption = "Select a SHORT output folder for the integrated DNMT3A figure"
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
            "DNMT3A_Integrated_Figure"
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

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")


# -----------------------------------------------------------------------------
# 4. Prepare master groups
# -----------------------------------------------------------------------------

clean_binary_group <- function(
        x,
        variable_name
) {

    y <- trimws(
        as.character(x)
    )

    low <- tolower(y)

    out <- rep(
        NA_character_,
        length(y)
    )

    out[
        low %in% c("low", "l")
    ] <- "Low"

    out[
        low %in% c("high", "h")
    ] <- "High"

    unexpected <- unique(
        y[
            !is.na(y) &
            nzchar(y) &
            is.na(out)
        ]
    )

    if (length(unexpected) > 0) {
        stop(
            paste0(
                "Unexpected values in ",
                variable_name,
                ": ",
                paste(unexpected, collapse = ", ")
            )
        )
    }

    out
}

master[
    ,
    Sample_Core := toupper(
        trimws(
            as.character(Sample_Core)
        )
    )
]

master[
    ,
    DNMT3A_Group := factor(
        clean_binary_group(
            DNMT3A_expression_Median_Group,
            "DNMT3A_expression_Median_Group"
        ),
        levels = c("Low", "High")
    )
]

master[
    ,
    cg21629895_Group := factor(
        clean_binary_group(
            cg21629895_Median_Group,
            "cg21629895_Median_Group"
        ),
        levels = c("High", "Low")
    )
]

master[
    ,
    Molecular_3Group := as.character(
        Molecular_3Group
    )
]

master[
    is.na(Molecular_3Group) |
    Molecular_3Group == "",
    Molecular_3Group := "Unclassified"
]


# -----------------------------------------------------------------------------
# 5. Build aligned expression matrix
# -----------------------------------------------------------------------------

if (!"Gene_Symbol" %in% names(rna_raw)) {
    stop("The RNA file must contain Gene_Symbol.")
}

tcga_sample_columns <- grep(
    "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
    names(rna_raw),
    value = TRUE
)

if (length(tcga_sample_columns) == 0) {
    stop(
        "No 15-character TCGA sample-core columns were detected in the RNA file."
    )
}

tcga_sample_columns <- tcga_sample_columns[
    substr(
        tcga_sample_columns,
        14,
        15
    ) == "01"
]

overlap_samples <- intersect(
    master$Sample_Core,
    tcga_sample_columns
)

if (length(overlap_samples) < 20) {
    stop(
        "Fewer than 20 overlapping primary tumor samples were detected."
    )
}

rna_dt <- rna_raw[
    !is.na(Gene_Symbol) &
    Gene_Symbol != "",
    c(
        "Gene_Symbol",
        overlap_samples
    ),
    with = FALSE
]

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

rna_long[
    ,
    Gene_Symbol := standardize_gene(
        Gene_Symbol
    )
]

rna_long <- rna_long[
    !is.na(Gene_Symbol) &
    is.finite(Expression)
]

rna_collapsed <- rna_long[
    ,
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

sample_order <- unique(
    master$Sample_Core[
        master$Sample_Core %in%
        colnames(expr_matrix)
    ]
)

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

if (!identical(
    analysis_master$Sample_Core,
    colnames(expr_matrix)
)) {
    stop(
        "Master/RNA sample alignment failed."
    )
}

gene_keep <- apply(
    expr_matrix,
    1,
    function(z) {

        z <- z[
            is.finite(z)
        ]

        length(z) >= 20 &&
        is.finite(stats::sd(z)) &&
        stats::sd(z) > 0
    }
)

expr_matrix <- expr_matrix[
    gene_keep,
    ,
    drop = FALSE
]

cat("\nAligned RNA matrix:\n")
cat("Genes:", nrow(expr_matrix), "\n")
cat("Samples:", ncol(expr_matrix), "\n")
cat("RNA unit: log2(norm_count + 1)\n")


# -----------------------------------------------------------------------------
# 6. limma contrasts shared by Panels A, C, D, E, and F
# -----------------------------------------------------------------------------

run_limma_contrast <- function(
        expression_matrix,
        group,
        positive_level,
        reference_level,
        analysis_label
) {

    group <- factor(
        as.character(group),
        levels = c(
            reference_level,
            positive_level
        )
    )

    keep <- !is.na(group)

    group_use <- droplevels(
        group[
            keep
        ]
    )

    if (nlevels(group_use) != 2) {
        stop(
            paste0(
                analysis_label,
                " does not contain two valid groups."
            )
        )
    }

    design <- stats::model.matrix(
        ~ 0 + group_use
    )

    colnames(design) <- levels(group_use)

    fit <- limma::lmFit(
        expression_matrix[
            ,
            keep,
            drop = FALSE
        ],
        design
    )

    contrast_matrix <- limma::makeContrasts(
        contrasts = paste0(
            positive_level,
            "-",
            reference_level
        ),
        levels = design
    )

    fit <- limma::contrasts.fit(
        fit,
        contrast_matrix
    )

    fit <- limma::eBayes(fit)

    result <- limma::topTable(
        fit,
        number = Inf,
        adjust.method = "BH",
        sort.by = "none"
    )

    result$Gene <- rownames(result)

    result <- data.table::as.data.table(
        result
    )

    result[
        ,
        Analysis := analysis_label
    ]

    list(
        result = result,
        fit = fit,
        group = group
    )
}

rna_limma <- run_limma_contrast(
    expression_matrix = expr_matrix,
    group = analysis_master$DNMT3A_Group,
    positive_level = "High",
    reference_level = "Low",
    analysis_label = "DNMT3A mRNA High vs Low"
)

meth_limma <- run_limma_contrast(
    expression_matrix = expr_matrix,
    group = analysis_master$cg21629895_Group,
    positive_level = "Low",
    reference_level = "High",
    analysis_label = "cg21629895 methylation Low vs High"
)


# -----------------------------------------------------------------------------
# 7. Panel A: Figure 2C-style concordant DEG heatmap
# -----------------------------------------------------------------------------

expr_deg <- rna_limma$result[
    ,
    .(
        Gene = standardize_gene(Gene),
        log2FC_DNMT3A = logFC,
        FDR_DNMT3A = adj.P.Val
    )
]

meth_deg <- meth_limma$result[
    ,
    .(
        Gene = standardize_gene(Gene),
        log2FC_cg21629895 = logFC,
        FDR_cg21629895 = adj.P.Val
    )
]

deg_integrated <- merge(
    expr_deg,
    meth_deg,
    by = "Gene",
    all = FALSE
)

deg_integrated[
    ,
    Significant_DNMT3A :=
        is.finite(FDR_DNMT3A) &
        FDR_DNMT3A < FDR_CUTOFF &
        is.finite(log2FC_DNMT3A) &
        abs(log2FC_DNMT3A) >= DEG_ABS_LOG2FC_CUTOFF
]

deg_integrated[
    ,
    Significant_cg21629895 :=
        is.finite(FDR_cg21629895) &
        FDR_cg21629895 < FDR_CUTOFF &
        is.finite(log2FC_cg21629895) &
        abs(log2FC_cg21629895) >= DEG_ABS_LOG2FC_CUTOFF
]

deg_integrated[
    ,
    Concordant_Direction :=
        sign(log2FC_DNMT3A) ==
        sign(log2FC_cg21629895)
]

concordant_genes <- deg_integrated[
    Significant_DNMT3A &
    Significant_cg21629895 &
    Concordant_Direction
]

if (nrow(concordant_genes) == 0) {
    stop(
        "No concordantly regulated genes were detected."
    )
}

concordant_genes[
    ,
    Combined_Absolute_Log2FC :=
        abs(log2FC_DNMT3A) +
        abs(log2FC_cg21629895)
]

data.table::setorder(
    concordant_genes,
    -Combined_Absolute_Log2FC
)

heatmap_gene_table <- head(
    concordant_genes,
    min(
        HEATMAP_TOP_N,
        nrow(concordant_genes)
    )
)

heatmap_genes <- heatmap_gene_table$Gene

expr_heatmap_matrix <- expr_matrix[
    heatmap_genes[
        heatmap_genes %in%
        rownames(expr_matrix)
    ],
    ,
    drop = FALSE
]

z_matrix <- t(
    apply(
        expr_heatmap_matrix,
        1,
        function(x) {

            x <- as.numeric(x)

            if (
                all(is.na(x)) ||
                stats::sd(
                    x,
                    na.rm = TRUE
                ) == 0
            ) {
                return(
                    rep(
                        0,
                        length(x)
                    )
                )
            }

            as.numeric(
                scale(x)
            )
        }
    )
)

colnames(z_matrix) <- colnames(expr_heatmap_matrix)
rownames(z_matrix) <- rownames(expr_heatmap_matrix)

z_matrix[
    z_matrix > 2.5
] <- 2.5

z_matrix[
    z_matrix < -2.5
] <- -2.5

sample_annotation <- data.table::data.table(
    Sample_Core = analysis_master$Sample_Core,
    DNMT3A_Group = analysis_master$DNMT3A_Group,
    cg21629895_Group = analysis_master$cg21629895_Group,
    Molecular_3Group = analysis_master$Molecular_3Group
)

molecular_levels <- c(
    "IDH-mutant/1p19q-codeleted",
    "IDH-mutant/1p19q-non-codeleted",
    "IDH-wildtype",
    "Unclassified"
)

sample_annotation[
    ,
    Molecular_3Group := factor(
        Molecular_3Group,
        levels = molecular_levels
    )
]

sample_annotation[
    ,
    DNMT3A_Order := data.table::fcase(
        DNMT3A_Group == "Low", 1L,
        DNMT3A_Group == "High", 2L,
        default = 9L
    )
]

sample_annotation[
    ,
    Methylation_Order := data.table::fcase(
        cg21629895_Group == "High", 1L,
        cg21629895_Group == "Low", 2L,
        default = 9L
    )
]

sample_annotation[
    ,
    Molecular_Order := as.integer(
        Molecular_3Group
    )
]

data.table::setorder(
    sample_annotation,
    DNMT3A_Order,
    Methylation_Order,
    Molecular_Order,
    Sample_Core
)

ordered_samples <- sample_annotation$Sample_Core

z_matrix <- z_matrix[
    ,
    ordered_samples,
    drop = FALSE
]

# Panel A annotation order from TOP to BOTTOM:
#   1. DNMT3A expression
#   2. cg21629895 methylation
#   3. Molecular subtype
#
# pheatmap draws annotation rows in reverse column order, so the data.frame
# columns are intentionally supplied in reverse here to obtain the requested
# visual top-to-bottom order.
annotation_col <- as.data.frame(
    sample_annotation[
        ,
        .(
            `Molecular subtype` = Molecular_3Group,
            `cg21629895 methylation` = cg21629895_Group,
            `DNMT3A expression` = DNMT3A_Group
        )
    ]
)

rownames(annotation_col) <- ordered_samples

MOLECULAR_COLORS <- c(
    "IDH-mutant/1p19q-codeleted" = "#4C78A8",
    "IDH-mutant/1p19q-non-codeleted" = "#E39C37",
    "IDH-wildtype" = "#C94C4C",
    "Unclassified" = "#9E9E9E"
)

annotation_colors <- list(
    "Molecular subtype" = MOLECULAR_COLORS,
    "cg21629895 methylation" = c(
        "High" = "#D8ECEA",
        "Low" = "#1F6F78"
    ),
    "DNMT3A expression" = c(
        "Low" = "#E6E6E6",
        "High" = "#3A3A3A"
    )
)

heatmap_colors <- grDevices::colorRampPalette(
    c(
        "#2166AC",
        "#F7F7F7",
        "#B2182B"
    )
)(
    101
)

heatmap_breaks <- seq(
    -2.5,
    2.5,
    length.out = 102
)

panel_A_heatmap_object <- pheatmap::pheatmap(
    z_matrix,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    border_color = NA,
    scale = "none",
    fontsize = 5.0,
    fontsize_row = 4.2,
    fontsize_col = 4.3,
    treeheight_row = 14,
    treeheight_col = 0,
    main = paste0(
        "Concordant DEGs: ",
        nrow(concordant_genes),
        " total; top ",
        nrow(z_matrix),
        " shown"
    ),
    silent = TRUE
)

panel_A_grob <- grid::grobTree(
    panel_A_heatmap_object$gtable,
    grid::textGrob(
        "A",
        x = grid::unit(0.005, "npc"),
        y = grid::unit(0.995, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(
            fontfamily = FONT_FAMILY,
            fontface = "bold",
            fontsize = FONT_PANEL
        )
    )
)

panel_A <- cowplot::ggdraw() +
    cowplot::draw_grob(
        panel_A_grob,
        x = 0,
        y = 0,
        width = 1,
        height = 1
    )


# -----------------------------------------------------------------------------
# 8. Immune/TME analysis
# -----------------------------------------------------------------------------
# The previous integrated TME heatmap is intentionally omitted from this
# figure. It is reserved for the subsequent immune/TME figure.


# -----------------------------------------------------------------------------
# 9. Hallmark GSEA for Panel B
# -----------------------------------------------------------------------------

get_hallmark_term2gene <- function() {

    hallmark <- tryCatch(
        msigdbr::msigdbr(
            db_species = "HS",
            species = "Homo sapiens",
            collection = "H"
        ),
        error = function(e) {
            msigdbr::msigdbr(
                species = "Homo sapiens",
                category = "H"
            )
        }
    )

    hallmark <- as.data.frame(
        hallmark
    )

    hallmark[
        ,
        c(
            "gs_name",
            "gene_symbol"
        )
    ] |>
        dplyr::rename(
            term = gs_name,
            gene = gene_symbol
        ) |>
        dplyr::filter(
            !is.na(term),
            !is.na(gene),
            term != "",
            gene != ""
        ) |>
        dplyr::distinct()
}

clean_hallmark_name <- function(x) {

    x <- gsub(
        "^HALLMARK_",
        "",
        x
    )

    x <- gsub(
        "_",
        " ",
        x
    )

    tools::toTitleCase(
        tolower(x)
    )
}

hallmark_term2gene <- get_hallmark_term2gene()

run_hallmark_gsea <- function(
        limma_table,
        label
) {

    rank_vector <- limma_table$t
    names(rank_vector) <- limma_table$Gene

    rank_vector <- rank_vector[
        is.finite(rank_vector) &
        !is.na(names(rank_vector)) &
        names(rank_vector) != ""
    ]

    rank_dt <- data.table::data.table(
        Gene = names(rank_vector),
        Statistic = as.numeric(rank_vector)
    )

    rank_dt <- rank_dt[
        order(
            -abs(Statistic)
        )
    ]

    rank_dt <- rank_dt[
        !duplicated(Gene)
    ]

    rank_vector <- rank_dt$Statistic
    names(rank_vector) <- rank_dt$Gene

    rank_vector <- sort(
        rank_vector,
        decreasing = TRUE
    )

    gsea_object <- clusterProfiler::GSEA(
        geneList = rank_vector,
        TERM2GENE = hallmark_term2gene,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 0,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        seed = TRUE,
        verbose = FALSE,
        by = "fgsea"
    )

    result <- as.data.frame(
        gsea_object
    )

    if (nrow(result) == 0) {
        stop(
            paste0(
                "No Hallmark GSEA results were returned for ",
                label,
                "."
            )
        )
    }

    result$Analysis <- label
    result$Pathway <- clean_hallmark_name(
        result$ID
    )

    data.table::as.data.table(
        result
    )
}

gsea_rna <- run_hallmark_gsea(
    rna_limma$result,
    "DNMT3A mRNA High vs Low"
)

gsea_meth <- run_hallmark_gsea(
    meth_limma$result,
    "cg21629895 methylation Low vs High"
)

gsea_integrated <- merge(
    gsea_rna[
        ,
        .(
            Pathway,
            RNA_NES = NES,
            RNA_FDR = p.adjust
        )
    ],
    gsea_meth[
        ,
        .(
            Pathway,
            Methylation_NES = NES,
            Methylation_FDR = p.adjust
        )
    ],
    by = "Pathway",
    all = FALSE
)

gsea_integrated[
    ,
    Both_FDR_005 :=
        RNA_FDR < 0.05 &
        Methylation_FDR < 0.05
]

gsea_integrated[
    ,
    Worst_FDR := pmax(
        RNA_FDR,
        Methylation_FDR
    )
]

gsea_integrated[
    ,
    Mean_NES := (
        RNA_NES +
        Methylation_NES
    ) / 2
]

# IMPORTANT: Positive shared programs only.
panel_C_data <- gsea_integrated[
    is.finite(RNA_NES) &
    is.finite(Methylation_NES) &
    RNA_NES > 0 &
    Methylation_NES > 0
][
    order(
        -Both_FDR_005,
        Worst_FDR,
        -Mean_NES
    )
][
    1:min(
        TOP_SHARED_HALLMARK,
        .N
    )
]

if (nrow(panel_C_data) == 0) {
    stop(
        "No shared positive Hallmark pathways were available."
    )
}

panel_C_data <- panel_C_data[
    order(
        Mean_NES
    )
]

panel_C_data[
    ,
    Pathway := factor(
        Pathway,
        levels = Pathway
    )
]

panel_C_long <- data.table::rbindlist(
    list(
        panel_C_data[
            ,
            .(
                Pathway,
                Comparison = "DNMT3A mRNA (High vs Low)",
                Value = RNA_NES,
                FDR = RNA_FDR
            )
        ],
        panel_C_data[
            ,
            .(
                Pathway,
                Comparison = "cg21629895 methylation (Low vs High)",
                Value = Methylation_NES,
                FDR = Methylation_FDR
            )
        ]
    )
)

panel_C_long[
    ,
    Significant := factor(
        ifelse(
            FDR < 0.05,
            "FDR < 0.05",
            "FDR >= 0.05"
        ),
        levels = c(
            "FDR < 0.05",
            "FDR >= 0.05"
        )
    )
]

hallmark_x_max <- max(
    panel_C_long$Value,
    na.rm = TRUE
) * 1.08

panel_C_base <- ggplot() +
    annotate(
        "rect",
        xmin = 0,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = MRNA_COLOR,
        alpha = 0.025
    ) +
    geom_segment(
        data = panel_C_data,
        aes(
            x = RNA_NES,
            xend = Methylation_NES,
            y = Pathway,
            yend = Pathway
        ),
        color = "#B8B8B8",
        linewidth = 0.55
    ) +
    geom_point(
        data = panel_C_long,
        aes(
            x = Value,
            y = Pathway,
            fill = Comparison,
            alpha = Significant
        ),
        shape = 21,
        size = 2.20,
        color = "#222222",
        stroke = 0.40
    ) +
    scale_fill_manual(
        values = c(
            "DNMT3A mRNA (High vs Low)" = MRNA_COLOR,
            "cg21629895 methylation (Low vs High)" = METH_COLOR
        ),
        labels = c(
            "DNMT3A mRNA (High vs Low)" = "DNMT3A mRNA\n(High vs Low)",
            "cg21629895 methylation (Low vs High)" = "cg21629895 methylation\n(Low vs High)"
        ),
        name = NULL
    ) +
    scale_alpha_manual(
        values = c(
            "FDR < 0.05" = 1.00,
            "FDR >= 0.05" = 0.30
        ),
        name = "Significance"
    ) +
    scale_x_continuous(
        limits = c(
            0,
            hallmark_x_max
        ),
        breaks = scales::pretty_breaks(
            n = 5
        ),
        expand = c(
            0,
            0
        )
    ) +
    labs(
        title = "Shared positive Hallmark enrichment programs",
        subtitle = "Positive shared programs only\nDNMT3A mRNA: High vs Low | cg21629895 methylation: Low vs High",
        x = "Normalized enrichment score (NES)",
        y = NULL
    ) +
    theme_classic(
        base_size = 6,
        base_family = FONT_FAMILY
    ) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = FONT_TITLE,
            color = INK,
            margin = margin(b = 2)
        ),
        plot.subtitle = element_text(
            size = 4.45,
            color = MUTED,
            lineheight = 0.94,
            margin = margin(b = 4)
        ),
        axis.title.x = element_text(
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
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(
            size = FONT_LEGEND,
            face = "bold"
        ),
        legend.text = element_text(
            size = FONT_LEGEND
        ),
        legend.key.size = unit(
            2.0,
            "mm"
        ),
        plot.margin = margin(
            4,
            5,
            6,
            4
        )
    )

# Panel C and Panel D are aligned to the same left/right plot geometry below.
panel_C <- NULL


# -----------------------------------------------------------------------------
# 10. VIPER for Panel C
# -----------------------------------------------------------------------------

load_dorothea_regulon <- function() {

    regulon_env <- new.env(
        parent = emptyenv()
    )

    suppressWarnings(
        utils::data(
            "dorothea_hs",
            package = "dorothea",
            envir = regulon_env
        )
    )

    if (!exists(
        "dorothea_hs",
        envir = regulon_env,
        inherits = FALSE
    )) {
        stop(
            "Could not load dorothea_hs."
        )
    }

    regulon <- as.data.frame(
        get(
            "dorothea_hs",
            envir = regulon_env,
            inherits = FALSE
        )
    )

    regulon |>
        dplyr::filter(
            confidence %in%
            DOROTHEA_LEVELS
        ) |>
        dplyr::filter(
            !is.na(tf),
            !is.na(target),
            !is.na(mor)
        ) |>
        dplyr::distinct()
}

regulon <- load_dorothea_regulon()

run_viper_activity <- function(
        expression_matrix,
        regulon_table
) {

    result <- tryCatch(
        dorothea::run_viper(
            input = expression_matrix,
            regulons = regulon_table,
            options = list(
                minsize = VIPER_MIN_TARGETS
            ),
            tidy = FALSE
        ),
        error = function(e1) {
            dorothea::run_viper(
                input = expression_matrix,
                regulons = regulon_table,
                tidy = FALSE
            )
        }
    )

    if (inherits(result, "ExpressionSet")) {
        result <- Biobase::exprs(
            result
        )
    }

    result <- as.matrix(
        result
    )

    storage.mode(result) <- "numeric"

    sample_ids <- colnames(
        expression_matrix
    )

    if (
        sum(
            colnames(result) %in%
            sample_ids
        ) >= 5
    ) {

        activity <- result

    } else if (
        sum(
            rownames(result) %in%
            sample_ids
        ) >= 5
    ) {

        activity <- t(result)

    } else {

        stop(
            "Could not determine VIPER output orientation."
        )
    }

    activity <- activity[
        ,
        intersect(
            sample_ids,
            colnames(activity)
        ),
        drop = FALSE
    ]

    variable_tf <- apply(
        activity,
        1,
        function(z) {

            z <- z[
                is.finite(z)
            ]

            length(z) >= 20 &&
            is.finite(stats::sd(z)) &&
            stats::sd(z) > 0
        }
    )

    activity[
        variable_tf,
        ,
        drop = FALSE
    ]
}

viper_activity <- run_viper_activity(
    expr_matrix,
    regulon
)

viper_sample_order <- analysis_master$Sample_Core[
    analysis_master$Sample_Core %in%
    colnames(viper_activity)
]

viper_activity <- viper_activity[
    ,
    viper_sample_order,
    drop = FALSE
]

viper_master <- analysis_master[
    match(
        viper_sample_order,
        Sample_Core
    )
]

viper_rna_limma <- run_limma_contrast(
    expression_matrix = viper_activity,
    group = viper_master$DNMT3A_Group,
    positive_level = "High",
    reference_level = "Low",
    analysis_label = "DNMT3A mRNA High vs Low"
)

viper_meth_limma <- run_limma_contrast(
    expression_matrix = viper_activity,
    group = viper_master$cg21629895_Group,
    positive_level = "Low",
    reference_level = "High",
    analysis_label = "cg21629895 methylation Low vs High"
)

viper_integrated <- merge(
    viper_rna_limma$result[
        ,
        .(
            TF = Gene,
            RNA_Effect = logFC,
            RNA_FDR = adj.P.Val
        )
    ],
    viper_meth_limma$result[
        ,
        .(
            TF = Gene,
            Methylation_Effect = logFC,
            Methylation_FDR = adj.P.Val
        )
    ],
    by = "TF",
    all = FALSE
)

viper_integrated[
    ,
    Both_FDR_005 :=
        RNA_FDR < 0.05 &
        Methylation_FDR < 0.05
]

viper_integrated[
    ,
    Worst_FDR := pmax(
        RNA_FDR,
        Methylation_FDR
    )
]

viper_integrated[
    ,
    Mean_Effect := (
        RNA_Effect +
        Methylation_Effect
    ) / 2
]

# IMPORTANT: Positive shared TF programs only.
panel_D_data <- viper_integrated[
    is.finite(RNA_Effect) &
    is.finite(Methylation_Effect) &
    RNA_Effect > 0 &
    Methylation_Effect > 0
][
    order(
        -Both_FDR_005,
        Worst_FDR,
        -Mean_Effect
    )
][
    1:min(
        TOP_SHARED_VIPER,
        .N
    )
]

if (nrow(panel_D_data) == 0) {
    stop(
        "No shared positive VIPER TF programs were available."
    )
}

panel_D_data <- panel_D_data[
    order(
        Mean_Effect
    )
]

panel_D_data[
    ,
    TF := factor(
        TF,
        levels = TF
    )
]

panel_D_long <- data.table::rbindlist(
    list(
        panel_D_data[
            ,
            .(
                TF,
                Comparison = "DNMT3A mRNA (High vs Low)",
                Value = RNA_Effect,
                FDR = RNA_FDR
            )
        ],
        panel_D_data[
            ,
            .(
                TF,
                Comparison = "cg21629895 methylation (Low vs High)",
                Value = Methylation_Effect,
                FDR = Methylation_FDR
            )
        ]
    )
)

panel_D_long[
    ,
    Significant := factor(
        ifelse(
            FDR < 0.05,
            "FDR < 0.05",
            "FDR >= 0.05"
        ),
        levels = c(
            "FDR < 0.05",
            "FDR >= 0.05"
        )
    )
]

viper_x_max <- max(
    panel_D_long$Value,
    na.rm = TRUE
) * 1.08

panel_D_base <- ggplot() +
    annotate(
        "rect",
        xmin = 0,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = MRNA_COLOR,
        alpha = 0.025
    ) +
    geom_segment(
        data = panel_D_data,
        aes(
            x = RNA_Effect,
            xend = Methylation_Effect,
            y = TF,
            yend = TF
        ),
        color = "#B8B8B8",
        linewidth = 0.55
    ) +
    geom_point(
        data = panel_D_long,
        aes(
            x = Value,
            y = TF,
            fill = Comparison,
            alpha = Significant
        ),
        shape = 21,
        size = 2.20,
        color = "#222222",
        stroke = 0.40
    ) +
    scale_fill_manual(
        values = c(
            "DNMT3A mRNA (High vs Low)" = MRNA_COLOR,
            "cg21629895 methylation (Low vs High)" = METH_COLOR
        ),
        labels = c(
            "DNMT3A mRNA (High vs Low)" = "DNMT3A mRNA\n(High vs Low)",
            "cg21629895 methylation (Low vs High)" = "cg21629895 methylation\n(Low vs High)"
        ),
        name = NULL
    ) +
    scale_alpha_manual(
        values = c(
            "FDR < 0.05" = 1.00,
            "FDR >= 0.05" = 0.30
        ),
        name = "Significance"
    ) +
    scale_x_continuous(
        limits = c(
            0,
            viper_x_max
        ),
        breaks = scales::pretty_breaks(
            n = 5
        ),
        expand = c(
            0,
            0
        )
    ) +
    labs(
        title = "Shared positive VIPER transcription-factor programs",
        subtitle = "Positive shared TF programs only\nDNMT3A mRNA: High vs Low | cg21629895 methylation: Low vs High",
        x = "VIPER activity difference",
        y = NULL
    ) +
    theme_classic(
        base_size = 6,
        base_family = FONT_FAMILY
    ) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = FONT_TITLE,
            color = INK,
            margin = margin(b = 2)
        ),
        plot.subtitle = element_text(
            size = 4.45,
            color = MUTED,
            lineheight = 0.94,
            margin = margin(b = 4)
        ),
        axis.title.x = element_text(
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
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(
            size = FONT_LEGEND,
            face = "bold"
        ),
        legend.text = element_text(
            size = FONT_LEGEND
        ),
        legend.key.size = unit(
            2.0,
            "mm"
        ),
        plot.margin = margin(
            4,
            5,
            6,
            4
        )
    )

aligned_cd <- cowplot::align_plots(
    panel_C_base,
    panel_D_base,
    align = "v",
    axis = "lr"
)

panel_C <- cowplot::ggdraw(
    aligned_cd[[1]]
) +
    cowplot::draw_label(
        "C",
        x = 0.004,
        y = 0.996,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )

panel_D <- cowplot::ggdraw(
    aligned_cd[[2]]
) +
    cowplot::draw_label(
        "D",
        x = 0.004,
        y = 0.996,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )



# -----------------------------------------------------------------------------
# Final manuscript wrappers for Panels B and C
# -----------------------------------------------------------------------------

# Build one shared centered legend for Panels B and C.
# Both panels use the same molecular comparison and significance coding.
shared_bc_legend <- cowplot::get_legend(
    panel_C_base +
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.box.just = "center",
            legend.justification = "center"
        )
)

panel_B_base_nolegend <- panel_C_base +
    theme(
        legend.position = "none"
    )

panel_C_base_nolegend <- panel_D_base +
    theme(
        legend.position = "none"
    )

# Align the plotting regions so the x-axis length, y-axis height,
# and outer panel geometry are the same in B and C.
aligned_bc_base <- cowplot::align_plots(
    panel_B_base_nolegend,
    panel_C_base_nolegend,
    align = "hv",
    axis = "tblr"
)

panel_B_final <- cowplot::ggdraw(
    aligned_bc_base[[1]]
) +
    cowplot::draw_label(
        "B",
        x = 0.006,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )

panel_C_final <- cowplot::ggdraw(
    aligned_bc_base[[2]]
) +
    cowplot::draw_label(
        "C",
        x = 0.006,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )


# -----------------------------------------------------------------------------
# 11. GO:BP ORA and GOChord for Panels D and E
# -----------------------------------------------------------------------------

get_go_bp_term2gene <- function() {

    go_df <- tryCatch(
        msigdbr::msigdbr(
            db_species = "HS",
            species = "Homo sapiens",
            collection = "C5",
            subcollection = "GO:BP"
        ),
        error = function(e) {
            msigdbr::msigdbr(
                species = "Homo sapiens",
                category = "C5",
                subcategory = "GO:BP"
            )
        }
    )

    go_df <- as.data.frame(
        go_df
    )

    data.table::data.table(
        term = as.character(
            go_df$gs_name
        ),
        gene = as.character(
            go_df$gene_symbol
        )
    )[
        !is.na(term) &
        !is.na(gene) &
        term != "" &
        gene != ""
    ][
        !duplicated(
            paste(
                term,
                gene,
                sep = "|||"
            )
        )
    ]
}

clean_go_process <- function(x) {

    x <- sub(
        "^GOBP_",
        "",
        x
    )

    x <- sub(
        "^GO_",
        "",
        x
    )

    x <- gsub(
        "_",
        " ",
        x
    )

    stringr::str_to_sentence(
        tolower(x)
    )
}

go_bp_term2gene <- get_go_bp_term2gene()

run_go_bp_ora_for_chord <- function(
        limma_table,
        analysis_label
) {

    deg <- data.table::copy(
        limma_table
    )

    deg <- deg[
        !is.na(Gene) &
        Gene != "" &
        is.finite(logFC) &
        is.finite(adj.P.Val)
    ]

    deg[
        ,
        Abs_logFC := abs(logFC)
    ]

    data.table::setorder(
        deg,
        Gene,
        -Abs_logFC
    )

    deg <- deg[
        !duplicated(Gene)
    ]

    universe <- unique(
        deg$Gene
    )

    positive_genes <- unique(
        deg[
            adj.P.Val < GO_DEG_FDR_CUTOFF &
            logFC >= GO_DEG_LOG2FC_CUTOFF,
            Gene
        ]
    )

    negative_genes <- unique(
        deg[
            adj.P.Val < GO_DEG_FDR_CUTOFF &
            logFC <= -GO_DEG_LOG2FC_CUTOFF,
            Gene
        ]
    )

    run_one <- function(
            genes,
            direction_label
    ) {

        if (length(genes) < 5) {
            return(
                data.table::data.table()
            )
        }

        fit <- clusterProfiler::enricher(
            gene = genes,
            universe = universe,
            TERM2GENE = as.data.frame(
                go_bp_term2gene
            ),
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 10,
            maxGSSize = 500
        )

        result <- as.data.frame(
            fit
        )

        if (nrow(result) == 0) {
            return(
                data.table::data.table()
            )
        }

        result <- data.table::as.data.table(
            result
        )

        result[
            ,
            Direction := direction_label
        ]

        result
    }

    positive_ora <- run_one(
        positive_genes,
        "Positive"
    )

    negative_ora <- run_one(
        negative_genes,
        "Negative"
    )

    all_ora <- data.table::rbindlist(
        list(
            positive_ora,
            negative_ora
        ),
        fill = TRUE
    )

    significant_ora <- all_ora[
        is.finite(p.adjust) &
        p.adjust < GO_ORA_FDR_CUTOFF
    ][
        order(
            p.adjust,
            -Count
        )
    ]

    if (nrow(significant_ora) == 0) {
        stop(
            paste0(
                "No GO:BP terms reached FDR < ",
                GO_ORA_FDR_CUTOFF,
                " for ",
                analysis_label,
                " using DEG FDR < ",
                GO_DEG_FDR_CUTOFF,
                " and |log2FC| >= ",
                GO_DEG_LOG2FC_CUTOFF,
                "."
            )
        )
    }

    positive_sig <- significant_ora[
        Direction == "Positive"
    ]

    negative_sig <- significant_ora[
        Direction == "Negative"
    ]

    # FINAL FIGURE RULE:
    # Use the positive DEG direction only so that Panels C-F all describe
    # the same biological direction.
    #
    # DNMT3A contrast: High - Low
    #   Positive log2FC = higher in DNMT3A-high tumors.
    #
    # cg21629895 contrast: Low - High
    #   Positive log2FC = higher in cg21629895-low tumors.
    #
    # Negative-direction ORA results are retained in the ALL output table
    # for transparency, but they are not used to construct the GOChord.
    if (nrow(positive_sig) == 0) {
        stop(
            paste0(
                "No positively directed GO:BP terms reached FDR < ",
                GO_ORA_FDR_CUTOFF,
                " for ",
                analysis_label,
                ". The final integrated figure is configured to display ",
                "positive-direction GOChord results only."
            )
        )
    }

    selected_direction <- "Positive"
    selected_terms <- positive_sig
    selected_genes <- positive_genes

    selected_terms <- selected_terms[
        1:min(
            GOCHORD_MAX_TERMS,
            .N
        )
    ]

    list(
        deg = deg,
        all_ora = all_ora,
        significant_ora = significant_ora,
        selected_direction = selected_direction,
        selected_terms = selected_terms,
        selected_genes = selected_genes
    )
}

prep_gochord_data <- function(
        ora_object,
        panel_label
) {

    deg <- data.table::copy(
        ora_object$deg
    )

    selected_terms <- data.table::copy(
        ora_object$selected_terms
    )

    selected_genes <- ora_object$selected_genes

    term_gene_all <- go_bp_term2gene[
        term %in% selected_terms$ID &
        gene %in% selected_genes
    ]

    term_gene_all <- merge(
        term_gene_all,
        deg[
            ,
            .(
                SYMBOL = Gene,
                logFC
            )
        ],
        by.x = "gene",
        by.y = "SYMBOL",
        all.x = TRUE,
        sort = FALSE
    )

    term_gene_all <- term_gene_all[
        is.finite(logFC)
    ]

    term_gene_display <- term_gene_all[
        order(
            term,
            -abs(logFC)
        )
    ][
        ,
        head(
            .SD,
            GOCHORD_MAX_GENES_PER_TERM
        ),
        by = term
    ]

    gene_priority <- term_gene_display[
        ,
        .(
            Max_Abs_logFC = max(
                abs(logFC),
                na.rm = TRUE
            ),
            N_Terms = data.table::uniqueN(
                term
            )
        ),
        by = gene
    ][
        order(
            -N_Terms,
            -Max_Abs_logFC
        )
    ]

    display_genes <- head(
        gene_priority$gene,
        GOCHORD_MAX_DISPLAY_GENES
    )

    term_gene_display <- term_gene_display[
        gene %in% display_genes
    ]

    selected_terms <- selected_terms[
        ID %in%
        unique(
            term_gene_display$term
        )
    ]

    gene_fc <- unique(
        term_gene_display[
            ,
            .(
                SYMBOL = gene,
                logFC
            )
        ]
    )

    symbol_map <- suppressMessages(
        clusterProfiler::bitr(
            unique(
                gene_fc$SYMBOL
            ),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db::org.Hs.eg.db
        )
    )

    symbol_map <- data.table::as.data.table(
        symbol_map
    )

    symbol_map <- symbol_map[
        !duplicated(SYMBOL)
    ]

    gene_fc_mapped <- merge(
        gene_fc,
        symbol_map,
        by = "SYMBOL",
        all = FALSE,
        sort = FALSE
    )

    gene_fc_mapped <- gene_fc_mapped[
        !duplicated(ENTREZID)
    ]

    term_gene_entrez <- merge(
        term_gene_display,
        symbol_map,
        by.x = "gene",
        by.y = "SYMBOL",
        all = FALSE,
        sort = FALSE
    )

    go_for_plot <- data.table::rbindlist(
        lapply(
            seq_len(
                nrow(selected_terms)
            ),
            function(i) {

                term_id <- selected_terms$ID[i]

                genes_here <- unique(
                    term_gene_entrez[
                        term == term_id,
                        ENTREZID
                    ]
                )

                if (length(genes_here) == 0) {
                    return(NULL)
                }

                data.table::data.table(
                    Category = "BP",
                    ID = term_id,
                    Term = clean_go_process(
                        term_id
                    ),
                    Genes = paste(
                        genes_here,
                        collapse = ", "
                    ),
                    adj_pval = selected_terms$p.adjust[i]
                )
            }
        ),
        fill = TRUE
    )

    if (nrow(go_for_plot) == 0) {
        stop(
            paste0(
                "GOChord preparation failed for ",
                panel_label,
                "."
            )
        )
    }

    # Use the signed gene log2FC directly in GOChord, matching the
    # uploaded reference Figure 5 + GO plot script.
    genelist <- data.frame(
        ID = gene_fc_mapped$ENTREZID,
        logFC = gene_fc_mapped$logFC,
        stringsAsFactors = FALSE
    )

    go_table_df <- as.data.frame(
        go_for_plot
    )

    circ <- GOplot::circle_dat(
        go_table_df,
        genelist
    )

    chord <- GOplot::chord_dat(
        circ,
        genelist,
        go_table_df$Term
    )

    chord_symbol <- chord

    rn <- rownames(
        chord_symbol
    )

    if (!is.null(rn)) {

        entrez_to_symbol <- stats::setNames(
            symbol_map$SYMBOL,
            symbol_map$ENTREZID
        )

        new_rn <- unname(
            entrez_to_symbol[
                rn
            ]
        )

        replace_idx <- is.na(new_rn) |
            !nzchar(new_rn)

        new_rn[
            replace_idx
        ] <- rn[
            replace_idx
        ]

        rownames(chord_symbol) <- make.unique(
            new_rn
        )
    }

    list(
        chord = chord_symbol,
        go_table = go_table_df,
        gene_fc = gene_fc_mapped,
        direction = ora_object$selected_direction,
        term_gene_all = term_gene_all,
        term_gene_display = term_gene_display
    )
}

go_expr_ora <- run_go_bp_ora_for_chord(
    rna_limma$result,
    "DNMT3A mRNA High vs Low"
)

go_meth_ora <- run_go_bp_ora_for_chord(
    meth_limma$result,
    "cg21629895 methylation Low vs High"
)

go_expr_data <- prep_gochord_data(
    go_expr_ora,
    "DNMT3A expression"
)

go_meth_data <- prep_gochord_data(
    go_meth_ora,
    "cg21629895 methylation"
)

all_lfc <- c(
    go_expr_data$gene_fc$logFC,
    go_meth_data$gene_fc$logFC
)

all_lfc <- all_lfc[
    is.finite(all_lfc)
]

if (length(all_lfc) == 0) {
    stop(
        "No finite GOChord gene log2FC values were available."
    )
}

# Because Panels E and F are now restricted to positive-direction DEGs,
# the displayed gene log2FC values should all be positive.
if (any(all_lfc < 0, na.rm = TRUE)) {
    warning(
        "Negative GOChord log2FC values were detected despite positive-direction filtering."
    )
}

# Use one common publication-friendly 0-based scale for Panels E and F.
# The upper bound is rounded upward to the nearest 0.5.
lfc_max_raw <- max(
    all_lfc,
    na.rm = TRUE
)

lfc_max_common <- ceiling(
    lfc_max_raw * 2
) / 2

if (
    !is.finite(lfc_max_common) ||
    lfc_max_common <= 0
) {
    lfc_max_common <- 1
}

lfc_limits <- c(
    0,
    lfc_max_common
)

lfc_breaks <- c(
    0,
    round(
        lfc_max_common / 2,
        1
    ),
    lfc_max_common
)

GOCHORD_LFC_COLORS <- c(
    "#8F6A27",
    "#B9882E",
    "#E3AD42"
)


make_chord_plot <- function(
        chord_mat,
        n_terms,
        title_text
) {

    GOplot::GOChord(
        chord_mat,
        space = 0.025,
        gene.order = "logFC",
        lfc.col = GOCHORD_LFC_COLORS,
        ribbon.col = GOCHORD_GO_COLORS[
            seq_len(n_terms)
        ],
        gene.space = 0.42,
        gene.size = 1.50,
        border.size = 0.08,
        process.label = 0
    ) +
        labs(
            title = title_text
        ) +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold",
                size = 5.8,
                margin = margin(
                    b = 1
                )
            ),
            plot.margin = margin(
                t = 2,
                r = 4,
                b = 2,
                l = 4,
                unit = "pt"
            ),
            legend.position = "none",
            aspect.ratio = 1
        )
}

make_term_legend_compact <- function(
        go_table
) {

    df <- data.table::data.table(
        id = seq_len(
            nrow(go_table)
        ),
        Term = stringr::str_wrap(
            go_table$Term,
            width = 46
        )
    )

    df[
        ,
        fill := GOCHORD_GO_COLORS[
            id
        ]
    ]

    n_terms <- nrow(df)

    df[
        ,
        y := seq(
            0.90,
            0.07,
            length.out = n_terms
        )
    ]

    ggplot() +
        annotate(
            "text",
            x = 0.00,
            y = 0.98,
            label = "GO Terms",
            hjust = 0,
            vjust = 1,
            fontface = "bold",
            size = 5.2 * (1 / 2.845276)
        ) +
        geom_rect(
            data = df,
            aes(
                xmin = 0.00,
                xmax = 0.050,
                ymin = y - 0.038,
                ymax = y + 0.038,
                fill = fill
            ),
            color = "#333333",
            linewidth = 0.15
        ) +
        geom_text(
            data = df,
            aes(
                x = 0.060,
                y = y,
                label = Term
            ),
            hjust = 0,
            vjust = 0.5,
            size = 3.95 * (1 / 2.845276),
            lineheight = 0.96,
            color = "#111111"
        ) +
        scale_fill_identity() +
        coord_cartesian(
            xlim = c(
                0,
                1
            ),
            ylim = c(
                0,
                1
            ),
            clip = "off"
        ) +
        theme_void() +
        theme(
            plot.margin = margin(
                t = 2,
                r = 0,
                b = 2,
                l = 0,
                unit = "pt"
            )
        )
}

make_abs_lfc_bar <- function() {

    p_data <- data.frame(
        x = seq(
            lfc_limits[1],
            lfc_limits[2],
            length.out = 200
        ),
        y = 1
    )

    ggplot(
        p_data,
        aes(
            x = x,
            y = y,
            fill = x
        )
    ) +
        geom_raster(
            height = 0.40
        ) +
        scale_fill_gradientn(
            colours = GOCHORD_LFC_COLORS,
            limits = lfc_limits
        ) +
        annotate(
            "text",
            x = mean(lfc_limits),
            y = 1.75,
            label = "Gene log2FC",
            size = 5.0 * (1 / 2.845276),
            fontface = "bold"
        ) +
        annotate(
            "text",
            x = lfc_breaks[1],
            y = 0.35,
            label = sprintf(
                "%.1f",
                lfc_breaks[1]
            ),
            size = 4.6 * (1 / 2.845276),
            hjust = 0
        ) +
        annotate(
            "text",
            x = lfc_breaks[2],
            y = 0.35,
            label = sprintf(
                "%.1f",
                lfc_breaks[2]
            ),
            size = 4.6 * (1 / 2.845276),
            hjust = 0.5
        ) +
        annotate(
            "text",
            x = lfc_breaks[3],
            y = 0.35,
            label = sprintf(
                "%.1f",
                lfc_breaks[3]
            ),
            size = 4.6 * (1 / 2.845276),
            hjust = 1
        ) +
        coord_cartesian(
            xlim = c(
                lfc_limits[1] -
                    abs(
                        diff(lfc_limits)
                    ) * 0.05,
                lfc_limits[2] +
                    abs(
                        diff(lfc_limits)
                    ) * 0.05
            ),
            ylim = c(
                0.10,
                2.05
            ),
            clip = "off"
        ) +
        theme_void() +
        theme(
            legend.position = "none",
            plot.margin = margin(
                t = 1,
                r = 0,
                b = 1,
                l = 0
            )
        )
}

build_gochord_panel <- function(
        go_data,
        panel_title,
        contrast_text,
        show_scale = TRUE
) {

    chord_plot <- make_chord_plot(
        go_data$chord,
        nrow(
            go_data$go_table
        ),
        panel_title
    )

    term_legend <- make_term_legend_compact(
        go_data$go_table
    )

    # Keep only one concise directional line in the figure.
    # Full DEG thresholds belong in the figure legend / Methods.
    subtitle_band <- cowplot::ggdraw() +
        cowplot::draw_label(
            contrast_text,
            x = 0.50,
            y = 0.50,
            hjust = 0.5,
            vjust = 0.5,
            size = 4.45,
            color = MUTED,
            fontfamily = FONT_FAMILY
        )

    chord_legend_row <- cowplot::plot_grid(
        chord_plot,
        term_legend,
        ncol = 2,
        rel_widths = c(
            0.50,
            0.50
        ),
        align = "h",
        axis = "tb"
    )

    if (isTRUE(show_scale)) {

        scale_row <- cowplot::plot_grid(
            NULL,
            make_abs_lfc_bar(),
            NULL,
            ncol = 3,
            rel_widths = c(
                0.27,
                0.46,
                0.27
            )
        )

        cowplot::plot_grid(
            subtitle_band,
            chord_legend_row,
            scale_row,
            ncol = 1,
            rel_heights = c(
                0.065,
                0.855,
                0.080
            )
        )

    } else {

        cowplot::plot_grid(
            subtitle_band,
            chord_legend_row,
            ncol = 1,
            rel_heights = c(
                0.07,
                0.93
            )
        )
    }
}

panel_E_base <- build_gochord_panel(
    go_expr_data,
    "DNMT3A mRNA",
    "High vs Low | positive direction = High",
    show_scale = FALSE
)

panel_F_base <- build_gochord_panel(
    go_meth_data,
    "cg21629895 methylation",
    "Low vs High | positive direction = Low",
    show_scale = TRUE
)

panel_E <- cowplot::ggdraw(
    panel_E_base
) +
    cowplot::draw_label(
        "E",
        x = 0.005,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )

panel_F <- cowplot::ggdraw(
    panel_F_base
) +
    cowplot::draw_label(
        "F",
        x = 0.005,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )



# -----------------------------------------------------------------------------
# 12. Panel F: Continuous Hallmark associations
# -----------------------------------------------------------------------------
# This panel follows the continuous-correlation model in the uploaded
# reference Figure 5 code, adapted to the present TCGA-LGG study.
#
# The RNA input is already log2(norm_count + 1), so no additional log
# transformation is performed. GSVA is applied to the normalized continuous
# expression matrix using a Gaussian kernel.

hallmark_gene_sets <- split(
    hallmark_term2gene$gene,
    hallmark_term2gene$term
)

hallmark_gene_sets <- lapply(
    hallmark_gene_sets,
    unique
)

run_hallmark_gsva <- function(
        expression_matrix,
        gene_sets
) {

    gsva_ns <- asNamespace("GSVA")

    if (exists(
        "gsvaParam",
        envir = gsva_ns,
        inherits = FALSE
    )) {

        param <- GSVA::gsvaParam(
            exprData = expression_matrix,
            geneSets = gene_sets,
            kcdf = "Gaussian",
            minSize = 10,
            maxSize = 500
        )

        score_matrix <- GSVA::gsva(
            param,
            verbose = FALSE
        )

    } else {

        score_matrix <- GSVA::gsva(
            expression_matrix,
            gene_sets,
            method = "gsva",
            kcdf = "Gaussian",
            min.sz = 10,
            max.sz = 500,
            verbose = FALSE
        )
    }

    as.matrix(
        score_matrix
    )
}

hallmark_gsva_scores <- run_hallmark_gsva(
    expr_matrix,
    hallmark_gene_sets
)

if (!identical(
    colnames(hallmark_gsva_scores),
    analysis_master$Sample_Core
)) {

    hallmark_gsva_scores <- hallmark_gsva_scores[
        ,
        analysis_master$Sample_Core,
        drop = FALSE
    ]
}

run_continuous_hallmark_cor <- function(
        score_matrix,
        continuous_variable,
        analysis_label
) {

    continuous_variable <- safe_numeric(
        continuous_variable
    )

    out_list <- lapply(
        seq_len(
            nrow(score_matrix)
        ),
        function(i) {

            pathway_score <- as.numeric(
                score_matrix[
                    i,
                    ,
                    drop = TRUE
                ]
            )

            keep <- is.finite(pathway_score) &
                is.finite(continuous_variable)

            if (sum(keep) < 20) {
                return(NULL)
            }

            test <- suppressWarnings(
                stats::cor.test(
                    pathway_score[keep],
                    continuous_variable[keep],
                    method = "spearman",
                    exact = FALSE
                )
            )

            data.frame(
                Pathway_ID = rownames(score_matrix)[i],
                Pathway = clean_hallmark_name(
                    rownames(score_matrix)[i]
                ),
                rho = unname(
                    test$estimate
                ),
                P_value = test$p.value,
                N = sum(keep),
                stringsAsFactors = FALSE
            )
        }
    )

    result <- data.table::rbindlist(
        out_list,
        fill = TRUE
    )

    result[
        ,
        FDR := stats::p.adjust(
            P_value,
            method = "BH"
        )
    ]

    result[
        ,
        Analysis := analysis_label
    ]

    result
}

continuous_rna <- run_continuous_hallmark_cor(
    hallmark_gsva_scores,
    analysis_master$DNMT3A_expression,
    "DNMT3A mRNA"
)

# Direction alignment for continuous methylation:
# Use the negative beta value so that positive Spearman rho corresponds to
# lower cg21629895 methylation being associated with higher Hallmark activity.
# This preserves the original continuous methylation information while aligning
# the displayed positive direction with the study's DNMT3A-high / cg21629895-low
# biological direction.
continuous_meth <- run_continuous_hallmark_cor(
    hallmark_gsva_scores,
    -safe_numeric(analysis_master$cg21629895),
    "cg21629895 low methylation"
)

select_continuous_top <- function(
        result_table,
        n_show = 15L
) {

    x <- data.table::copy(
        result_table
    )

    x <- x[
        is.finite(FDR) &
        is.finite(rho)
    ]

    # data.table::setorder() accepts column names, not expressions.
    # Create an explicit absolute-correlation column for secondary ranking.
    x[
        ,
        Abs_rho := abs(rho)
    ]

    data.table::setorder(
        x,
        FDR,
        -Abs_rho
    )

    x <- head(
        x,
        n_show
    )

    x[
        ,
        Abs_rho := NULL
    ]

    x
}

# Match Panel D to the exact Hallmark pathways displayed in Panel B.
# Panel B is the categorical shared-positive Hallmark result.
# Panel D now evaluates those SAME pathways using continuous associations,
# rather than independently selecting a new top-15 pathway set.

panel_B_shared_pathways <- as.character(
    panel_C_data$Pathway
)

continuous_rna_top <- data.table::copy(
    continuous_rna[
        Pathway %in% panel_B_shared_pathways
    ]
)

continuous_meth_top <- data.table::copy(
    continuous_meth[
        Pathway %in% panel_B_shared_pathways
    ]
)

# Preserve the Panel B pathway order in both continuous subplots.
continuous_rna_top[
    ,
    Pathway := factor(
        Pathway,
        levels = panel_B_shared_pathways
    )
]

continuous_meth_top[
    ,
    Pathway := factor(
        Pathway,
        levels = panel_B_shared_pathways
    )
]

continuous_combined <- data.table::rbindlist(
    list(
        continuous_rna_top,
        continuous_meth_top
    ),
    fill = TRUE
)

continuous_combined[
    ,
    neglog10FDR := -log10(
        pmax(
            FDR,
            .Machine$double.xmin
        )
    )
]

continuous_color_max <- max(
    continuous_combined$neglog10FDR,
    na.rm = TRUE
)

continuous_color_max <- max(
    1,
    ceiling(
        continuous_color_max
    )
)

# Use one shared symmetric Spearman-rho axis for BOTH continuous subplots.
# This makes zero the visual center and allows direct comparison between
# DNMT3A expression and the low-methylation direction.
# Use the same slightly expanded x-axis in both continuous panels.
# The full signed rho range is retained so negative correlations are not hidden,
# while the labeled ticks follow the requested 0.0, 0.2, 0.4 presentation.
continuous_rho_limit <- max(
    0.44,
    max(
        abs(
            continuous_combined$rho
        ),
        na.rm = TRUE
    ) + 0.02
)

continuous_rho_limit <- ceiling(
    continuous_rho_limit * 100
) / 100

continuous_rho_breaks <- c(
    0.0,
    0.2,
    0.4
)

build_continuous_panel <- function(
        plot_data,
        title_text,
        show_y_labels = TRUE
) {

    plot_data <- data.table::copy(
        plot_data
    )

    # Compute the color-mapping variable inside each subplot.
    # This was missing in V9.0 and caused:
    # "object 'neglog10FDR' not found" during geom_point() rendering.
    plot_data[
        ,
        neglog10FDR := -log10(
            pmax(
                FDR,
                .Machine$double.xmin
            )
        )
    ]

    # Preserve the exact shared pathway ordering inherited from Panel B.
    plot_data[
        ,
        Pathway := factor(
            as.character(Pathway),
            levels = levels(plot_data$Pathway)
        )
    ]

    p <- ggplot(
        plot_data,
        aes(
            x = rho,
            y = Pathway,
            color = neglog10FDR
        )
    ) +
        geom_vline(
            xintercept = 0,
            color = "#666666",
            linewidth = 0.35,
            linetype = "dashed"
        ) +
        geom_point(
            size = 1.80,
            alpha = 0.95
        ) +
        scale_color_gradientn(
            colours = c(
                "#2C105C",
                "#711F81",
                "#B63679",
                "#ED6925",
                "#F6D746"
            ),
            limits = c(
                0,
                continuous_color_max
            ),
            breaks = scales::pretty_breaks(
                n = 4
            )(
                c(
                    0,
                    continuous_color_max
                )
            ),
            oob = scales::squish,
            name = expression(-log[10]("FDR"))
        ) +
        # Match the uploaded reference continuous-panel x-axis.
        scale_x_continuous(
            limits = c(
                -0.25,
                0.55
            ),
            breaks = seq(
                -0.2,
                0.4,
                0.2
            ),
            labels = function(x) sprintf("%.1f", x),
            expand = c(
                0,
                0
            )
        ) +
        labs(
            title = title_text,
            x = expression("Spearman " * rho),
            y = NULL
        ) +
        theme_classic(
            base_size = 6,
            base_family = FONT_FAMILY
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 6.0,
                hjust = 0.5,
                margin = margin(
                    b = 2
                ),
                color = INK
            ),
            axis.title.x = element_text(
                face = "bold",
                size = 5.2,
                margin = margin(
                    t = 2
                ),
                color = INK
            ),
            axis.text.x = element_text(
                size = 5.0,
                color = INK
            ),
            axis.text.y = element_text(
                size = 5.0,
                color = INK,
                lineheight = 0.85
            ),
            axis.line = element_line(
                linewidth = 0.35,
                color = INK
            ),
            axis.ticks = element_line(
                linewidth = 0.30,
                color = INK
            ),
            legend.position = "none",
            plot.margin = margin(
                t = 2,
                r = 3,
                b = 2,
                l = 2
            )
        )

    if (!isTRUE(show_y_labels)) {
        p <- p +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )
    }

    p
}


continuous_plot_rna <- build_continuous_panel(
    continuous_rna_top,
    "DNMT3A mRNA",
    show_y_labels = TRUE
)

continuous_plot_meth <- build_continuous_panel(
    continuous_meth_top,
    "cg21629895 low methylation",
    show_y_labels = TRUE
)

# Reference-style shared continuous color bar.
dummy_continuous_colorbar <- ggplot(
    data.frame(
        x = seq(
            0,
            continuous_color_max,
            length.out = 100
        ),
        y = 1,
        z = seq(
            0,
            continuous_color_max,
            length.out = 100
        )
    ),
    aes(
        x = x,
        y = y,
        color = z
    )
) +
    geom_point() +
    scale_color_gradientn(
        colours = c(
            "#2C105C",
            "#711F81",
            "#B63679",
            "#ED6925",
            "#F6D746"
        ),
        limits = c(
            0,
            continuous_color_max
        ),
        breaks = scales::pretty_breaks(
            n = 4
        )(
            c(
                0,
                continuous_color_max
            )
        ),
        name = expression(-log[10]("FDR"))
    ) +
    theme_void() +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(
            size = 5.0,
            face = "bold",
            vjust = 0.8
        ),
        legend.text = element_text(
            size = 5.0
        ),
        legend.key.height = unit(
            1.4,
            "mm"
        ),
        legend.key.width = unit(
            10,
            "mm"
        ),
        legend.margin = margin(
            t = 0,
            b = 0
        )
    )

continuous_colorbar_grob <- cowplot::get_legend(
    dummy_continuous_colorbar
)

continuous_plots_row <- cowplot::plot_grid(
    continuous_plot_rna,
    continuous_plot_meth,
    ncol = 2,
    rel_widths = c(
        1,
        1
    ),
    align = "hv",
    axis = "tblr"
)

continuous_colorbar_row <- cowplot::ggdraw() +
    cowplot::draw_grob(
        continuous_colorbar_grob,
        x = 0.50,
        y = 0.50,
        hjust = 0.50,
        vjust = 0.50
    )

continuous_core <- cowplot::plot_grid(
    continuous_plots_row,
    continuous_colorbar_row,
    ncol = 1,
    rel_heights = c(
        0.92,
        0.08
    )
)

panel_D_final <- cowplot::ggdraw(
    continuous_core
) +
    cowplot::draw_label(
        "D",
        x = 0.006,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )


# -----------------------------------------------------------------------------
# 13. Assemble final integrated figure
# -----------------------------------------------------------------------------

# GOChord subplots are intentionally NOT given separate manuscript panel
# letters. They are two subplots within one GO biological-process panel (E),
# analogous to the two cohort subplots in the uploaded reference GO panel.

go_expr_subpanel <- panel_E_base
go_meth_subpanel <- panel_F_base

go_pair_vertical <- cowplot::plot_grid(
    go_expr_subpanel,
    go_meth_subpanel,
    ncol = 1,
    rel_heights = c(
        1,
        1
    ),
    align = "v",
    axis = "lr"
)

panel_E_final <- cowplot::ggdraw(
    go_pair_vertical
) +
    cowplot::draw_label(
        "E",
        x = 0.006,
        y = 0.995,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    )

# B and C are matched half-width panels with one shared centered legend.
bc_plot_row <- cowplot::plot_grid(
    panel_B_final,
    panel_C_final,
    ncol = 2,
    rel_widths = c(
        1,
        1
    ),
    align = "hv",
    axis = "tblr"
)

row_BC <- cowplot::plot_grid(
    bc_plot_row,
    shared_bc_legend,
    ncol = 1,
    rel_heights = c(
        0.90,
        0.10
    )
)

# Bottom row now matches the organizational model of the uploaded reference:
#
# LEFT (D)
#   DNMT3A continuous | cg21629895 continuous
#
# RIGHT (E)
#   DNMT3A GOChord
#   cg21629895 GOChord
#
# Fine-tune the bottom-row horizontal geometry.
# Panel E is shifted slightly to the right so that its outer right edge aligns
# visually with the right edge of Panel C. This also opens the D-E gutter
# without changing the internal geometry of Panels D or E.
panel_E_bottom_aligned <- cowplot::ggdraw() +
    cowplot::draw_plot(
        panel_E_final,
        x = 0.075,
        y = 0,
        width = 0.925,
        height = 1
    )

bottom_row <- cowplot::plot_grid(
    panel_D_final,
    panel_E_bottom_aligned,
    ncol = 2,
    rel_widths = c(
        0.50,
        0.50
    ),
    align = "hv",
    axis = "tblr"
)

final_integrated_figure <- cowplot::plot_grid(
    panel_A,
    row_BC,
    bottom_row,
    ncol = 1,
    rel_heights = c(
        0.84,
        0.92,
        1.18
    ),
    align = "v",
    axis = "lr"
)


# -----------------------------------------------------------------------------
# 14. Save helpers
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# 15. Save final figure and individual panels
# -----------------------------------------------------------------------------

combined_pdf <- file.path(
    OUTPUT_DIR,
    "DNMT3A_Integrated_Figure_FINAL_V9_3_FINAL.pdf"
)

combined_png <- file.path(
    OUTPUT_DIR,
    "DNMT3A_Integrated_Figure_FINAL_V9_3_FINAL_600dpi.png"
)

save_pdf(
    final_integrated_figure,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

save_png(
    final_integrated_figure,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

save_pdf(
    panel_A,
    file.path(
        OUTPUT_DIR,
        "Integrated_Panel_A_Concordant_DEG_Heatmap.pdf"
    ),
    180,
    72
)

save_pdf(
    panel_B_final,
    file.path(
        OUTPUT_DIR,
        "Integrated_Panel_B_Positive_Hallmark.pdf"
    ),
    88,
    72
)

save_pdf(
    panel_C_final,
    file.path(
        OUTPUT_DIR,
        "Integrated_Panel_C_Positive_VIPER.pdf"
    ),
    88,
    72
)

save_pdf(
    panel_D_final,
    file.path(
        OUTPUT_DIR,
        "Integrated_Panel_D_Continuous_Hallmark_Associations.pdf"
    ),
    80,
    108
)

save_pdf(
    panel_E_final,
    file.path(
        OUTPUT_DIR,
        "Integrated_Panel_E_GO_BP_Programs.pdf"
    ),
    94,
    108
)


# -----------------------------------------------------------------------------
# 16. Save analysis tables
# -----------------------------------------------------------------------------

data.table::fwrite(
    gsea_integrated,
    file.path(
        OUTPUT_DIR,
        "Integrated_Hallmark_GSEA_ALL.csv"
    )
)

data.table::fwrite(
    panel_C_data,
    file.path(
        OUTPUT_DIR,
        "Integrated_Hallmark_Positive_Shared_DISPLAYED.csv"
    )
)

data.table::fwrite(
    panel_D_data,
    file.path(
        OUTPUT_DIR,
        "Integrated_VIPER_Positive_Shared_DISPLAYED.csv"
    )
)

data.table::fwrite(
    continuous_rna,
    file.path(
        OUTPUT_DIR,
        "Integrated_Continuous_Hallmark_DNMT3A_expression_ALL.csv"
    )
)

data.table::fwrite(
    continuous_meth,
    file.path(
        OUTPUT_DIR,
        "Integrated_Continuous_Hallmark_cg21629895_LOW_DIRECTION_ALL.csv"
    )
)

data.table::fwrite(
    continuous_rna_top,
    file.path(
        OUTPUT_DIR,
        "Integrated_Continuous_Hallmark_DNMT3A_expression_DISPLAYED.csv"
    )
)

data.table::fwrite(
    continuous_meth_top,
    file.path(
        OUTPUT_DIR,
        "Integrated_Continuous_Hallmark_cg21629895_LOW_DIRECTION_DISPLAYED.csv"
    )
)


# -----------------------------------------------------------------------------
# 17. Console report
# -----------------------------------------------------------------------------

cat(
    "\n============================================================\n",
    "FINAL V9.3 integrated figure completed.\n",
    "Final canvas: 180 mm x 248 mm.\n",
    "\nLayout:\n",
    "  A       = full-width concordant DEG heatmap\n",
    "  B | C   = matched Hallmark and VIPER panels\n",
    "  D | E   = continuous Hallmark panel beside combined GOChord panel\n",
    "\nPanel D:\n",
    "  Left  = continuous DNMT3A mRNA Hallmark associations\n",
    "  Right = continuous cg21629895 LOW-methylation-direction associations for the SAME Hallmark pathways shown in B\n",
    "  These are arranged like the TCGA / BeatAML pair in the reference figure.\n",
    "\nPanel E:\n",
    "  Upper = DNMT3A-high positive-direction GOChord\n",
    "  Lower = cg21629895-low positive-direction GOChord\n",
    "  Both are subplots within ONE GO biological-process panel.\n",
    "\nContinuous analysis:\n",
    "  GSVA input = log2(norm_count + 1), no additional log transform\n",
    "  Association = Spearman rho\n",
    "  Multiple testing = BH FDR separately for each continuous variable\n",
    "\nImmune/TME analysis is intentionally excluded from this figure.\n",
    "Saved to:\n",
    normalizePath(
        OUTPUT_DIR,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n============================================================\n",
    sep = ""
)
