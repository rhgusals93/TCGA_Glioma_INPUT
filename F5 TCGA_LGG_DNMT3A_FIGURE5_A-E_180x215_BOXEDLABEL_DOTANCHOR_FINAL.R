# ==============================================================================
# TCGA-LGG DNMT3A FIGURE 5 — FINAL INTEGRATED A–E
#
# FINAL PANEL STRUCTURE
#   A. Shared concordant RPPA feature heatmap
#   B. RPPA differential effect-size concordance
#   C. Concordant RPPA paired-effect summary
#   D. DNMT3A-centered RPPA STRING core network
#   E. DNMT3A-centered VIPER TF STRING core network
#
# FINAL MERGED FIGURE
#   Width  = 180 mm
#   Height = 215 mm
#   600 dpi PNG + vector PDF
#
# INPUTS SELECTED MANUALLY
#   1. Completed DNMT3A master CSV
#   2. TCGA-LGG RPPA file
#   3-5. RPPA STRING interaction / coordinate / degree TSVs
#   6-8. VIPER STRING interaction / coordinate / degree TSVs
#
# Notes
#   - RPPA A-C calculations are unchanged from the supplied RPPA script.
#   - STRING core definition is unchanged from the supplied STRING script:
#       STRING score >= 0.400 and shortest-path distance <= 2 from DNMT3A.
#   - D1 and D2 are relabeled D and E in the merged manuscript figure.
#   - Associations are descriptive/associational, not causal.
# ==============================================================================


# ==============================================================================
# TCGA-LGG DNMT3A RPPA Figure
# FINAL CLEAN MAIN FIGURE - Panels A-B only
#
# Data source:
#   UCSC Xena legacy TCGA-LGG RPPA
#   Unit: normalized RPPA value
#
# Final panel structure:
#   A. Shared concordant RPPA feature heatmap
#      - pheatmap-style sample-level heatmap
#      - top annotations:
#          DNMT3A expression
#          cg21629895 methylation
#          Molecular subtype
#
#   B. RPPA differential effect-size concordance
#      - DNMT3A High vs Low
#      - cg21629895 Low vs High
#
# Supplementary:
#   Molecular-subtype-adjusted RPPA associations
#   - standardized RPPA ~ standardized biomarker + Molecular_3Group
#
# IMPORTANT:
# - Continuous dual-association panel has been removed.
# - STRING PPI will be prepared separately as a later panel.
# - No blank space is reserved for STRING PPI in this figure.
# - Protein-type annotation is simplified to:
#       "Phosphorylation"
#       "General"
# - High/Low groups are taken directly from the completed DNMT3A master file.
# - No median cutoff is recalculated.
# - Only primary tumor RPPA samples ending in "-01" are retained.
# - Associations are descriptive/associational, not causal.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Package setup
# ------------------------------------------------------------------------------

cran_packages <- c(
    "data.table",
    "ggplot2",
    "cowplot",
    "scales",
    "ggrepel",
    "sandwich",
    "lmtest",
    "pheatmap",
    "RColorBrewer",
    "openxlsx"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(
            pkg,
            dependencies = TRUE
        )
    }
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(scales)
    library(pheatmap)
    library(openxlsx)
})


# ------------------------------------------------------------------------------
# 1. Global publication settings
# ------------------------------------------------------------------------------

set.seed(1234)

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 215
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"
FONT_PANEL <- 6.8
FONT_TITLE <- 6.2
FONT_AXIS_TITLE <- 5.6
FONT_AXIS_TEXT <- 5.2
FONT_LEGEND <- 5.0
FONT_LABEL <- 4.4
FONT_HEATMAP <- 5.0
FONT_CELL <- 5.0

INK <- "#222222"
MUTED <- "#666666"

MRNA_COLOR <- "#D1495B"
METH_COLOR <- "#3B82C4"

TYPE_COLORS <- c(
    "Phosphorylation" = "#7B61A8",
    "General" = "#707070"
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

cat("\n1. Select the completed DNMT3A master CSV file.\n")
file_1_path <- file.choose()

cat("\n2. Select the TCGA-LGG RPPA file from UCSC Xena.\n")
file_2_path <- file.choose()


read_any_table <- function(path) {

    lower_path <- tolower(path)

    if (grepl("\\.csv$", lower_path)) {
        return(
            data.table::fread(
                path,
                sep = ",",
                check.names = FALSE
            )
        )
    }

    data.table::fread(
        path,
        sep = "\t",
        check.names = FALSE
    )
}

file_1 <- read_any_table(file_1_path)
file_2 <- read_any_table(file_2_path)

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


# ------------------------------------------------------------------------------
# 3. Identify master and RPPA files automatically
# ------------------------------------------------------------------------------

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
        master_signature %in%
            names(x)
    )
}

is_rppa_like <- function(x) {

    sample_like_columns <- names(x)[
        grepl(
            "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
            names(x)
        )
    ]

    length(sample_like_columns) > 100
}

if (
    is_master_like(file_1) &&
    is_rppa_like(file_2)
) {

    master <- file_1
    rppa_raw <- file_2

} else if (
    is_rppa_like(file_1) &&
    is_master_like(file_2)
) {

    cat("\nThe selected files were reversed and were reassigned automatically.\n")

    master <- file_2
    rppa_raw <- file_1

} else {

    stop(
        paste0(
            "The selected files could not be identified as the required master/RPPA pair.\n",
            "The master file must contain the DNMT3A/cg21629895 grouping columns and Molecular_3Group."
        )
    )
}


# ------------------------------------------------------------------------------
# 4. Output directory
# ------------------------------------------------------------------------------

choose_output_dir <- function() {

    selected <- NULL

    if (.Platform$OS.type == "windows") {

        selected <- tryCatch(
            utils::choose.dir(
                default = Sys.getenv("USERPROFILE"),
                caption = "Select a SHORT output folder for the DNMT3A RPPA figure"
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
            "DNMT3A_RPPA_Figure_Revised"
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


# ------------------------------------------------------------------------------
# 5. Prepare master variables
# ------------------------------------------------------------------------------

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
        levels = c(
            "Low",
            "High"
        )
    )
]

master[
    ,
    cg21629895_Group := factor(
        clean_binary_group(
            cg21629895_Median_Group,
            "cg21629895_Median_Group"
        ),
        levels = c(
            "High",
            "Low"
        )
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


# ------------------------------------------------------------------------------
# 6. Prepare RPPA matrix
# ------------------------------------------------------------------------------

feature_column <- names(rppa_raw)[1]

data.table::setnames(
    rppa_raw,
    feature_column,
    "RPPA_Feature"
)

rppa_raw[
    ,
    RPPA_Feature := as.character(
        RPPA_Feature
    )
]

sample_columns <- names(rppa_raw)[
    grepl(
        "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
        names(rppa_raw)
    )
]

primary_rppa_samples <- sample_columns[
    substr(
        sample_columns,
        14,
        15
    ) == "01"
]

overlap_samples <- intersect(
    master$Sample_Core,
    primary_rppa_samples
)

if (length(overlap_samples) < 20) {
    stop(
        paste0(
            "Too few overlapping primary tumor RPPA samples were detected: ",
            length(overlap_samples)
        )
    )
}

rppa_dt <- rppa_raw[
    !is.na(RPPA_Feature) &
        RPPA_Feature != "",
    c(
        "RPPA_Feature",
        overlap_samples
    ),
    with = FALSE
]

if (
    anyDuplicated(
        rppa_dt$RPPA_Feature
    ) > 0
) {

    cat(
        "\nDuplicate RPPA feature names were detected and will be collapsed by mean.\n"
    )

    rppa_long_tmp <- data.table::melt(
        rppa_dt,
        id.vars = "RPPA_Feature",
        variable.name = "Sample_Core",
        value.name = "RPPA_Value"
    )

    rppa_long_tmp[
        ,
        RPPA_Value := safe_numeric(
            RPPA_Value
        )
    ]

    rppa_collapsed <- rppa_long_tmp[
        is.finite(RPPA_Value),
        .(
            RPPA_Value = mean(
                RPPA_Value,
                na.rm = TRUE
            )
        ),
        by = .(
            RPPA_Feature,
            Sample_Core
        )
    ]

    rppa_wide <- data.table::dcast(
        rppa_collapsed,
        RPPA_Feature ~ Sample_Core,
        value.var = "RPPA_Value"
    )

} else {

    rppa_wide <- copy(
        rppa_dt
    )
}

feature_names <- rppa_wide$RPPA_Feature

rppa_matrix <- as.matrix(
    rppa_wide[
        ,
        -1,
        with = FALSE
    ]
)

rownames(
    rppa_matrix
) <- feature_names

storage.mode(
    rppa_matrix
) <- "numeric"

sample_order <- unique(
    master$Sample_Core[
        master$Sample_Core %in%
            colnames(rppa_matrix)
    ]
)

rppa_matrix <- rppa_matrix[
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

if (
    !identical(
        analysis_master$Sample_Core,
        colnames(rppa_matrix)
    )
) {
    stop(
        "Master/RPPA sample alignment failed."
    )
}

cat("\nAligned RPPA matrix:\n")
cat("Features:", nrow(rppa_matrix), "\n")
cat("Samples:", ncol(rppa_matrix), "\n")
cat("Unit: normalized RPPA value\n")


# ------------------------------------------------------------------------------
# 7. Simplified protein-type annotation
# ------------------------------------------------------------------------------

# Display annotation is intentionally concise:
#   Phosphorylation = phosphosite-specific antibody/readout
#   General         = all other RPPA protein readouts

feature_type <- data.table::data.table(
    RPPA_Feature = rownames(rppa_matrix)
)

feature_type[
    ,
    Protein_Type := ifelse(
        grepl(
            "_p[STY][0-9]+",
            RPPA_Feature,
            ignore.case = FALSE
        ),
        "Phosphorylation",
        "General"
    )
]

feature_type[
    ,
    Protein_Type := factor(
        Protein_Type,
        levels = c(
            "General",
            "Phosphorylation"
        )
    )
]


# Clean display labels:
# Remove terminal RPPA antibody/source suffixes such as -R-V, -M-C, -R-E,
# -G-C, etc. Original RPPA_Feature names remain unchanged internally.
feature_type[
    ,
    Display_Feature := sub(
        "-[A-Za-z]-[A-Za-z]$",
        "",
        RPPA_Feature
    )
]


# ------------------------------------------------------------------------------
# 8. Hedges' g helper
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


# ------------------------------------------------------------------------------
# 9. Differential RPPA statistics
# ------------------------------------------------------------------------------

compare_rppa <- function(
        comparison_name,
        group_vector,
        group1,
        group0
) {

    output <- vector(
        "list",
        nrow(rppa_matrix)
    )

    for (i in seq_len(nrow(rppa_matrix))) {

        values <- safe_numeric(
            rppa_matrix[
                i,
                ,
                drop = TRUE
            ]
        )

        x1 <- values[
            group_vector == group1
        ]

        x0 <- values[
            group_vector == group0
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

        output[[i]] <- data.table::data.table(
            RPPA_Feature = rownames(rppa_matrix)[i],
            Comparison = comparison_name,
            Group1 = group1,
            Group0 = group0,
            N_Group1 = sum(is.finite(x1)),
            N_Group0 = sum(is.finite(x0)),
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

    result <- data.table::rbindlist(
        output,
        fill = TRUE
    )

    result[
        ,
        FDR := stats::p.adjust(
            P_Value,
            method = "BH"
        )
    ]

    result
}

stats_expression <- compare_rppa(
    "DNMT3A High vs Low",
    analysis_master$DNMT3A_Group,
    "High",
    "Low"
)

stats_methylation <- compare_rppa(
    "cg21629895 Low vs High",
    analysis_master$cg21629895_Group,
    "Low",
    "High"
)

differential_wide <- merge(
    stats_expression[
        ,
        .(
            RPPA_Feature,
            Hedges_g_Expression = Hedges_g,
            FDR_Expression = FDR
        )
    ],
    stats_methylation[
        ,
        .(
            RPPA_Feature,
            Hedges_g_Methylation = Hedges_g,
            FDR_Methylation = FDR
        )
    ],
    by = "RPPA_Feature",
    all = TRUE,
    sort = FALSE
)

differential_wide <- merge(
    differential_wide,
    feature_type,
    by = "RPPA_Feature",
    all.x = TRUE,
    sort = FALSE
)

differential_wide[
    ,
    Both_Significant :=
        is.finite(FDR_Expression) &
        is.finite(FDR_Methylation) &
        FDR_Expression < 0.05 &
        FDR_Methylation < 0.05
]

differential_wide[
    ,
    Concordance_Class := data.table::fcase(
        Both_Significant &
            Hedges_g_Expression > 0 &
            Hedges_g_Methylation > 0,
        "Concordant positive",

        Both_Significant &
            Hedges_g_Expression < 0 &
            Hedges_g_Methylation < 0,
        "Concordant negative",

        Both_Significant &
            sign(Hedges_g_Expression) !=
            sign(Hedges_g_Methylation),
        "Discordant significant",

        default = "Partial / not significant"
    )
]

differential_wide[
    ,
    Shared_Strength := pmin(
        abs(Hedges_g_Expression),
        abs(Hedges_g_Methylation)
    )
]


# ------------------------------------------------------------------------------
# 10. Select shared concordant RPPA features for Panel A
# ------------------------------------------------------------------------------

shared_concordant <- differential_wide[
    Concordance_Class %in% c(
        "Concordant positive",
        "Concordant negative"
    )
]

if (nrow(shared_concordant) > 0) {

    shared_concordant[
        ,
        Joint_FDR := pmax(
            FDR_Expression,
            FDR_Methylation,
            na.rm = TRUE
        )
    ]

    shared_concordant[
        ,
        Rank_Score := -log10(
            pmax(
                Joint_FDR,
                1e-300
            )
        ) * Shared_Strength
    ]

    selected_A <- shared_concordant[
        order(
            -Rank_Score
        )
    ][
        seq_len(
            min(
                .N,
                30
            )
        )
    ]

} else {

    warning(
        paste0(
            "No RPPA features were concordantly significant at FDR < 0.05 in both comparisons. ",
            "Panel A will display the top 30 features ranked by shared absolute effect size."
        )
    )

    selected_A <- differential_wide[
        order(
            -Shared_Strength
        )
    ][
        seq_len(
            min(
                .N,
                30
            )
        )
    ]
}


# ------------------------------------------------------------------------------
# 11. Panel A: pheatmap-style shared concordant RPPA heatmap
# ------------------------------------------------------------------------------

mat_A <- rppa_matrix[
    selected_A$RPPA_Feature,
    ,
    drop = FALSE
]

z_matrix <- t(
    apply(
        mat_A,
        1,
        function(x) {

            x <- safe_numeric(x)

            if (
                all(!is.finite(x)) ||
                !is.finite(
                    stats::sd(
                        x,
                        na.rm = TRUE
                    )
                ) ||
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

original_heatmap_features <- rownames(mat_A)

display_heatmap_features <- feature_type[
    match(
        original_heatmap_features,
        RPPA_Feature
    ),
    Display_Feature
]

rownames(
    z_matrix
) <- display_heatmap_features

colnames(
    z_matrix
) <- colnames(mat_A)

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


# pheatmap displays annotation rows in reverse column order.
# Columns are therefore supplied in reverse to obtain the visual order:
#   DNMT3A expression
#   cg21629895 methylation
#   Molecular subtype

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

rownames(
    annotation_col
) <- ordered_samples




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
    fontsize_row = 5.0,
    fontsize_col = FONT_HEATMAP,
    treeheight_row = 14,
    treeheight_col = 0,
    main = paste0(
        "Shared concordant RPPA features (top ",
        nrow(z_matrix),
        ")"
    ),
    silent = TRUE
)

panel_A_grob <- grid::grobTree(
    panel_A_heatmap_object$gtable,
    grid::textGrob(
        "A",
        x = grid::unit(
            0.005,
            "npc"
        ),
        y = grid::unit(
            0.995,
            "npc"
        ),
        just = c(
            "left",
            "top"
        ),
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


# ------------------------------------------------------------------------------
# 12. Panel B: RPPA differential effect-size concordance
# ------------------------------------------------------------------------------

cor_B <- suppressWarnings(
    stats::cor.test(
        differential_wide$Hedges_g_Expression,
        differential_wide$Hedges_g_Methylation,
        method = "spearman",
        exact = FALSE
    )
)

rho_B <- unname(
    cor_B$estimate
)

p_B <- cor_B$p.value

# Publication-style selective labeling for Panel B.
# Apply the SAME rule to all three significant classes:
#   1. Both comparisons must be significant (BH FDR < 0.05).
#   2. Rank proteins by Shared_Strength.
#   3. Display the top 5 proteins within each class.
#
# Shared_Strength is defined identically for every class as:
#   min(abs(Hedges_g_Expression), abs(Hedges_g_Methylation))

TOP_LABELS_PER_CLASS <- 5

label_B_positive <- differential_wide[
    Both_Significant == TRUE &
        Concordance_Class == "Concordant positive" &
        is.finite(Hedges_g_Expression) &
        is.finite(Hedges_g_Methylation)
][
    order(-Shared_Strength)
][
    seq_len(min(.N, TOP_LABELS_PER_CLASS))
]

label_B_negative <- differential_wide[
    Both_Significant == TRUE &
        Concordance_Class == "Concordant negative" &
        is.finite(Hedges_g_Expression) &
        is.finite(Hedges_g_Methylation)
][
    order(-Shared_Strength)
][
    seq_len(min(.N, TOP_LABELS_PER_CLASS))
]

label_B_discordant <- differential_wide[
    Both_Significant == TRUE &
        Concordance_Class == "Discordant significant" &
        is.finite(Hedges_g_Expression) &
        is.finite(Hedges_g_Methylation)
][
    order(-Shared_Strength)
][
    seq_len(min(.N, TOP_LABELS_PER_CLASS))
]

label_B <- data.table::rbindlist(
    list(
        label_B_positive,
        label_B_negative,
        label_B_discordant
    ),
    fill = TRUE
)

if (!"Display_Feature" %in% names(label_B)) {
    stop(
        "Display_Feature is missing from Panel B labels. Check the feature_type merge in Section 9."
    )
}

limit_B <- max(
    abs(
        c(
            differential_wide$Hedges_g_Expression,
            differential_wide$Hedges_g_Methylation
        )
    ),
    na.rm = TRUE
)

if (
    !is.finite(limit_B) ||
    limit_B <= 0
) {
    limit_B <- 0.5
}

limit_B <- max(
    0.25,
    ceiling(
        limit_B * 10
    ) / 10
)

plot_limit_B <- limit_B * 1.12

panel_B_base <- ggplot(
    differential_wide,
    aes(
        x = Hedges_g_Expression,
        y = Hedges_g_Methylation
    )
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
    geom_point(
        data = differential_wide[
            Both_Significant == FALSE
        ],
        aes(
            color = Concordance_Class
        ),
        size = 0.66,
        alpha = 0.22,
        stroke = 0.20
    ) +
    geom_point(
        data = differential_wide[
            Both_Significant == TRUE
        ],
        aes(
            color = Concordance_Class
        ),
        size = 0.88,
        alpha = 0.88,
        stroke = 0.25
    ) +
    # Fig. 4-like boxed labels keep protein names above the point cloud.
    ggrepel::geom_label_repel(
        data = label_B_positive,
        aes(label = Display_Feature),
        size = pt_to_mm(4.7),
        family = FONT_FAMILY,
        color = CONCORDANCE_COLORS["Concordant positive"],
        fill = scales::alpha("white", 0.96),
        label.size = 0.13,
        label.r = grid::unit(0.45, "mm"),
        label.padding = grid::unit(0.42, "mm"),
        box.padding = 0.42,
        point.padding = 0.06,
        min.segment.length = 0,
        segment.size = 0.22,
        segment.color = scales::alpha(CONCORDANCE_COLORS["Concordant positive"], 0.78),
        arrow = grid::arrow(length = grid::unit(0.75, "mm"), type = "closed"),
        max.overlaps = Inf,
        force = 5.4,
        force_pull = 0.10,
        max.time = 45,
        max.iter = 100000,
        seed = 20260913,
        show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
        data = label_B_negative,
        aes(label = Display_Feature),
        size = pt_to_mm(4.7),
        family = FONT_FAMILY,
        color = CONCORDANCE_COLORS["Concordant negative"],
        fill = scales::alpha("white", 0.96),
        label.size = 0.13,
        label.r = grid::unit(0.45, "mm"),
        label.padding = grid::unit(0.42, "mm"),
        box.padding = 0.42,
        point.padding = 0.06,
        min.segment.length = 0,
        segment.size = 0.22,
        segment.color = scales::alpha(CONCORDANCE_COLORS["Concordant negative"], 0.78),
        arrow = grid::arrow(length = grid::unit(0.75, "mm"), type = "closed"),
        max.overlaps = Inf,
        force = 5.4,
        force_pull = 0.10,
        max.time = 45,
        max.iter = 100000,
        seed = 20260914,
        show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
        data = label_B_discordant,
        aes(label = Display_Feature),
        size = pt_to_mm(4.7),
        family = FONT_FAMILY,
        color = CONCORDANCE_COLORS["Discordant significant"],
        fill = scales::alpha("white", 0.96),
        label.size = 0.13,
        label.r = grid::unit(0.45, "mm"),
        label.padding = grid::unit(0.42, "mm"),
        box.padding = 0.42,
        point.padding = 0.06,
        min.segment.length = 0,
        segment.size = 0.22,
        segment.color = scales::alpha(CONCORDANCE_COLORS["Discordant significant"], 0.76),
        arrow = grid::arrow(length = grid::unit(0.75, "mm"), type = "closed"),
        max.overlaps = Inf,
        force = 5.4,
        force_pull = 0.10,
        max.time = 45,
        max.iter = 100000,
        seed = 20260915,
        show.legend = FALSE
    ) +
    # Exact-dot anchors for labeled proteins.
    # A small white halo separates the target dot from nearby text/segments,
    # while the colored center preserves the concordance class.
    geom_point(
        data = label_B,
        aes(
            x = Hedges_g_Expression,
            y = Hedges_g_Methylation
        ),
        shape = 21,
        size = 1.55,
        fill = "white",
        color = "white",
        stroke = 0.25,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    geom_point(
        data = label_B,
        aes(
            x = Hedges_g_Expression,
            y = Hedges_g_Methylation,
            fill = Concordance_Class
        ),
        shape = 21,
        size = 1.05,
        color = "#202020",
        stroke = 0.22,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    scale_fill_manual(
        values = CONCORDANCE_COLORS,
        guide = "none"
    ) +
    scale_color_manual(
        values = CONCORDANCE_COLORS,
        breaks = c(
            "Concordant positive",
            "Concordant negative",
            "Discordant significant",
            "Partial / not significant"
        ),
        name = NULL,
        guide = guide_legend(
            nrow = 2,
            byrow = TRUE
        )
    ) +
    scale_x_continuous(
        limits = c(
            -plot_limit_B,
            plot_limit_B
        ),
        breaks = scales::pretty_breaks(
            n = 5
        ),
        expand = expansion(
            mult = c(
                0.02,
                0.02
            )
        )
    ) +
    scale_y_continuous(
        limits = c(
            -plot_limit_B,
            plot_limit_B
        ),
        breaks = scales::pretty_breaks(
            n = 5
        ),
        expand = expansion(
            mult = c(
                0.02,
                0.02
            )
        )
    ) +
    labs(
        title = "RPPA effect-size concordance",
        subtitle = paste0(
            "Spearman rho = ",
            sprintf(
                "%.2f",
                rho_B
            ),
            "; P ",
            ifelse(
                p_B < 0.001,
                "< 0.001",
                paste0(
                    "= ",
                    sprintf(
                        "%.3f",
                        p_B
                    )
                )
            ),
            "; n = ",
            nrow(differential_wide)
        ),
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
            hjust = 0
        ),
        plot.subtitle = element_text(
            size = 5.0,
            color = MUTED,
            hjust = 0,
            margin = margin(
                b = 3
            )
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
        legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(
            size = 4.7
        ),
        legend.key.width = grid::unit(
            3.0,
            "mm"
        ),
        legend.spacing.x = grid::unit(
            1.2,
            "mm"
        ),
        plot.margin = margin(
            4,
            8,
            6,
            4
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
# 15. Panel C: complete concordant RPPA paired-effect summary + final layout
# ------------------------------------------------------------------------------

# Panel C displays ALL concordant significant proteins from Panel B.
# It uses one continuous y-axis so that the plotting region visually aligns
# with Panel B. A small gap separates positive and negative groups.
#
# Color encodes concordance direction:
#   red  = concordant positive
#   blue = concordant negative
#
# Shape encodes the comparison:
#   circle   = DNMT3A mRNA High vs Low
#   triangle = cg21629895 Low vs High
#
# Legends are removed from the panel and placed in a separate annotation row
# below the lower panels for a cleaner manuscript-style presentation.

summary_C <- differential_wide[
    Both_Significant == TRUE &
        Concordance_Class %in% c(
            "Concordant positive",
            "Concordant negative"
        ) &
        is.finite(Hedges_g_Expression) &
        is.finite(Hedges_g_Methylation)
]

summary_C[
    ,
    Direction_Rank := data.table::fcase(
        Concordance_Class == "Concordant positive", 1L,
        Concordance_Class == "Concordant negative", 2L,
        default = 9L
    )
]

summary_C[
    ,
    Mean_Abs_Effect := (
        abs(Hedges_g_Expression) +
        abs(Hedges_g_Methylation)
    ) / 2
]

data.table::setorder(
    summary_C,
    Direction_Rank,
    -Mean_Abs_Effect
)

# Assign compact numeric y positions with a small visual gap between groups.
n_positive_C <- sum(
    summary_C$Concordance_Class == "Concordant positive"
)

n_negative_C <- sum(
    summary_C$Concordance_Class == "Concordant negative"
)

gap_C <- 2.2

summary_C[
    Concordance_Class == "Concordant positive",
    Y_Pos := rev(
        seq_len(.N)
    ) + n_negative_C + gap_C
]

summary_C[
    Concordance_Class == "Concordant negative",
    Y_Pos := rev(
        seq_len(.N)
    )
]

summary_C_long <- data.table::melt(
    summary_C[
        ,
        .(
            RPPA_Feature,
            Display_Feature,
            Concordance_Class,
            Y_Pos,
            Hedges_g_Expression,
            Hedges_g_Methylation
        )
    ],
    id.vars = c(
        "RPPA_Feature",
        "Display_Feature",
        "Concordance_Class",
        "Y_Pos"
    ),
    measure.vars = c(
        "Hedges_g_Expression",
        "Hedges_g_Methylation"
    ),
    variable.name = "Effect_Type",
    value.name = "Hedges_g"
)

summary_C_long[
    ,
    Effect_Type := factor(
        Effect_Type,
        levels = c(
            "Hedges_g_Expression",
            "Hedges_g_Methylation"
        ),
        labels = c(
            "DNMT3A expression",
            "cg21629895 methylation"
        )
    )
]

# Protein labels are rendered manually as y-axis labels.
label_breaks_C <- summary_C$Y_Pos
label_values_C <- summary_C$Display_Feature

x_limit_C <- max(
    abs(
        c(
            summary_C$Hedges_g_Expression,
            summary_C$Hedges_g_Methylation
        )
    ),
    na.rm = TRUE
)

if (
    !is.finite(x_limit_C) ||
    x_limit_C <= 0
) {
    x_limit_C <- 0.5
}

x_limit_C <- ceiling(
    x_limit_C * 10
) / 10

# Panel C is split into two directional columns to improve readability while
# retaining ALL concordant significant proteins.

summary_C_positive <- summary_C[
    Concordance_Class == "Concordant positive"
][order(-Mean_Abs_Effect)]

summary_C_negative <- summary_C[
    Concordance_Class == "Concordant negative"
][order(-Mean_Abs_Effect)]

make_C_direction_plot <- function(dt, direction_name, direction_color) {

    dt <- data.table::copy(dt)
    dt[, C_Row := rev(seq_len(.N))]

    dt_long <- data.table::melt(
        dt[, .(
            RPPA_Feature,
            Display_Feature,
            Concordance_Class,
            C_Row,
            Hedges_g_Expression,
            Hedges_g_Methylation
        )],
        id.vars = c(
            "RPPA_Feature",
            "Display_Feature",
            "Concordance_Class",
            "C_Row"
        ),
        measure.vars = c(
            "Hedges_g_Expression",
            "Hedges_g_Methylation"
        ),
        variable.name = "Effect_Type",
        value.name = "Hedges_g"
    )

    dt_long[
        ,
        Effect_Type := factor(
            Effect_Type,
            levels = c(
                "Hedges_g_Expression",
                "Hedges_g_Methylation"
            ),
            labels = c(
                "DNMT3A expression",
                "cg21629895 methylation"
            )
        )
    ]

    ggplot(dt, aes(y = C_Row)) +
        geom_hline(
            yintercept = dt$C_Row,
            linewidth = 0.11,
            color = "#F3F3F3"
        ) +
        geom_segment(
            aes(
                x = Hedges_g_Expression,
                xend = Hedges_g_Methylation,
                yend = C_Row
            ),
            linewidth = 0.40,
            color = direction_color,
            alpha = 0.62
        ) +
        geom_point(
            data = dt_long,
            aes(
                x = Hedges_g,
                y = C_Row,
                shape = Effect_Type
            ),
            color = direction_color,
            size = 1.22,
            stroke = 0.22
        ) +
        geom_vline(
            xintercept = 0,
            linewidth = 0.28,
            linetype = "dashed",
            color = "#A8A8A8"
        ) +
        scale_shape_manual(
            values = c(
                "DNMT3A expression" = 16,
                "cg21629895 methylation" = 17
            ),
            guide = "none"
        ) +
        scale_y_continuous(
            breaks = dt$C_Row,
            labels = dt$Display_Feature,
            limits = c(0.4, max(dt$C_Row) + 0.6),
            expand = c(0, 0)
        ) +
        scale_x_continuous(
            limits = c(-x_limit_C, x_limit_C),
            breaks = scales::pretty_breaks(n = 4),
            expand = expansion(mult = c(0.04, 0.04))
        ) +
        labs(
            title = direction_name,
            x = "Hedges' g",
            y = NULL
        ) +
        theme_classic(
            base_size = 6,
            base_family = FONT_FAMILY
        ) +
        theme(
            plot.title = element_text(
                size = 5.2,
                face = "bold",
                color = direction_color,
                hjust = 0
            ),
            axis.title.x = element_text(
                size = 5.0,
                face = "bold",
                color = INK
            ),
            axis.text.x = element_text(
                size = 5.0,
                color = INK
            ),
            axis.text.y = element_text(
                size = 5.0,
                color = INK,
                margin = margin(r = 3.0)
            ),
            axis.line = element_line(
                linewidth = 0.32,
                color = INK
            ),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(
                linewidth = 0.28,
                color = INK
            ),
            legend.position = "none",
            plot.margin = margin(3, 4, 4, 4)
        )
}

panel_C_positive <- make_C_direction_plot(
    summary_C_positive,
    paste0("Positive (n = ", nrow(summary_C_positive), ")"),
    CONCORDANCE_COLORS["Concordant positive"]
)

panel_C_negative <- make_C_direction_plot(
    summary_C_negative,
    paste0("Negative (n = ", nrow(summary_C_negative), ")"),
    CONCORDANCE_COLORS["Concordant negative"]
)

panel_C_body <- cowplot::plot_grid(
    panel_C_positive,
    panel_C_negative,
    nrow = 1,
    rel_widths = c(0.88, 1.12),
    align = "hv",
    axis = "tb"
)

# Separate header row prevents the panel title from overlapping the
# Positive/Negative subsection titles.
panel_C_header <- cowplot::ggdraw() +
    cowplot::draw_label(
        "C",
        x = 0.002,
        y = 0.50,
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
        size = FONT_PANEL,
        color = INK
    ) +
    cowplot::draw_label(
        "Concordant RPPA effect-size summary",
        x = 0.055,
        y = 0.50,
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
        size = FONT_TITLE,
        color = INK
    )

panel_C <- cowplot::plot_grid(
    panel_C_header,
    panel_C_body,
    ncol = 1,
    rel_heights = c(
        0.075,
        0.925
    )
)




# ==============================================================================
# TCGA-LGG DNMT3A STRING NETWORK
# DNMT3A-CENTERED CORE NETWORK FOR MANUSCRIPT FIGURE
#
# Main design:
#   - DNMT3A is placed at the center.
#   - First-order DNMT3A neighbors are placed on the inner ring.
#   - Second-order DNMT3A neighbors are placed on the outer ring.
#   - Nodes beyond two steps from DNMT3A are excluded from the main figure.
#   - Core membership is defined using the same STRING threshold used for the
#     downloaded network (score >= 0.400).
#   - All STRING edges >= 0.400 among retained nodes are shown in the plot.
#
# This produces a compact, interpretable "core network" rather than a dense
# all-node STRING cloud.
#
# IMPORTANT:
#   The pruning rule is algorithmic and reproducible:
#       core = shortest-path distance <= 2 from DNMT3A
#       using STRING edges with score >= 0.400
#
# Defaults:
#   STRING plotting threshold = 0.400
#   Core-selection threshold = 0.400
#   Core radius = 2 steps
#
# Six STRING files are selected manually:
#   D1 RPPA: interaction, coordinates, node degree
#   D2 VIPER: interaction, coordinates, node degree
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

packages <- c(
    "data.table",
    "ggplot2",
    "ggrepel",
    "igraph",
    "cowplot",
    "scales"
)

missing_packages <- packages[
    !vapply(
        packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing_packages) > 0) {
    install.packages(
        missing_packages,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
    )
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(igraph)
    library(cowplot)
    library(scales)
})


# ------------------------------------------------------------------------------
# 1. Global settings
# ------------------------------------------------------------------------------

set.seed(20260914)

PLOT_SCORE <- 0.400
CORE_SELECTION_SCORE <- 0.400
CORE_ORDER <- 2

FIGURE_DPI <- 600

PANEL_WIDTH_MM <- 96
PANEL_HEIGHT_MM <- 90

COMBINED_WIDTH_MM <- 180
COMBINED_HEIGHT_MM <- 101

FONT_FAMILY <- "sans"

# Restrained Nature-style palette.
COLOR_CENTER <- "#9E3D3F"
COLOR_FIRST <- "#D48A7A"
COLOR_SECOND <- "#E8D8C7"

EDGE_DIRECT <- "#A95A5C"
EDGE_OTHER <- "#AEB4BA"

TEXT_COLOR <- "#202124"
MUTED_COLOR <- "#6B6F73"
NODE_BORDER <- "#FFFFFF"


# ------------------------------------------------------------------------------
# 2. Manual file selection
# ------------------------------------------------------------------------------

choose_input_file <- function(message_text) {

    cat("\n", message_text, "\n", sep = "")

    selected_file <- file.choose()

    if (!file.exists(selected_file)) {
        stop(
            paste0(
                "Selected file does not exist: ",
                selected_file
            )
        )
    }

    normalizePath(
        selected_file,
        winslash = "/",
        mustWork = TRUE
    )
}


cat("\n============================================================\n")
cat("SELECT SIX STRING EXPORT FILES\n")
cat("============================================================\n")

D1_INTERACTIONS_FILE <- choose_input_file(
    "1. Select D1 RPPA STRING short interaction TSV."
)

D1_COORDINATES_FILE <- choose_input_file(
    "2. Select D1 RPPA STRING network coordinate TSV."
)

D1_DEGREES_FILE <- choose_input_file(
    "3. Select D1 RPPA STRING node degree TSV."
)

D2_INTERACTIONS_FILE <- choose_input_file(
    "4. Select D2 VIPER STRING short interaction TSV."
)

D2_COORDINATES_FILE <- choose_input_file(
    "5. Select D2 VIPER STRING network coordinate TSV."
)

D2_DEGREES_FILE <- choose_input_file(
    "6. Select D2 VIPER STRING node degree TSV."
)


# ------------------------------------------------------------------------------
# 3. Output directory
# ------------------------------------------------------------------------------

# Use the SAME output folder already selected for the RPPA A-C analysis.
# This prevents duplicate nested folders and simplifies final Figure 5 export.

if (!exists("OUTPUT_DIR") || !dir.exists(OUTPUT_DIR)) {
    stop("The shared Figure 5 output directory is unavailable.")
}

cat("\nShared Figure 5 output directory:\n")
cat(OUTPUT_DIR, "\n")


# ------------------------------------------------------------------------------
# 4. Robust STRING readers
# ------------------------------------------------------------------------------

read_string_table <- function(path) {

    dt <- data.table::fread(
        path,
        sep = "\t",
        check.names = FALSE,
        data.table = TRUE
    )

    names(dt) <- trimws(
        sub(
            "^\ufeff",
            "",
            names(dt)
        )
    )

    dt
}


first_existing_column <- function(candidates, available_names) {

    hit <- candidates[
        candidates %in% available_names
    ]

    if (length(hit) == 0) {
        return(NA_character_)
    }

    hit[1]
}


standardize_interactions <- function(dt) {

    node1_col <- first_existing_column(
        c(
            "#node1",
            "node1",
            "preferredName_A",
            "preferredName1"
        ),
        names(dt)
    )

    node2_col <- first_existing_column(
        c(
            "node2",
            "preferredName_B",
            "preferredName2"
        ),
        names(dt)
    )

    score_col <- first_existing_column(
        c(
            "combined_score",
            "score",
            "combinedScore"
        ),
        names(dt)
    )

    if (
        is.na(node1_col) ||
        is.na(node2_col) ||
        is.na(score_col)
    ) {
        stop(
            paste0(
                "Interaction file is missing required columns.\n",
                "Detected columns: ",
                paste(names(dt), collapse = ", ")
            )
        )
    }

    out <- dt[
        ,
        .(
            node1 = as.character(get(node1_col)),
            node2 = as.character(get(node2_col)),
            combined_score = as.numeric(get(score_col))
        )
    ]

    if (
        nrow(out) > 0 &&
        max(out$combined_score, na.rm = TRUE) > 1.5
    ) {
        out[
            ,
            combined_score := combined_score / 1000
        ]
    }

    out <- out[
        is.finite(combined_score) &
            !is.na(node1) &
            !is.na(node2) &
            node1 != "" &
            node2 != ""
    ]

    unique(out)
}


standardize_nodes <- function(dt) {

    node_col <- first_existing_column(
        c(
            "#node",
            "node",
            "preferredName",
            "preferred_name",
            "protein"
        ),
        names(dt)
    )

    if (is.na(node_col)) {
        stop(
            paste0(
                "Node file is missing a recognizable node column.\n",
                "Detected columns: ",
                paste(names(dt), collapse = ", ")
            )
        )
    }

    unique(
        dt[
            ,
            .(
                node = as.character(get(node_col))
            )
        ][
            !is.na(node) &
                node != ""
        ]
    )
}


standardize_degrees <- function(dt) {

    node_col <- first_existing_column(
        c(
            "#node",
            "node",
            "preferredName",
            "preferred_name",
            "protein"
        ),
        names(dt)
    )

    degree_col <- first_existing_column(
        c(
            "node_degree",
            "degree",
            "NodeDegree"
        ),
        names(dt)
    )

    if (
        is.na(node_col) ||
        is.na(degree_col)
    ) {
        stop(
            paste0(
                "Node-degree file is missing required columns.\n",
                "Detected columns: ",
                paste(names(dt), collapse = ", ")
            )
        )
    }

    out <- dt[
        ,
        .(
            node = as.character(get(node_col)),
            node_degree = as.numeric(get(degree_col))
        )
    ]

    out[
        is.na(node_degree),
        node_degree := 0
    ]

    unique(
        out,
        by = "node"
    )
}


# ------------------------------------------------------------------------------
# 5. Read data
# ------------------------------------------------------------------------------

d1_edges_all <- standardize_interactions(
    read_string_table(
        D1_INTERACTIONS_FILE
    )
)

d1_nodes_all <- standardize_nodes(
    read_string_table(
        D1_COORDINATES_FILE
    )
)

d1_degrees_all <- standardize_degrees(
    read_string_table(
        D1_DEGREES_FILE
    )
)


d2_edges_all <- standardize_interactions(
    read_string_table(
        D2_INTERACTIONS_FILE
    )
)

d2_nodes_all <- standardize_nodes(
    read_string_table(
        D2_COORDINATES_FILE
    )
)

d2_degrees_all <- standardize_degrees(
    read_string_table(
        D2_DEGREES_FILE
    )
)


# ------------------------------------------------------------------------------
# 6. DNMT3A-centered core selection
# ------------------------------------------------------------------------------

build_dnmt3a_core <- function(
        edges_all,
        nodes_all,
        degrees_all,
        network_name
) {

    if (!"DNMT3A" %in% nodes_all$node) {
        stop(
            paste0(
                network_name,
                ": DNMT3A is not present in the STRING node list."
            )
        )
    }

    selection_edges <- edges_all[
        combined_score >= CORE_SELECTION_SCORE
    ]

    selection_graph <- igraph::graph_from_data_frame(
        d = as.data.frame(selection_edges),
        directed = FALSE,
        vertices = data.frame(
            name = nodes_all$node,
            stringsAsFactors = FALSE
        )
    )

    dist_to_dnmt3a <- igraph::distances(
        selection_graph,
        v = "DNMT3A",
        to = igraph::V(selection_graph),
        mode = "all"
    )

    dist_vector <- as.numeric(
        dist_to_dnmt3a[1, ]
    )

    names(dist_vector) <- colnames(
        dist_to_dnmt3a
    )

    keep_nodes <- names(
        dist_vector[
            is.finite(dist_vector) &
                dist_vector <= CORE_ORDER
        ]
    )

    if (length(keep_nodes) < 2) {
        stop(
            paste0(
                network_name,
                ": the high-confidence DNMT3A-centered core contains fewer than two nodes. ",
                "No DNMT3A-centered core was recovered at the selected STRING threshold of ",
                CORE_SELECTION_SCORE,
                "."
            )
        )
    }

    node_table <- data.table::data.table(
        node = keep_nodes,
        Core_Distance = as.numeric(
            dist_vector[keep_nodes]
        )
    )

    node_table <- merge(
        node_table,
        degrees_all,
        by = "node",
        all.x = TRUE,
        sort = FALSE
    )

    node_table[
        is.na(node_degree),
        node_degree := 0
    ]

    node_table[
        ,
        Shell := data.table::fcase(
            Core_Distance == 0,
            "DNMT3A",
            Core_Distance == 1,
            "First-order",
            Core_Distance == 2,
            "Second-order",
            default = "Other"
        )
    ]

    plot_edges <- edges_all[
        combined_score >= PLOT_SCORE &
            node1 %in% keep_nodes &
            node2 %in% keep_nodes
    ]

    list(
        nodes = node_table,
        edges = plot_edges
    )
}


D1_core <- build_dnmt3a_core(
    edges_all = d1_edges_all,
    nodes_all = d1_nodes_all,
    degrees_all = d1_degrees_all,
    network_name = "D1 RPPA"
)

D2_core <- build_dnmt3a_core(
    edges_all = d2_edges_all,
    nodes_all = d2_nodes_all,
    degrees_all = d2_degrees_all,
    network_name = "D2 VIPER"
)


# ------------------------------------------------------------------------------
# 7. Concentric DNMT3A-centered layout
# ------------------------------------------------------------------------------

make_core_layout <- function(core_object) {

    nodes <- data.table::copy(
        core_object$nodes
    )

    shell0 <- nodes[
        Core_Distance == 0
    ]

    shell1 <- nodes[
        Core_Distance == 1
    ][
        order(
            -node_degree,
            node
        )
    ]

    shell2 <- nodes[
        Core_Distance == 2
    ][
        order(
            -node_degree,
            node
        )
    ]


    shell0[
        ,
        `:=`(
            x = 0,
            y = 0
        )
    ]


    place_ring <- function(dt, radius, angle_offset = 0) {

        if (nrow(dt) == 0) {
            return(dt)
        }

        if (nrow(dt) == 1) {

            dt[
                ,
                `:=`(
                    x = radius,
                    y = 0
                )
            ]

            return(dt)
        }

        angles <- seq(
            0,
            2 * pi,
            length.out = nrow(dt) + 1
        )[
            seq_len(
                nrow(dt)
            )
        ] + angle_offset

        dt[
            ,
            `:=`(
                x = radius * cos(angles),
                y = radius * sin(angles)
            )
        ]

        dt
    }


    shell1 <- place_ring(
        shell1,
        radius = 1.38,
        angle_offset = 0
    )

    shell2 <- place_ring(
        shell2,
        radius = 2.62,
        angle_offset = pi / max(
            nrow(shell2),
            1
        )
    )


    layout <- data.table::rbindlist(
        list(
            shell0,
            shell1,
            shell2
        ),
        fill = TRUE
    )


    degree_values <- layout$node_degree

    if (
        length(unique(degree_values)) <= 1
    ) {

        layout[
            ,
            Node_Size := 9.5
        ]

    } else {

        layout[
            ,
            Node_Size :=
                scales::rescale(
                    sqrt(
                        node_degree + 1
                    ),
                    to = c(
                        8.0,
                        11.8
                    )
                )
        ]
    }


    layout[
        ,
        Node_Fill := data.table::fcase(
            Shell == "DNMT3A",
            COLOR_CENTER,
            Shell == "First-order",
            COLOR_FIRST,
            Shell == "Second-order",
            COLOR_SECOND,
            default = "#DDDDDD"
        )
    ]


    edges <- merge(
        core_object$edges,
        layout[
            ,
            .(
                node1 = node,
                x_from = x,
                y_from = y
            )
        ],
        by = "node1",
        all.x = TRUE,
        sort = FALSE
    )

    edges <- merge(
        edges,
        layout[
            ,
            .(
                node2 = node,
                x_to = x,
                y_to = y
            )
        ],
        by = "node2",
        all.x = TRUE,
        sort = FALSE
    )


    if (nrow(edges) > 0) {

        # Direct DNMT3A associations receive a subtle warm accent.
        # All other retained associations remain neutral gray.
        edges[
            ,
            Edge_Class := ifelse(
                node1 == "DNMT3A" | node2 == "DNMT3A",
                "Direct DNMT3A",
                "Other retained"
            )
        ]

        score_range <- range(
            edges$combined_score,
            na.rm = TRUE
        )

        if (
            !all(is.finite(score_range)) ||
            diff(score_range) == 0
        ) {
            edges[
                ,
                Edge_Width := 0.32
            ]

            edges[
                ,
                Edge_Alpha := 0.42
            ]
        } else {
            edges[
                ,
                Edge_Width :=
                    scales::rescale(
                        combined_score,
                        to = c(
                            0.12,
                            0.62
                        )
                    )
            ]

            edges[
                ,
                Edge_Alpha :=
                    scales::rescale(
                        combined_score,
                        to = c(
                            0.16,
                            0.58
                        )
                    )
            ]
        }
    }


    list(
        nodes = layout,
        edges = edges
    )
}


D1_layout <- make_core_layout(
    D1_core
)

D2_layout <- make_core_layout(
    D2_core
)


# ------------------------------------------------------------------------------
# 8. Core-network plotting function
# ------------------------------------------------------------------------------

make_core_plot <- function(
        network_object,
        panel_label,
        title_text
) {

    nodes <- data.table::copy(
        network_object$nodes
    )

    edges <- data.table::copy(
        network_object$edges
    )

    # Use slightly smaller text for long gene symbols while keeping all labels
    # fully inside the nodes.
    nodes[
        ,
        Label_Size := data.table::fcase(
            nchar(node) >= 7,
            3.6 / 2.845,
            nchar(node) >= 6,
            3.9 / 2.845,
            default = 4.2 / 2.845
        )
    ]

    p <- ggplot() +

        # Neutral retained associations first.
        geom_segment(
            data = edges[
                Edge_Class == "Other retained"
            ],
            aes(
                x = x_from,
                y = y_from,
                xend = x_to,
                yend = y_to,
                linewidth = Edge_Width,
                alpha = Edge_Alpha
            ),
            color = EDGE_OTHER,
            lineend = "round",
            show.legend = FALSE
        ) +

        # Direct DNMT3A associations are subtly emphasized.
        geom_segment(
            data = edges[
                Edge_Class == "Direct DNMT3A"
            ],
            aes(
                x = x_from,
                y = y_from,
                xend = x_to,
                yend = y_to,
                linewidth = Edge_Width,
                alpha = Edge_Alpha
            ),
            color = EDGE_DIRECT,
            lineend = "round",
            show.legend = FALSE
        ) +

        scale_linewidth_identity() +
        scale_alpha_identity() +

        # Second-order nodes.
        geom_point(
            data = nodes[
                Shell == "Second-order"
            ],
            aes(
                x = x,
                y = y,
                size = Node_Size,
                fill = I(Node_Fill)
            ),
            shape = 21,
            color = NODE_BORDER,
            stroke = 0.45,
            show.legend = FALSE
        ) +

        # First-order nodes.
        geom_point(
            data = nodes[
                Shell == "First-order"
            ],
            aes(
                x = x,
                y = y,
                size = Node_Size,
                fill = I(Node_Fill)
            ),
            shape = 21,
            color = NODE_BORDER,
            stroke = 0.55,
            show.legend = FALSE
        ) +

        # DNMT3A center.
        geom_point(
            data = nodes[
                Shell == "DNMT3A"
            ],
            aes(
                x = x,
                y = y,
                size = Node_Size,
                fill = I(Node_Fill)
            ),
            shape = 21,
            color = "#5A2021",
            stroke = 0.85,
            show.legend = FALSE
        ) +

        scale_size_identity() +

        # Gene symbols are centered inside the circles.
        geom_text(
            data = nodes[
                Shell != "DNMT3A"
            ],
            aes(
                x = x,
                y = y,
                label = node,
                size = Label_Size
            ),
            family = FONT_FAMILY,
            color = TEXT_COLOR,
            fontface = "plain",
            lineheight = 0.88,
            show.legend = FALSE
        ) +

        geom_text(
            data = nodes[
                Shell == "DNMT3A"
            ],
            aes(
                x = x,
                y = y,
                label = node
            ),
            family = FONT_FAMILY,
            color = "#FFFFFF",
            fontface = "bold",
            size = 4.25 / 2.845,
            lineheight = 0.88,
            show.legend = FALSE
        ) +

        coord_equal(
            xlim = c(
                -3.05,
                3.05
            ),
            ylim = c(
                -3.05,
                3.05
            ),
            clip = "off"
        ) +

        labs(
            title = title_text,
            subtitle = paste0(
                nrow(nodes),
                " nodes · ",
                nrow(edges),
                " STRING associations"
            )
        ) +

        theme_void(
            base_family = FONT_FAMILY
        ) +

        theme(
            plot.title = element_text(
                size = 7.4,
                face = "bold",
                color = TEXT_COLOR,
                hjust = 0,
                margin = margin(
                    b = 1.5
                )
            ),
            plot.subtitle = element_text(
                size = 5.0,
                color = MUTED_COLOR,
                hjust = 0,
                margin = margin(
                    b = 3
                )
            ),
            plot.margin = margin(
                7,
                7,
                6,
                7
            )
        )

    cowplot::ggdraw(
        p
    ) +
        cowplot::draw_label(
            panel_label,
            x = 0.004,
            y = 0.996,
            hjust = 0,
            vjust = 1,
            fontface = "bold",
            size = 8.2,
            color = TEXT_COLOR
        )
}



# ------------------------------------------------------------------------------
# 9. Shared annotation legend
# ------------------------------------------------------------------------------

make_network_legend <- function() {

    ggplot() +
        geom_point(aes(0.10, 0.5), shape = 21, size = 4.8,
                   fill = COLOR_CENTER, color = "#5A2021", stroke = 0.7) +
        annotate("text", x = 0.20, y = 0.5, label = "DNMT3A",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = TEXT_COLOR) +
        geom_point(aes(0.82, 0.5), shape = 21, size = 4.8,
                   fill = COLOR_FIRST, color = NODE_BORDER, stroke = 0.5) +
        annotate("text", x = 0.92, y = 0.5, label = "First-order neighbor",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = TEXT_COLOR) +
        geom_point(aes(1.96, 0.5), shape = 21, size = 4.8,
                   fill = COLOR_SECOND, color = NODE_BORDER, stroke = 0.45) +
        annotate("text", x = 2.06, y = 0.5, label = "Second-order neighbor",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = TEXT_COLOR) +
        annotate("segment", x = 3.30, xend = 3.54, y = 0.5, yend = 0.5,
                 color = EDGE_DIRECT, linewidth = 0.50, alpha = 0.72, lineend = "round") +
        annotate("text", x = 3.62, y = 0.5, label = "Direct DNMT3A association",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = TEXT_COLOR) +
        annotate("segment", x = 5.34, xend = 5.58, y = 0.5, yend = 0.5,
                 color = EDGE_OTHER, linewidth = 0.39, alpha = 0.55, lineend = "round") +
        annotate("text", x = 5.66, y = 0.5, label = "Other retained association",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = TEXT_COLOR) +
        coord_cartesian(xlim = c(-0.15, 7.05), ylim = c(0, 1), clip = "off") +
        theme_void(base_family = FONT_FAMILY) +
        theme(plot.margin = margin(0, 0, 0, 0))
}


# ------------------------------------------------------------------------------
# 10. Build D1 and D2 core panels
# ------------------------------------------------------------------------------

panel_D1 <- make_core_plot(
    network_object = D1_layout,
    panel_label = "D1",
    title_text = "DNMT3A-centered RPPA core network"
)

panel_D2 <- make_core_plot(
    network_object = D2_layout,
    panel_label = "D2",
    title_text = "DNMT3A-centered VIPER TF core network"
)


panel_row <- cowplot::plot_grid(
    panel_D1,
    panel_D2,
    nrow = 1,
    rel_widths = c(
        1,
        1
    ),
    align = "hv"
)

shared_legend <- make_network_legend()

figure_D <- cowplot::plot_grid(
    panel_row,
    shared_legend,
    ncol = 1,
    rel_heights = c(
        0.89,
        0.11
    ),
    align = "v"
)





# ==============================================================================
# FINAL INTEGRATED FIGURE 5 ASSEMBLY: A–E
# ==============================================================================

# ------------------------------------------------------------------------------
# Relabel the two network panels as manuscript Panels D and E
# ------------------------------------------------------------------------------

panel_D <- make_core_plot(
    network_object = D1_layout,
    panel_label = "D",
    title_text = "DNMT3A-centered RPPA core network"
)

panel_E <- make_core_plot(
    network_object = D2_layout,
    panel_label = "E",
    title_text = "DNMT3A-centered VIPER TF core network"
)


# ------------------------------------------------------------------------------
# Compact RPPA legend for Panels B-C
# ------------------------------------------------------------------------------

make_rppa_legend <- function() {

    ggplot() +
        geom_point(aes(0.10, 0.5), size = 2.0,
                   color = CONCORDANCE_COLORS["Concordant positive"]) +
        annotate("text", x = 0.18, y = 0.5, label = "Concordant positive",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        geom_point(aes(1.08, 0.5), size = 2.0,
                   color = CONCORDANCE_COLORS["Concordant negative"]) +
        annotate("text", x = 1.16, y = 0.5, label = "Concordant negative",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        geom_point(aes(2.06, 0.5), size = 2.0,
                   color = CONCORDANCE_COLORS["Discordant significant"]) +
        annotate("text", x = 2.14, y = 0.5, label = "Discordant significant",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        geom_point(aes(3.18, 0.5), size = 2.0,
                   color = CONCORDANCE_COLORS["Partial / not significant"]) +
        annotate("text", x = 3.26, y = 0.5, label = "Partial / not significant",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        geom_point(aes(4.48, 0.5), shape = 16, size = 1.95, color = "#4C4C4C") +
        annotate("text", x = 4.56, y = 0.5, label = "DNMT3A expression",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        geom_point(aes(5.50, 0.5), shape = 17, size = 2.1, color = "#4C4C4C") +
        annotate("text", x = 5.58, y = 0.5, label = "cg21629895 methylation",
                 hjust = 0, family = FONT_FAMILY, size = 4.2 / 2.845, color = INK) +
        coord_cartesian(xlim = c(-0.15, 6.85), ylim = c(0, 1), clip = "off") +
        theme_void(base_family = FONT_FAMILY) +
        theme(plot.margin = margin(0, 0, 0, 0))
}

rppa_shared_legend <- make_rppa_legend()
network_shared_legend <- make_network_legend()


# ------------------------------------------------------------------------------
# Final Figure 5 geometry: 180 x 215 mm
#
# A     full width
#
# B C   paired row
#       one-line RPPA annotation strip
#
# D E   paired network row
#       one-line STRING annotation strip
#
# Blank vertical space is intentionally minimized. The two annotation strips
# are kept compact and the three major visual blocks retain balanced spacing.
# ------------------------------------------------------------------------------

row_BC <- cowplot::plot_grid(
    panel_B,
    panel_C,
    nrow = 1,
    rel_widths = c(
        0.48,
        0.52
    ),
    align = "hv",
    axis = "tblr"
)

row_DE <- cowplot::plot_grid(
    panel_D,
    panel_E,
    nrow = 1,
    rel_widths = c(
        1,
        1
    ),
    align = "hv"
)

FIGURE5_WIDTH_MM <- 180
FIGURE5_HEIGHT_MM <- 215
FIGURE5_DPI <- 600

# Build each section first so the vertical spacing is controlled locally.
middle_section <- cowplot::plot_grid(
    row_BC,
    rppa_shared_legend,
    ncol = 1,
    rel_heights = c(
        0.947,
        0.053
    ),
    align = "v",
    axis = "lr"
)

bottom_section <- cowplot::plot_grid(
    row_DE,
    network_shared_legend,
    ncol = 1,
    rel_heights = c(
        0.947,
        0.053
    ),
    align = "v",
    axis = "lr"
)

# Main vertical structure:
#   A             35.5%
#   B/C section   30.0%
#   D/E section   34.5%
#
# This keeps the major rows visually balanced while giving the networks enough
# height and reclaiming space previously consumed by two-line legends.
figure5_AE <- cowplot::plot_grid(
    panel_A,
    middle_section,
    bottom_section,
    ncol = 1,
    rel_heights = c(
        0.350,
        0.315,
        0.335
    ),
    align = "v",
    axis = "lr"
)


# ------------------------------------------------------------------------------
# Windows-safe output handling
# ------------------------------------------------------------------------------

if (.Platform$OS.type == "windows" && nchar(OUTPUT_DIR) > 170) {

    cat(
        "\nThe selected Figure 5 output path is long and may cause PDF errors:\n",
        OUTPUT_DIR,
        "\n\nPlease select a SHORT output folder, for example:\n",
        "C:/Users/hko/Desktop/Fig5\n\n",
        sep = ""
    )

    short_dir <- utils::choose.dir(
        default = file.path(
            Sys.getenv("USERPROFILE"),
            "Desktop"
        ),
        caption = "Choose a SHORT output folder for final Figure 5"
    )

    if (
        is.na(short_dir) ||
        !nzchar(short_dir)
    ) {
        stop("Short Figure 5 output-folder selection was cancelled.")
    }

    OUTPUT_DIR <- normalizePath(
        short_dir,
        winslash = "/",
        mustWork = TRUE
    )
}


safe_save_pdf <- function(
        plot_object,
        filename,
        width_mm,
        height_mm
) {

    out_dir <- dirname(filename)

    if (!dir.exists(out_dir)) {
        dir.create(
            out_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    full_path <- file.path(
        normalizePath(
            out_dir,
            winslash = "/",
            mustWork = TRUE
        ),
        basename(filename)
    )

    if (
        .Platform$OS.type == "windows" &&
        nchar(full_path) > 240
    ) {
        stop(
            paste0(
                "PDF path is too long (",
                nchar(full_path),
                " characters).\nChoose a shorter output folder.\n",
                full_path
            )
        )
    }

    pdf_open <- FALSE

    tryCatch(
        {
            grDevices::pdf(
                file = full_path,
                width = width_mm / 25.4,
                height = height_mm / 25.4,
                family = "sans",
                useDingbats = FALSE
            )

            pdf_open <- TRUE
            print(plot_object)
            grDevices::dev.off()
            pdf_open <- FALSE
        },
        error = function(e) {

            if (
                pdf_open &&
                grDevices::dev.cur() > 1
            ) {
                try(
                    grDevices::dev.off(),
                    silent = TRUE
                )
            }

            stop(
                paste0(
                    "PDF save failed:\n",
                    full_path,
                    "\n\n",
                    conditionMessage(e)
                )
            )
        }
    )
}


safe_save_png <- function(
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
        dpi = FIGURE5_DPI,
        bg = "white",
        limitsize = FALSE
    )
}


# ------------------------------------------------------------------------------
# Save FINAL MERGED Figure 5
# ------------------------------------------------------------------------------

merged_pdf <- file.path(
    OUTPUT_DIR,
    "Figure5_A-E_180x215_BOXEDLABEL_DOTANCHOR_FINAL.pdf"
)

merged_png <- file.path(
    OUTPUT_DIR,
    "Figure5_A-E_180x215_BOXEDLABEL_DOTANCHOR_FINAL_600dpi.png"
)

safe_save_pdf(
    figure5_AE,
    merged_pdf,
    FIGURE5_WIDTH_MM,
    FIGURE5_HEIGHT_MM
)

safe_save_png(
    figure5_AE,
    merged_png,
    FIGURE5_WIDTH_MM,
    FIGURE5_HEIGHT_MM
)


# ------------------------------------------------------------------------------
# Save individual manuscript panels
# ------------------------------------------------------------------------------

panel_dir <- file.path(
    OUTPUT_DIR,
    "Figure5_Panels"
)

dir.create(
    panel_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

safe_save_pdf(
    panel_A,
    file.path(panel_dir, "Figure5_A_RPPA_Heatmap.pdf"),
    180,
    75
)

safe_save_pdf(
    panel_B,
    file.path(panel_dir, "Figure5_B_RPPA_Concordance.pdf"),
    88,
    72
)

safe_save_pdf(
    panel_C,
    file.path(panel_dir, "Figure5_C_RPPA_EffectSummary.pdf"),
    92,
    72
)

safe_save_pdf(
    panel_D,
    file.path(panel_dir, "Figure5_D_RPPA_STRING.pdf"),
    90,
    72
)

safe_save_pdf(
    panel_E,
    file.path(panel_dir, "Figure5_E_VIPER_STRING.pdf"),
    90,
    72
)


# ------------------------------------------------------------------------------
# Save retained STRING core tables for reproducibility
# ------------------------------------------------------------------------------

data.table::fwrite(
    D1_layout$nodes,
    file.path(
        OUTPUT_DIR,
        "Figure5_D_RPPA_STRING_Core_Nodes.csv"
    )
)

data.table::fwrite(
    D1_layout$edges,
    file.path(
        OUTPUT_DIR,
        "Figure5_D_RPPA_STRING_Core_Edges.csv"
    )
)

data.table::fwrite(
    D2_layout$nodes,
    file.path(
        OUTPUT_DIR,
        "Figure5_E_VIPER_STRING_Core_Nodes.csv"
    )
)

data.table::fwrite(
    D2_layout$edges,
    file.path(
        OUTPUT_DIR,
        "Figure5_E_VIPER_STRING_Core_Edges.csv"
    )
)


# ------------------------------------------------------------------------------
# Final console report
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FIGURE 5 A-E COMPLETE\n")
cat("============================================================\n")

cat("\nFinal layout:\n")
cat("  A = shared concordant RPPA heatmap\n")
cat("  B = RPPA differential effect-size concordance\n")
cat("  C = concordant RPPA paired-effect summary\n")
cat("  D = DNMT3A-centered RPPA STRING core network\n")
cat("  E = DNMT3A-centered VIPER TF STRING core network\n")

cat("\nMerged figure geometry:\n")
cat("  180 mm x 215 mm\n")
cat("  PDF + 600-dpi PNG\n")
cat("  RPPA B-C annotation: single horizontal row\n")
cat("  STRING D-E annotation: single horizontal row\n")
cat("  Node-size/STRING-degree note removed\n")
cat("  Panel B uses boxed labels, short arrow connectors, and exact labeled-dot anchors\n")

cat("\nSTRING core definition:\n")
cat("  score >= ", CORE_SELECTION_SCORE, "\n", sep = "")
cat("  shortest-path distance from DNMT3A <= ", CORE_ORDER, "\n", sep = "")

cat("\nPanel D network:\n")
cat("  nodes = ", nrow(D1_layout$nodes), "\n", sep = "")
cat("  associations = ", nrow(D1_layout$edges), "\n", sep = "")

cat("\nPanel E network:\n")
cat("  nodes = ", nrow(D2_layout$nodes), "\n", sep = "")
cat("  associations = ", nrow(D2_layout$edges), "\n", sep = "")

cat("\nMerged PDF:\n")
cat(merged_pdf, "\n")

cat("\nMerged PNG:\n")
cat(merged_png, "\n")

cat("\n============================================================\n")
