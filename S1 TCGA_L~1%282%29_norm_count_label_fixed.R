
# =============================================================================
# TCGA-LGG DNMT3A Supplementary Figure 1
# All 10 additional promoter CpG–DNMT3A expression correlations
#
# Purpose:
#   Main Figure 1 already displays cg21629895 vs DNMT3A expression.
#   Supplementary Figure 1 therefore shows the remaining 10 promoter-associated
#   CpGs from the predefined 11-CpG DNMT3A promoter set.
#
# Promoter CpGs shown here:
#   cg15998962
#   cg03463641
#   cg02208653
#   cg21708767
#   cg13344237
#   cg06112956
#   cg13076778
#   cg26470599
#   cg22731525
#   cg07150430
#
# Main-Figure CpG excluded from this supplementary panel:
#   cg21629895
#
# Analysis:
#   - Pearson correlation between each CpG beta value and DNMT3A expression.
#   - Overall linear-model trend is fit using all complete samples but drawn
#     only across each panel's visible x-range for a cleaner presentation.
#   - Confidence ribbons are intentionally omitted because Pearson r and P value
#     already provide the inferential summary.
#   - Points are colored by Molecular_3Group, matching Main Figure 1.
#   - BH-adjusted P values across the 10 supplementary correlations are saved
#     in the statistics output.
#
# Outputs:
#   - Supplementary_Figure_1_All_10_CpG_Correlations.pdf
#   - Supplementary_Figure_1_All_10_CpG_Correlations_600dpi.png
#   - Supplementary_Figure_1_Correlation_Statistics.csv
#   - Supplementary_Figure_1_Correlation_Data.xlsx
#   - Individual panel PDFs/PNGs
#
# Input:
#   1) Final DNMT3A master CSV selected manually with file.choose().
#   2) HM450 annotation CSV selected manually with file.choose().
# Output:
#   Save location selected manually by the user.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install and load packages
# -----------------------------------------------------------------------------

required_packages <- c(
    "data.table",
    "ggplot2",
    "patchwork",
    "openxlsx",
    "scales",
    "cowplot"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
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
    library(patchwork)
    library(openxlsx)
    library(scales)
    library(cowplot)
})


# -----------------------------------------------------------------------------
# 1. Global publication settings
# -----------------------------------------------------------------------------

set.seed(20260914)

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 215
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"

FONT_AXIS_TEXT <- 5.0
FONT_AXIS_TITLE <- 5.2
FONT_STAT <- 4.7
FONT_PANEL_TAG <- 7.0
FONT_TITLE <- 5.6
FONT_LEGEND <- 5.0

POINT_SIZE <- 0.70
LINE_WIDTH <- 0.38

# Same molecular-group palette used in Main Figure 1.
GROUP_COLORS <- c(
    "IDHmut-codel" = "#4C78A8",
    "IDHmut-noncodel" = "#E39C37",
    "IDHwt" = "#C94C4C",
    "Unclassified" = "#9E9E9E"
)


# -----------------------------------------------------------------------------
# 2. Select input master file
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE 1: 10 CpG–EXPRESSION CORRELATIONS\n")
cat("============================================================\n")

cat("\nSelect the final TCGA-LGG DNMT3A master CSV file.\n")

input_file <- file.choose()

if (!file.exists(input_file)) {
    stop(
        paste0(
            "Selected input file does not exist: ",
            input_file
        )
    )
}

df <- data.table::fread(
    input_file,
    check.names = FALSE,
    data.table = TRUE
)

cat("\nInput file loaded successfully.\n")
cat("Rows: ", nrow(df), "\n", sep = "")
cat("Columns: ", ncol(df), "\n", sep = "")


# -----------------------------------------------------------------------------
# 2B. Select HM450 probe annotation file manually
# -----------------------------------------------------------------------------

cat(
    "\nSelect the DNMT3A HM450 probe annotation CSV file.\n",
    "Required columns: Name, UCSC_RefGene_Group, Relation_to_Island\n",
    sep = ""
)

annotation_file <- file.choose()

if (!file.exists(annotation_file)) {
    stop(
        paste0(
            "Selected annotation file does not exist: ",
            annotation_file
        )
    )
}

annotation_df <- data.table::fread(
    annotation_file,
    check.names = FALSE,
    data.table = TRUE
)

required_annotation_cols <- c(
    "Name",
    "UCSC_RefGene_Group",
    "Relation_to_Island"
)

missing_annotation_cols <- setdiff(
    required_annotation_cols,
    names(annotation_df)
)

if (length(missing_annotation_cols) > 0) {
    stop(
        paste0(
            "Annotation file is missing required columns:\n",
            paste(
                missing_annotation_cols,
                collapse = "\n"
            )
        )
    )
}

cat("\nAnnotation file loaded successfully.\n")
cat("Rows: ", nrow(annotation_df), "\n", sep = "")


# -----------------------------------------------------------------------------
# 3. Select output directory manually
#
# The user chooses the folder directly.
# The script creates one subfolder inside the selected location:
#   DNMT3A_SuppFig1_10CpG_Correlations
# -----------------------------------------------------------------------------

cat(
    "\nNow select the folder where you want to save all Supplementary Figure 1 outputs.\n"
)

if (.Platform$OS.type == "windows") {

    selected_output <- utils::choose.dir(
        default = Sys.getenv("USERPROFILE"),
        caption = "Choose where to save Supplementary Figure 1 outputs"
    )

    if (
        is.na(selected_output) ||
        !nzchar(selected_output)
    ) {
        stop("Output-folder selection was cancelled.")
    }

} else {

    cat(
        "\nNon-Windows system detected.\n",
        "Please enter the full output folder path and press Enter:\n",
        sep = ""
    )

    selected_output <- readline()

    if (
        !nzchar(selected_output) ||
        !dir.exists(selected_output)
    ) {
        stop("A valid output folder was not selected.")
    }
}

output_dir <- file.path(
    selected_output,
    "DNMT3A_SuppFig1_10CpG_Correlations"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
    stop(
        paste0(
            "Could not create output directory: ",
            output_dir
        )
    )
}

cat("\nAll outputs will be saved to:\n")
cat(output_dir, "\n")


# -----------------------------------------------------------------------------
# 4. Define the promoter CpG set
# -----------------------------------------------------------------------------

expression_col <- "DNMT3A_expression"
subtype_col <- "Molecular_3Group"

promoter_probes_11 <- c(
    "cg15998962",
    "cg03463641",
    "cg02208653",
    "cg21708767",
    "cg13344237",
    "cg06112956",
    "cg13076778",
    "cg21629895",
    "cg26470599",
    "cg22731525",
    "cg07150430"
)

main_figure_cpg <- "cg21629895"

supplementary_probes_10 <- setdiff(
    promoter_probes_11,
    main_figure_cpg
)

cat("\nSupplementary Figure 1 CpGs:\n")
print(
    supplementary_probes_10
)


# -----------------------------------------------------------------------------
# 4B. Build concise CpG annotation labels
#
# UCSC_RefGene_Group may contain repeated transcript annotations, for example:
#   TSS1500;TSS200;TSS200
# These are reduced to unique labels:
#   TSS1500/TSS200
#
# Relation_to_Island labels are converted to manuscript-friendly text:
#   S_Shore -> S shore
#   N_Shore -> N shore
#   S_Shelf -> S shelf
#   N_Shelf -> N shelf
#   OpenSea -> Open sea
# -----------------------------------------------------------------------------

collapse_unique_semicolon <- function(x) {

    if (
        is.na(x) ||
        !nzchar(
            trimws(
                as.character(x)
            )
        )
    ) {
        return("NA")
    }

    parts <- trimws(
        unlist(
            strsplit(
                as.character(x),
                ";",
                fixed = TRUE
            )
        )
    )

    parts <- parts[
        nzchar(parts)
    ]

    paste(
        unique(parts),
        collapse = "/"
    )
}

format_island_relation <- function(x) {

    if (
        is.na(x) ||
        !nzchar(
            trimws(
                as.character(x)
            )
        )
    ) {
        return("NA")
    }

    out <- as.character(x)

    out <- gsub(
        "_",
        " ",
        out,
        fixed = TRUE
    )

    out <- gsub(
        "OpenSea",
        "Open sea",
        out,
        fixed = TRUE
    )

    out <- gsub(
        "S Shore",
        "S shore",
        out,
        fixed = TRUE
    )

    out <- gsub(
        "N Shore",
        "N shore",
        out,
        fixed = TRUE
    )

    out <- gsub(
        "S Shelf",
        "S shelf",
        out,
        fixed = TRUE
    )

    out <- gsub(
        "N Shelf",
        "N shelf",
        out,
        fixed = TRUE
    )

    out
}

probe_annotation <- annotation_df[
    Name %in% supplementary_probes_10,
    .(
        CpG = Name,
        UCSC_RefGene_Group = vapply(
            UCSC_RefGene_Group,
            collapse_unique_semicolon,
            FUN.VALUE = character(1)
        ),
        Relation_to_Island = vapply(
            Relation_to_Island,
            format_island_relation,
            FUN.VALUE = character(1)
        )
    )
]

probe_annotation <- unique(
    probe_annotation,
    by = "CpG"
)

missing_probe_annotations <- setdiff(
    supplementary_probes_10,
    probe_annotation$CpG
)

if (length(missing_probe_annotations) > 0) {
    stop(
        paste0(
            "The following Supplementary Figure 1 probes were not found in the annotation file:\n",
            paste(
                missing_probe_annotations,
                collapse = "\n"
            )
        )
    )
}

probe_annotation[
    ,
    Plot_Label := paste0(
        CpG,
        "\n",
        UCSC_RefGene_Group,
        ", ",
        Relation_to_Island
    )
]

cat("\nCpG annotations used in Supplementary Figure 1:\n")
print(
    probe_annotation[
        match(
            supplementary_probes_10,
            CpG
        )
    ]
)


# -----------------------------------------------------------------------------
# 5. Verify required columns
# -----------------------------------------------------------------------------

required_columns <- c(
    expression_col,
    subtype_col,
    supplementary_probes_10
)

missing_columns <- setdiff(
    required_columns,
    names(df)
)

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "Required columns are missing from the master file:\n",
            paste(
                missing_columns,
                collapse = "\n"
            )
        )
    )
}


# -----------------------------------------------------------------------------
# 6. Convert expression and CpGs to numeric
# -----------------------------------------------------------------------------

df[
    ,
    (expression_col) := suppressWarnings(
        as.numeric(
            get(expression_col)
        )
    )
]

for (probe in supplementary_probes_10) {

    df[
        ,
        (probe) := suppressWarnings(
            as.numeric(
                get(probe)
            )
        )
    ]
}


# -----------------------------------------------------------------------------
# 7. Prepare molecular subtype factor
# -----------------------------------------------------------------------------

molecular_levels <- c(
    "IDH-mutant/1p19q-codeleted",
    "IDH-mutant/1p19q-non-codeleted",
    "IDH-wildtype"
)

molecular_labels <- c(
    "IDHmut-codel",
    "IDHmut-noncodel",
    "IDHwt"
)

df[
    ,
    Molecular_3Group_Plot := factor(
        get(subtype_col),
        levels = molecular_levels,
        labels = molecular_labels
    )
]

# Keep all samples in the overall Pearson analysis.
# Samples without one of the three predefined groups remain in the analysis and
# are colored as Unclassified in the figure.
df[
    ,
    Molecular_3Group_Scatter := as.character(
        Molecular_3Group_Plot
    )
]

df[
    is.na(Molecular_3Group_Scatter),
    Molecular_3Group_Scatter := "Unclassified"
]

df[
    ,
    Molecular_3Group_Scatter := factor(
        Molecular_3Group_Scatter,
        levels = c(
            "IDHmut-codel",
            "IDHmut-noncodel",
            "IDHwt",
            "Unclassified"
        )
    )
]


# -----------------------------------------------------------------------------
# 8. Helper: format P values
# -----------------------------------------------------------------------------

format_p <- function(p) {

    if (is.na(p)) {
        return("P = NA")
    }

    if (p < 2.2e-16) {
        return("P < 2.2e-16")
    }

    if (p < 0.001) {
        return(
            paste0(
                "P = ",
                format(
                    p,
                    scientific = TRUE,
                    digits = 2
                )
            )
        )
    }

    paste0(
        "P = ",
        format(
            p,
            digits = 3,
            nsmall = 3
        )
    )
}


# -----------------------------------------------------------------------------
# 9. Publication theme
# -----------------------------------------------------------------------------

publication_theme <- function() {

    theme_classic(
        base_family = FONT_FAMILY,
        base_size = FONT_AXIS_TEXT
    ) +
        theme(
            axis.title = element_text(
                family = FONT_FAMILY,
                face = "plain",
                size = FONT_AXIS_TITLE,
                color = "black"
            ),
            axis.text = element_text(
                family = FONT_FAMILY,
                size = FONT_AXIS_TEXT,
                color = "black"
            ),
            axis.line = element_line(
                linewidth = 0.32,
                color = "black"
            ),
            axis.ticks = element_line(
                linewidth = 0.28,
                color = "black"
            ),
            axis.ticks.length = grid::unit(
                1.2,
                "mm"
            ),
            legend.title = element_blank(),
            legend.text = element_text(
                family = FONT_FAMILY,
                size = FONT_LEGEND
            ),
            legend.key.width = grid::unit(
                3.0,
                "mm"
            ),
            plot.title = element_text(
                family = FONT_FAMILY,
                face = "plain",
                size = FONT_TITLE,
                hjust = 0.5,
                lineheight = 0.92,
                margin = margin(
                    b = 1.5,
                    unit = "mm"
                )
            ),
            plot.tag = element_text(
                family = FONT_FAMILY,
                face = "bold",
                size = FONT_PANEL_TAG
            ),
            plot.margin = margin(
                2.2,
                2.0,
                2.8,
                2.4,
                unit = "mm"
            )
        )
}


# -----------------------------------------------------------------------------
# 10. Run all 10 Pearson correlations
# -----------------------------------------------------------------------------

run_probe_correlation <- function(
        data,
        probe,
        expression_col
) {

    analysis_df <- data[
        !is.na(get(probe)) &
            !is.na(get(expression_col))
    ]

    if (nrow(analysis_df) < 3) {

        return(
            data.table::data.table(
                CpG = probe,
                N = nrow(analysis_df),
                Pearson_R = NA_real_,
                P_Value = NA_real_,
                CI_Lower = NA_real_,
                CI_Upper = NA_real_
            )
        )
    }

    test <- stats::cor.test(
        analysis_df[[probe]],
        analysis_df[[expression_col]],
        method = "pearson",
        use = "complete.obs"
    )

    data.table::data.table(
        CpG = probe,
        N = nrow(analysis_df),
        Pearson_R = unname(
            test$estimate
        ),
        P_Value = test$p.value,
        CI_Lower = test$conf.int[1],
        CI_Upper = test$conf.int[2]
    )
}


correlation_results <- data.table::rbindlist(
    lapply(
        supplementary_probes_10,
        function(probe) {
            run_probe_correlation(
                df,
                probe,
                expression_col
            )
        }
    ),
    fill = TRUE
)

# Multiple-testing correction across the 10 supplementary CpG-expression tests.
correlation_results[
    ,
    BH_FDR := stats::p.adjust(
        P_Value,
        method = "BH"
    )
]

correlation_results[
    ,
    Significant_BH_0.05 := ifelse(
        !is.na(BH_FDR) &
            BH_FDR < 0.05,
        "Yes",
        "No"
    )
]


correlation_results <- merge(
    correlation_results,
    probe_annotation[
        ,
        .(
            CpG,
            UCSC_RefGene_Group,
            Relation_to_Island
        )
    ],
    by = "CpG",
    all.x = TRUE,
    sort = FALSE
)

# Restore the predefined Supplementary Figure order after merge.
correlation_results[
    ,
    CpG_Order := match(
        CpG,
        supplementary_probes_10
    )
]

data.table::setorder(
    correlation_results,
    CpG_Order
)

correlation_results[
    ,
    CpG_Order := NULL
]

cat("\nPearson correlation results:\n")
print(
    correlation_results
)


# -----------------------------------------------------------------------------
# 10B. Adaptive x-axis limits and breaks for each CpG
#
# CLEAN VISUALIZATION RULE:
#   The statistical correlation always uses ALL complete samples.
#   The plotting window is chosen independently for each CpG to avoid a few
#   extreme methylation values stretching the entire x-axis.
#
#   Primary rule:
#       Tukey inner-fence window = Q1 - 1.5*IQR to Q3 + 1.5*IQR
#
#   Safety rule:
#       The window is not allowed to be narrower than the 1st–99th percentile.
#
#   Final limits are constrained to biologically valid beta-value range [0, 1].
#
# This changes only the visual zoom, not Pearson r, P value, N, or the regression
# fit used for the displayed trend line.
# -----------------------------------------------------------------------------

get_adaptive_x_axis <- function(x) {

    x <- x[
        is.finite(x)
    ]

    if (length(x) < 3) {
        return(
            list(
                limits = range(
                    x,
                    na.rm = TRUE
                ),
                breaks = waiver()
            )
        )
    }

    q01 <- stats::quantile(
        x,
        0.01,
        na.rm = TRUE,
        names = FALSE
    )

    q99 <- stats::quantile(
        x,
        0.99,
        na.rm = TRUE,
        names = FALSE
    )

    q1 <- stats::quantile(
        x,
        0.25,
        na.rm = TRUE,
        names = FALSE
    )

    q3 <- stats::quantile(
        x,
        0.75,
        na.rm = TRUE,
        names = FALSE
    )

    iqr_value <- q3 - q1

    tukey_lower <- q1 - 1.5 * iqr_value
    tukey_upper <- q3 + 1.5 * iqr_value

    lower <- min(
        q01,
        tukey_lower
    )

    upper <- max(
        q99,
        tukey_upper
    )

    lower <- max(
        0,
        lower
    )

    upper <- min(
        1,
        upper
    )

    span <- upper - lower

    if (
        !is.finite(span) ||
        span <= 0
    ) {

        lower <- max(
            0,
            min(
                x,
                na.rm = TRUE
            )
        )

        upper <- min(
            1,
            max(
                x,
                na.rm = TRUE
            )
        )

        span <- upper - lower
    }

    # Small visual padding only.
    padding <- span * 0.025

    lower <- max(
        0,
        lower - padding
    )

    upper <- min(
        1,
        upper + padding
    )

    axis_breaks <- scales::pretty_breaks(
        n = 4
    )(
        c(
            lower,
            upper
        )
    )

    axis_breaks <- axis_breaks[
        axis_breaks >= lower &
            axis_breaks <= upper
    ]

    list(
        limits = c(
            lower,
            upper
        ),
        breaks = axis_breaks
    )
}


# -----------------------------------------------------------------------------
# 11. Create each supplementary scatter panel
# -----------------------------------------------------------------------------

make_scatter_panel <- function(
        data,
        probe,
        correlation_row,
        panel_letter
) {

    plot_df <- data[
        !is.na(get(probe)) &
            !is.na(get(expression_col))
    ]

    r_value <- correlation_row$Pearson_R[1]
    p_value <- correlation_row$P_Value[1]
    n_value <- correlation_row$N[1]

    annotation_text <- paste0(
        "Pearson r = ",
        sprintf(
            "%.3f",
            r_value
        ),
        "\n",
        format_p(
            p_value
        ),
        ", n = ",
        n_value
    )

    annotation_row <- probe_annotation[
        CpG == probe
    ]

    probe_title <- annotation_row$Plot_Label[1]

    x_axis <- get_adaptive_x_axis(
        plot_df[[probe]]
    )

    # -------------------------------------------------------------------------
    # Clean trend line:
    #
    # Fit the linear model using ALL complete samples, exactly as before.
    # Then predict only across the visible x-axis interval.
    #
    # This prevents the regression line / confidence ribbon from visually
    # dominating sparse tails or extending through large empty regions.
    # No confidence ribbon is drawn in the figure; inferential information is
    # already reported by Pearson r and P value.
    # -------------------------------------------------------------------------

    lm_formula <- stats::as.formula(
        paste0(
            "`",
            expression_col,
            "` ~ `",
            probe,
            "`"
        )
    )

    lm_fit <- stats::lm(
        lm_formula,
        data = as.data.frame(
            plot_df
        )
    )

    trend_x <- seq(
        x_axis$limits[1],
        x_axis$limits[2],
        length.out = 200
    )

    trend_df <- data.frame(
        x_value = trend_x
    )

    prediction_input <- setNames(
        data.frame(
            trend_x
        ),
        probe
    )

    trend_df$y_value <- as.numeric(
        stats::predict(
            lm_fit,
            newdata = prediction_input
        )
    )

    ggplot(
        plot_df,
        aes(
            x = .data[[probe]],
            y = .data[[expression_col]]
        )
    ) +
        geom_point(
            aes(
                color = Molecular_3Group_Scatter
            ),
            size = POINT_SIZE,
            alpha = 0.58,
            show.legend = TRUE
        ) +
        geom_line(
            data = trend_df,
            aes(
                x = x_value,
                y = y_value
            ),
            inherit.aes = FALSE,
            linewidth = 0.34,
            color = "#222222",
            lineend = "round"
        ) +
        scale_color_manual(
            values = GROUP_COLORS,
            drop = FALSE
        ) +
        annotate(
            "text",
            x = -Inf,
            y = Inf,
            label = annotation_text,
            hjust = -0.02,
            vjust = 1.04,
            family = FONT_FAMILY,
            size = FONT_STAT / 2.845,
            lineheight = 0.92,
            color = "black"
        ) +
        scale_x_continuous(
            breaks = x_axis$breaks
        ) +
        coord_cartesian(
            xlim = x_axis$limits,
            clip = "off"
        ) +
        labs(
            title = probe_title,
            x = "CpG methylation\n(beta value)",
            y = "DNMT3A expression\n(log2(norm_count + 1))",
            tag = panel_letter
        ) +
        publication_theme() +
        theme(
            legend.position = "none",
            plot.title = element_text(
                family = FONT_FAMILY,
                face = "plain",
                size = FONT_TITLE,
                hjust = 0.5,
                lineheight = 0.90,
                margin = margin(
                    b = 1.6,
                    unit = "mm"
                )
            )
        )
}

panel_letters <- LETTERS[
    seq_along(
        supplementary_probes_10
    )
]

panel_list <- vector(
    "list",
    length(
        supplementary_probes_10
    )
)

names(
    panel_list
) <- supplementary_probes_10

for (i in seq_along(supplementary_probes_10)) {

    probe <- supplementary_probes_10[i]

    correlation_row <- correlation_results[
        CpG == probe
    ]

    panel_list[[i]] <- make_scatter_panel(
        data = df,
        probe = probe,
        correlation_row = correlation_row,
        panel_letter = panel_letters[i]
    )
}


# -----------------------------------------------------------------------------
# 11B. Create ONE shared molecular-subtype legend
#
# The legend is intentionally shown only once at the very bottom of the merged
# Supplementary Figure 1. This avoids repeated legends between rows and gives a
# much cleaner Nature-style presentation.
# -----------------------------------------------------------------------------

legend_source <- ggplot(
    data.frame(
        x = 1:4,
        y = 1:4,
        Group = factor(
            c(
                "IDHmut-codel",
                "IDHmut-noncodel",
                "IDHwt",
                "Unclassified"
            ),
            levels = c(
                "IDHmut-codel",
                "IDHmut-noncodel",
                "IDHwt",
                "Unclassified"
            )
        )
    ),
    aes(
        x = x,
        y = y,
        color = Group
    )
) +
    geom_point(
        size = 1.6,
        alpha = 1
    ) +
    scale_color_manual(
        values = GROUP_COLORS,
        drop = FALSE
    ) +
    guides(
        color = guide_legend(
            title = NULL,
            nrow = 1,
            byrow = TRUE,
            override.aes = list(
                size = 1.6,
                alpha = 1
            )
        )
    ) +
    theme_void(
        base_family = FONT_FAMILY
    ) +
    theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(
            family = FONT_FAMILY,
            size = FONT_LEGEND,
            color = "black"
        ),
        legend.key.width = grid::unit(
            3.2,
            "mm"
        ),
        legend.key.height = grid::unit(
            2.0,
            "mm"
        ),
        legend.spacing.x = grid::unit(
            1.2,
            "mm"
        ),
        legend.margin = margin(
            t = 0,
            r = 0,
            b = 0,
            l = 0,
            unit = "mm"
        )
    )

shared_legend <- cowplot::get_legend(
    legend_source
)


# -----------------------------------------------------------------------------
# 12. Assemble Supplementary Figure 1
#
# Final publication layout:
#   180 mm x 215 mm
#   3 columns x 4 rows
#
# Row 1: A B C
# Row 2: D E F
# Row 3: G H I
# Row 4: J [blank] [blank]
#
# IMPORTANT:
# The panel grid and the single shared legend are assembled with cowplot rather
# than nesting patchwork layouts. This avoids the patchwork wrap_dims() error:
# "Need 5 panels, but together nrow and ncol only provide 2."
# -----------------------------------------------------------------------------

blank_panel <- grid::nullGrob()

figure_grid <- cowplot::plot_grid(
    panel_list[[1]],
    panel_list[[2]],
    panel_list[[3]],
    panel_list[[4]],
    panel_list[[5]],
    panel_list[[6]],
    panel_list[[7]],
    panel_list[[8]],
    panel_list[[9]],
    panel_list[[10]],
    blank_panel,
    blank_panel,
    ncol = 3,
    nrow = 4,
    align = "hv",
    axis = "tblr",
    rel_widths = c(
        1,
        1,
        1
    ),
    rel_heights = c(
        1,
        1,
        1,
        1
    )
)

supplementary_figure_1 <- cowplot::plot_grid(
    figure_grid,
    shared_legend,
    ncol = 1,
    rel_heights = c(
        0.955,
        0.045
    ),
    align = "v",
    axis = "lr"
)


# -----------------------------------------------------------------------------
# 13. Robust file-saving helpers
# -----------------------------------------------------------------------------

safe_save_pdf <- function(
        plot_object,
        filename,
        width_mm,
        height_mm
) {

    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4

    result <- tryCatch(
        {

            grDevices::pdf(
                file = filename,
                width = width_in,
                height = height_in,
                family = "sans",
                useDingbats = FALSE,
                onefile = TRUE
            )

            print(
                plot_object
            )

            grDevices::dev.off()

            TRUE
        },
        error = function(e) {

            if (grDevices::dev.cur() > 1) {
                try(
                    grDevices::dev.off(),
                    silent = TRUE
                )
            }

            message(
                "\nPDF save failed:\n",
                filename,
                "\nError: ",
                conditionMessage(e)
            )

            FALSE
        }
    )

    if (!result) {

        fallback_file <- file.path(
            tempdir(),
            basename(filename)
        )

        grDevices::pdf(
            file = fallback_file,
            width = width_in,
            height = height_in,
            family = "sans",
            useDingbats = FALSE,
            onefile = TRUE
        )

        print(
            plot_object
        )

        grDevices::dev.off()

        message(
            "\nFallback PDF saved successfully:\n",
            fallback_file
        )
    }
}


safe_save_png <- function(
        plot_object,
        filename,
        width_mm,
        height_mm
) {

    tryCatch(
        {

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

        },
        error = function(e) {

            fallback_file <- file.path(
                tempdir(),
                basename(filename)
            )

            ggplot2::ggsave(
                filename = fallback_file,
                plot = plot_object,
                width = width_mm,
                height = height_mm,
                units = "mm",
                dpi = FIGURE_DPI,
                bg = "white",
                limitsize = FALSE
            )

            message(
                "\nFallback PNG saved successfully:\n",
                fallback_file
            )
        }
    )
}


# -----------------------------------------------------------------------------
# 14. Save combined Supplementary Figure 1
# -----------------------------------------------------------------------------

combined_pdf <- file.path(
    output_dir,
    "SuppFig1_10CpG_Correlations_180x215_CLEAN_FIXED.pdf"
)

combined_png <- file.path(
    output_dir,
    "SuppFig1_10CpG_Correlations_180x215_CLEAN_FIXED_600dpi.png"
)

safe_save_pdf(
    supplementary_figure_1,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

safe_save_png(
    supplementary_figure_1,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)


# -----------------------------------------------------------------------------
# 15. Save individual scatter panels
# -----------------------------------------------------------------------------

individual_panel_dir <- file.path(
    output_dir,
    "Individual_Panels"
)

dir.create(
    individual_panel_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

for (i in seq_along(supplementary_probes_10)) {

    probe <- supplementary_probes_10[i]

    safe_save_pdf(
        panel_list[[i]],
        file.path(
            individual_panel_dir,
            paste0(
                "SuppFig1_",
                panel_letters[i],
                "_",
                probe,
                ".pdf"
            )
        ),
        86,
        62
    )

    safe_save_png(
        panel_list[[i]],
        file.path(
            individual_panel_dir,
            paste0(
                "SuppFig1_",
                panel_letters[i],
                "_",
                probe,
                "_600dpi.png"
            )
        ),
        86,
        62
    )
}


# -----------------------------------------------------------------------------
# 16. Save correlation statistics CSV
# -----------------------------------------------------------------------------

statistics_csv <- file.path(
    output_dir,
    "Supplementary_Figure_1_Correlation_Statistics.csv"
)

data.table::fwrite(
    correlation_results,
    statistics_csv,
    na = ""
)


# -----------------------------------------------------------------------------
# 17. Prepare sample-level long-form correlation data
# -----------------------------------------------------------------------------

id_columns <- intersect(
    c(
        "Patient_ID",
        "Sample_Core",
        "Molecular_3Group"
    ),
    names(df)
)

correlation_long <- data.table::melt(
    df,
    id.vars = c(
        id_columns,
        expression_col
    ),
    measure.vars = supplementary_probes_10,
    variable.name = "CpG",
    value.name = "Methylation_Beta",
    variable.factor = FALSE
)

data.table::setnames(
    correlation_long,
    expression_col,
    "DNMT3A_expression"
)

correlation_long <- correlation_long[
    !is.na(Methylation_Beta) &
        !is.na(DNMT3A_expression)
]


# -----------------------------------------------------------------------------
# 18. Save publication-friendly Excel workbook
# -----------------------------------------------------------------------------

excel_file <- file.path(
    output_dir,
    "Supplementary_Figure_1_Correlation_Data.xlsx"
)

wb <- openxlsx::createWorkbook()

# Styles.
title_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 11,
    textDecoration = "bold",
    halign = "left"
)

header_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    textDecoration = "bold",
    fgFill = "#E7E7E7",
    halign = "center",
    valign = "center",
    border = "Bottom",
    borderColour = "#000000"
)

body_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    valign = "center"
)

number_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    numFmt = "0.000"
)

p_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    numFmt = "0.000E+00"
)

# README.
openxlsx::addWorksheet(
    wb,
    "README",
    gridLines = FALSE
)

readme_df <- data.frame(
    Item = c(
        "Figure",
        "Cohort",
        "Expression variable",
        "Methylation variable",
        "Correlation method",
        "Multiple-testing correction",
        "Main-Figure CpG excluded",
        "Supplementary CpGs",
        "Interpretation"
    ),
    Value = c(
        "Supplementary Figure 1",
        "TCGA-LGG",
        "DNMT3A_expression",
        "HumanMethylation450 beta value",
        "Pearson correlation",
        "Benjamini-Hochberg across 10 CpG-expression correlations",
        "cg21629895",
        paste(
            supplementary_probes_10,
            collapse = ", "
        ),
        "Correlation represents association, not causation."
    ),
    stringsAsFactors = FALSE
)

openxlsx::writeData(
    wb,
    "README",
    readme_df,
    headerStyle = header_style
)

openxlsx::setColWidths(
    wb,
    "README",
    cols = 1,
    widths = 31
)

openxlsx::setColWidths(
    wb,
    "README",
    cols = 2,
    widths = 95
)

# Correlation summary.
openxlsx::addWorksheet(
    wb,
    "Correlation_Summary",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Correlation_Summary",
    "Supplementary Figure 1 | All 10 promoter CpG-expression correlations",
    startRow = 1,
    startCol = 1
)

openxlsx::mergeCells(
    wb,
    "Correlation_Summary",
    cols = 1:ncol(correlation_results),
    rows = 1
)

openxlsx::addStyle(
    wb,
    "Correlation_Summary",
    title_style,
    rows = 1,
    cols = 1:ncol(correlation_results),
    gridExpand = TRUE
)

openxlsx::writeData(
    wb,
    "Correlation_Summary",
    correlation_results,
    startRow = 3,
    startCol = 1,
    headerStyle = header_style
)

openxlsx::addStyle(
    wb,
    "Correlation_Summary",
    body_style,
    rows = 4:(3 + nrow(correlation_results)),
    cols = 1:ncol(correlation_results),
    gridExpand = TRUE
)

openxlsx::addStyle(
    wb,
    "Correlation_Summary",
    number_style,
    rows = 4:(3 + nrow(correlation_results)),
    cols = c(
        3,
        5,
        6
    ),
    gridExpand = TRUE
)

openxlsx::addStyle(
    wb,
    "Correlation_Summary",
    p_style,
    rows = 4:(3 + nrow(correlation_results)),
    cols = c(
        4,
        7
    ),
    gridExpand = TRUE
)

openxlsx::setColWidths(
    wb,
    "Correlation_Summary",
    cols = 1:ncol(correlation_results),
    widths = "auto"
)

openxlsx::freezePane(
    wb,
    "Correlation_Summary",
    firstActiveRow = 4
)

# Long-form data.
openxlsx::addWorksheet(
    wb,
    "Sample_Level_Data",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Sample_Level_Data",
    correlation_long,
    headerStyle = header_style,
    withFilter = TRUE
)

openxlsx::freezePane(
    wb,
    "Sample_Level_Data",
    firstRow = TRUE
)

openxlsx::setColWidths(
    wb,
    "Sample_Level_Data",
    cols = 1:ncol(correlation_long),
    widths = "auto"
)

openxlsx::saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
)


# -----------------------------------------------------------------------------
# 19. QC output
# -----------------------------------------------------------------------------

probe_qc <- data.table::rbindlist(
    lapply(
        supplementary_probes_10,
        function(probe) {

            x_axis_qc <- get_adaptive_x_axis(
                df[[probe]]
            )

            ann <- probe_annotation[
                CpG == probe
            ]

            data.table::data.table(
                CpG = probe,
                UCSC_RefGene_Group = ann$UCSC_RefGene_Group[1],
                Relation_to_Island = ann$Relation_to_Island[1],
                Plot_X_Min = x_axis_qc$limits[1],
                Plot_X_Max = x_axis_qc$limits[2],
                N_Total = nrow(df),
                N_Methylation_Available = sum(
                    !is.na(
                        df[[probe]]
                    )
                ),
                N_Expression_Available = sum(
                    !is.na(
                        df[[expression_col]]
                    )
                ),
                N_Complete_Pairs = sum(
                    !is.na(
                        df[[probe]]
                    ) &
                        !is.na(
                            df[[expression_col]]
                        )
                ),
                Mean_Beta = mean(
                    df[[probe]],
                    na.rm = TRUE
                ),
                Median_Beta = median(
                    df[[probe]],
                    na.rm = TRUE
                ),
                Min_Beta = min(
                    df[[probe]],
                    na.rm = TRUE
                ),
                Max_Beta = max(
                    df[[probe]],
                    na.rm = TRUE
                )
            )
        }
    ),
    fill = TRUE
)

qc_csv <- file.path(
    output_dir,
    "Supplementary_Figure_1_Probe_QC.csv"
)

data.table::fwrite(
    probe_qc,
    qc_csv,
    na = ""
)


# -----------------------------------------------------------------------------
# 20. Final console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE 1 COMPLETE\n")
cat("============================================================\n")

cat("\nMain Figure 1 CpG excluded:\n")
cat("  ", main_figure_cpg, "\n", sep = "")

cat("\n10 supplementary promoter CpGs analyzed:\n")
print(
    supplementary_probes_10
)

cat("\nCorrelation method:\n")
cat("  Pearson correlation\n")

cat("\nMultiple-testing correction:\n")
cat("  Benjamini-Hochberg across the 10 supplementary correlations\n")

cat("\nFigure geometry:\n")
cat("  180 mm x 215 mm\n")
cat("  3 columns x 4 rows\n")
cat("  Final row: J left-aligned + two blank positions\n")
cat("  One shared legend only, placed at the bottom\n")
cat("  cowplot assembly used to avoid patchwork wrap_dims errors\n")

cat("\nCombined figure:\n")
cat("  ", combined_pdf, "\n", sep = "")
cat("  ", combined_png, "\n", sep = "")

cat("\nStatistics:\n")
cat("  ", statistics_csv, "\n", sep = "")

cat("\nExcel workbook:\n")
cat("  ", excel_file, "\n", sep = "")

cat("\nQC:\n")
cat("  ", qc_csv, "\n", sep = "")

cat("\nFinal correlation summary:\n")
print(
    correlation_results
)

cat("\n============================================================\n")
