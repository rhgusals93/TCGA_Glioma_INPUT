# =============================================================================
# TCGA-LGG DNMT3A Main Figure 1
# Nature-style 180-mm publication layout with subtype colors and sample-size labels
#
# Panels:
#   A. DNMT3A mRNA expression across three molecular groups
#   B. cg21629895 methylation across three molecular groups
#   C. Mean methylation across 11 promoter-associated CpGs across groups
#   D. cg21629895 methylation vs DNMT3A mRNA expression
#   E. Mean methylation across 11 promoter-associated CpGs vs DNMT3A expression
#
# Publication format:
#   - Final combined figure width: 180 mm
#   - Final analytical font sizes: approximately 5-7 pt
#   - Vector PDF + 600-dpi PNG
#
# Important Windows fixes:
#   - cairo_pdf is NOT used.
#   - PDF files are written with the standard grDevices::pdf device.
#   - All figure text uses the generic "sans" family to avoid invalid-font-type
#     errors on Windows PDF devices.
#   - The user can select a short output directory to avoid Windows/OneDrive
#     path-length and output-stream errors.
#
# Current data context:
#   - TCGA-LGG legacy HiSeqV2 DNMT3A expression
#   - TCGA-LGG legacy HumanMethylation450 beta values
#   - 11-probe DNMT3A promoter-associated mean methylation
#   - Molecular_3Group from the completed current master file
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install and load packages
# -----------------------------------------------------------------------------

required_packages <- c(
    "data.table",
    "ggplot2",
    "patchwork"
)

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

library(data.table)
library(ggplot2)
library(patchwork)


# -----------------------------------------------------------------------------
# 1. Global publication settings
# -----------------------------------------------------------------------------

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 148
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"

FONT_AXIS_TEXT <- 5.6
FONT_AXIS_TITLE <- 6.2
FONT_STAT <- 5.2
FONT_PANEL_TAG <- 6.8
FONT_LEGEND <- 5.2

POINT_SIZE <- 0.90
JITTER_SIZE <- 0.58
LINE_WIDTH <- 0.35


# -----------------------------------------------------------------------------
# 2. Select input file
# -----------------------------------------------------------------------------

cat("\nSelect the final TCGA-LGG DNMT3A master CSV file.\n")

input_file <- file.choose()

df <- data.table::fread(
    input_file,
    check.names = FALSE
)

cat("\nInput file loaded successfully.\n")
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n")


# -----------------------------------------------------------------------------
# 3. Select a SHORT output directory
#
# A short output path is recommended on Windows because long OneDrive paths can
# cause PDF device errors such as:
#   cairo error 'error while writing to output stream'
#
# This script does not use cairo_pdf, but a short path is still safer.
# -----------------------------------------------------------------------------

if (.Platform$OS.type == "windows") {

    cat(
        "\nSelect a SHORT output folder for Figure 1.\n",
        "A Desktop folder is recommended.\n",
        sep = ""
    )

    selected_output <- tryCatch(
        choose.dir(
            default = Sys.getenv("USERPROFILE"),
            caption = "Select a SHORT output folder for DNMT3A Figure 1"
        ),
        error = function(e) NA_character_
    )

} else {

    selected_output <- NA_character_
}


if (
    length(selected_output) == 0 ||
    is.na(selected_output) ||
    selected_output == ""
) {

    desktop_dir <- file.path(
        Sys.getenv("USERPROFILE"),
        "Desktop"
    )

    if (
        .Platform$OS.type == "windows" &&
        dir.exists(desktop_dir)
    ) {

        output_dir <- file.path(
            desktop_dir,
            "DNMT3A_Fig1_11Probe_Legacy"
        )

    } else {

        output_dir <- file.path(
            getwd(),
            "DNMT3A_Fig1_11Probe_Legacy"
        )
    }

} else {

    output_dir <- file.path(
        selected_output,
        "DNMT3A_Fig1_11Probe_Legacy"
    )
}


dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat(
    "Output path length:",
    nchar(normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
    "characters\n"
)


# -----------------------------------------------------------------------------
# 4. Define required variables
# -----------------------------------------------------------------------------

expression_col <- "DNMT3A_expression"
single_cpg_col <- "cg21629895"
promoter_mean_col <- "DNMT3A_Promoter_Mean_Methylation_11Probe"
subtype_col <- "Molecular_3Group"

required_cols <- c(
    expression_col,
    single_cpg_col,
    promoter_mean_col,
    subtype_col
)

missing_cols <- setdiff(
    required_cols,
    names(df)
)

if (length(missing_cols) > 0) {
    stop(
        paste0(
            "Required columns are missing: ",
            paste(missing_cols, collapse = ", ")
        )
    )
}


# -----------------------------------------------------------------------------
# 5. Prepare analysis data
# -----------------------------------------------------------------------------

df[[expression_col]] <- suppressWarnings(
    as.numeric(df[[expression_col]])
)

df[[single_cpg_col]] <- suppressWarnings(
    as.numeric(df[[single_cpg_col]])
)

df[[promoter_mean_col]] <- suppressWarnings(
    as.numeric(df[[promoter_mean_col]])
)

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

df[, Molecular_3Group_Plot := factor(
    get(subtype_col),
    levels = molecular_levels,
    labels = molecular_labels
)]

# Consistent molecular-group palette used across all Figure 1 panels.
GROUP_COLORS <- c(
    "IDHmut-codel" = "#4C78A8",
    "IDHmut-noncodel" = "#E39C37",
    "IDHwt" = "#C94C4C",
    "Unclassified" = "#9E9E9E"
)

# Keep all samples in correlation panels. Samples without a molecular subtype
# are shown as Unclassified but remain included in the overall Pearson analysis.
df[, Molecular_3Group_Scatter := as.character(Molecular_3Group_Plot)]
df[is.na(Molecular_3Group_Scatter), Molecular_3Group_Scatter := "Unclassified"]

df[, Molecular_3Group_Scatter := factor(
    Molecular_3Group_Scatter,
    levels = c(
        "IDHmut-codel",
        "IDHmut-noncodel",
        "IDHwt",
        "Unclassified"
    )
)]

cat("\nMolecular group counts:\n")
print(
    table(
        df$Molecular_3Group_Plot,
        useNA = "ifany"
    )
)


cat("\nCurrent Figure 1 input QC:\n")
cat("Expression column:", expression_col, "\n")
cat("Single CpG column:", single_cpg_col, "\n")
cat("Promoter mean column:", promoter_mean_col, "\n")

if ("Sample_Core" %in% names(df)) {
    cat(
        "Unique Sample_Core:",
        data.table::uniqueN(df$Sample_Core),
        "\n"
    )
}


# -----------------------------------------------------------------------------
# 6. Publication theme
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
                linewidth = 0.35,
                color = "black"
            ),
            axis.ticks = element_line(
                linewidth = 0.30,
                color = "black"
            ),
            axis.ticks.length = grid::unit(
                1.3,
                "mm"
            ),
            legend.title = element_text(
                family = FONT_FAMILY,
                size = FONT_LEGEND
            ),
            legend.text = element_text(
                family = FONT_FAMILY,
                size = FONT_LEGEND
            ),
            plot.tag = element_text(
                family = FONT_FAMILY,
                face = "bold",
                size = FONT_PANEL_TAG
            ),
            plot.margin = margin(
                3.5,
                4.0,
                4.5,
                4.0,
                unit = "mm"
            )
        )
}


# -----------------------------------------------------------------------------
# 7. Helper functions
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


run_group_statistics <- function(
        data,
        value_col,
        variable_label
) {

    analysis_df <- data[
        !is.na(Molecular_3Group_Plot) &
            !is.na(get(value_col))
    ]

    kw <- stats::kruskal.test(
        analysis_df[[value_col]] ~
            analysis_df$Molecular_3Group_Plot
    )

    pw <- stats::pairwise.wilcox.test(
        x = analysis_df[[value_col]],
        g = analysis_df$Molecular_3Group_Plot,
        p.adjust.method = "BH",
        exact = FALSE
    )

    pw_long <- as.data.frame(
        as.table(
            pw$p.value
        ),
        stringsAsFactors = FALSE
    )

    names(pw_long) <- c(
        "Group_1",
        "Group_2",
        "BH_Adjusted_P"
    )

    pw_long <- pw_long[
        !is.na(pw_long$BH_Adjusted_P),
        ,
        drop = FALSE
    ]

    group_summary <- analysis_df[
        ,
        .(
            N = .N,
            Median = median(
                get(value_col),
                na.rm = TRUE
            ),
            Q1 = quantile(
                get(value_col),
                0.25,
                na.rm = TRUE
            ),
            Q3 = quantile(
                get(value_col),
                0.75,
                na.rm = TRUE
            ),
            Mean = mean(
                get(value_col),
                na.rm = TRUE
            ),
            SD = sd(
                get(value_col),
                na.rm = TRUE
            )
        ),
        by = Molecular_3Group_Plot
    ]

    overall <- data.frame(
        Variable = variable_label,
        Test = "Kruskal-Wallis",
        Statistic = unname(kw$statistic),
        DF = unname(kw$parameter),
        P_Value = kw$p.value,
        stringsAsFactors = FALSE
    )

    list(
        overall = overall,
        pairwise = pw_long,
        summary = group_summary
    )
}


run_correlation <- function(
        data,
        x_col,
        y_col,
        comparison_label
) {

    analysis_df <- data[
        !is.na(get(x_col)) &
            !is.na(get(y_col))
    ]

    test <- stats::cor.test(
        analysis_df[[x_col]],
        analysis_df[[y_col]],
        method = "pearson"
    )

    data.frame(
        Comparison = comparison_label,
        N = nrow(analysis_df),
        Pearson_R = unname(test$estimate),
        P_Value = test$p.value,
        CI_Lower = test$conf.int[1],
        CI_Upper = test$conf.int[2],
        stringsAsFactors = FALSE
    )
}


make_boxplot <- function(
        data,
        value_col,
        y_label,
        panel_letter,
        kw_p
) {

    plot_df <- data[
        !is.na(Molecular_3Group_Plot) &
            !is.na(get(value_col))
    ]

    # Count samples actually used in this panel.
    group_counts <- plot_df[
        ,
        .N,
        by = Molecular_3Group_Plot
    ]

    # Build compact, publication-friendly x-axis labels.
    axis_labels <- setNames(
        paste0(
            c(
                "IDHmut-codel",
                "IDHmut-\nnoncodel",
                "IDHwt"
            ),
            "\n(n = ",
            group_counts$N[
                match(
                    c(
                        "IDHmut-codel",
                        "IDHmut-noncodel",
                        "IDHwt"
                    ),
                    as.character(
                        group_counts$Molecular_3Group_Plot
                    )
                )
            ],
            ")"
        ),
        c(
            "IDHmut-codel",
            "IDHmut-noncodel",
            "IDHwt"
        )
    )

    ggplot(
        plot_df,
        aes(
            x = Molecular_3Group_Plot,
            y = .data[[value_col]],
            fill = Molecular_3Group_Plot
        )
    ) +
        geom_violin(
            trim = FALSE,
            scale = "width",
            linewidth = 0.30,
            color = "black",
            alpha = 0.52
        ) +
        geom_boxplot(
            width = 0.17,
            outlier.shape = NA,
            linewidth = 0.34,
            fill = "white",
            color = "black",
            alpha = 0.88
        ) +
        geom_jitter(
            aes(
                color = Molecular_3Group_Plot
            ),
            width = 0.070,
            height = 0,
            size = JITTER_SIZE,
            alpha = 0.62,
            shape = 16,
            show.legend = FALSE
        ) +
        scale_fill_manual(
            values = GROUP_COLORS[c(
                "IDHmut-codel",
                "IDHmut-noncodel",
                "IDHwt"
            )],
            drop = FALSE,
            guide = "none"
        ) +
        scale_color_manual(
            values = GROUP_COLORS[c(
                "IDHmut-codel",
                "IDHmut-noncodel",
                "IDHwt"
            )],
            drop = FALSE,
            guide = "none"
        ) +
        scale_x_discrete(
            labels = axis_labels
        ) +
        labs(
            x = NULL,
            y = y_label,
            subtitle = paste0(
                "Kruskal-Wallis ",
                format_p(kw_p)
            ),
            tag = panel_letter
        ) +
        publication_theme() +
        theme(
            axis.text.x = element_text(
                size = 5.0,
                angle = 0,
                lineheight = 0.90,
                margin = margin(
                    t = 1.0,
                    unit = "mm"
                ),
                hjust = 0.5,
                vjust = 1
            ),
            plot.subtitle = element_text(
                family = FONT_FAMILY,
                size = FONT_STAT,
                hjust = 0.5,
                margin = margin(
                    b = 1.6,
                    unit = "mm"
                )
            ),
            legend.position = "none"
        )
}

make_scatter <- function(
        data,
        x_col,
        y_col,
        x_label,
        y_label,
        panel_letter,
        correlation_result
) {

    plot_df <- data[
        !is.na(get(x_col)) &
            !is.na(get(y_col))
    ]

    r_value <- correlation_result$Pearson_R[1]
    p_value <- correlation_result$P_Value[1]
    n_value <- correlation_result$N[1]

    annotation_text <- paste0(
        "Pearson r = ",
        sprintf("%.3f", r_value),
        "\n",
        format_p(p_value),
        ", n = ",
        n_value
    )

    ggplot(
        plot_df,
        aes(
            x = .data[[x_col]],
            y = .data[[y_col]]
        )
    ) +
        geom_point(
            aes(
                color = Molecular_3Group_Scatter
            ),
            size = POINT_SIZE,
            alpha = 0.64
        ) +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = TRUE,
            linewidth = 0.42,
            color = "black",
            fill = "grey80",
            alpha = 0.24
        ) +
        scale_color_manual(
            values = GROUP_COLORS,
            drop = FALSE,
            name = NULL
        ) +
        annotate(
            "text",
            x = -Inf,
            y = Inf,
            label = annotation_text,
            hjust = -0.04,
            vjust = 1.08,
            family = FONT_FAMILY,
            size = FONT_STAT / 2.845,
            lineheight = 0.95
        ) +
        coord_cartesian(
            clip = "off"
        ) +
        labs(
            x = x_label,
            y = y_label,
            tag = panel_letter
        ) +
        publication_theme() +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.key.width = grid::unit(
                3.0,
                "mm"
            ),
            legend.spacing.x = grid::unit(
                1.0,
                "mm"
            ),
            legend.margin = margin(
                t = 0.8,
                r = 0,
                b = 0,
                l = 0,
                unit = "mm"
            ),
            plot.margin = margin(
                4.0,
                5.0,
                6.0,
                5.0,
                unit = "mm"
            )
        ) +
        guides(
            color = guide_legend(
                override.aes = list(
                    size = 1.6,
                    alpha = 1
                ),
                nrow = 1,
                byrow = TRUE
            )
        )
}

# -----------------------------------------------------------------------------
# 8. Robust file saving functions
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

        message(
            "\nRetrying PDF save using a temporary short path:\n",
            fallback_file
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

            ggsave(
                filename = filename,
                plot = plot_object,
                width = width_mm,
                height = height_mm,
                units = "mm",
                dpi = FIGURE_DPI,
                bg = "white"
            )

        },
        error = function(e) {

            fallback_file <- file.path(
                tempdir(),
                basename(filename)
            )

            message(
                "\nPNG save failed at the selected output path.\n",
                "Retrying with a temporary short path:\n",
                fallback_file
            )

            ggsave(
                filename = fallback_file,
                plot = plot_object,
                width = width_mm,
                height = height_mm,
                units = "mm",
                dpi = FIGURE_DPI,
                bg = "white"
            )

            message(
                "\nFallback PNG saved successfully:\n",
                fallback_file
            )
        }
    )
}


save_panel <- function(
        plot_object,
        file_stub,
        width_mm,
        height_mm
) {

    pdf_file <- file.path(
        output_dir,
        paste0(
            file_stub,
            ".pdf"
        )
    )

    png_file <- file.path(
        output_dir,
        paste0(
            file_stub,
            ".png"
        )
    )

    safe_save_pdf(
        plot_object,
        pdf_file,
        width_mm,
        height_mm
    )

    safe_save_png(
        plot_object,
        png_file,
        width_mm,
        height_mm
    )
}


# -----------------------------------------------------------------------------
# 9. Statistical analyses
# -----------------------------------------------------------------------------

stats_expression <- run_group_statistics(
    df,
    expression_col,
    "DNMT3A mRNA expression"
)

stats_cg21629895 <- run_group_statistics(
    df,
    single_cpg_col,
    "cg21629895 methylation"
)

stats_promoter_mean <- run_group_statistics(
    df,
    promoter_mean_col,
    "Mean methylation across 10 DNMT3A promoter-associated CpGs"
)

cor_single <- run_correlation(
    df,
    single_cpg_col,
    expression_col,
    "cg21629895 methylation vs DNMT3A mRNA expression"
)

cor_promoter_mean <- run_correlation(
    df,
    promoter_mean_col,
    expression_col,
    "11-CpG promoter-associated mean methylation vs DNMT3A mRNA expression"
)


# -----------------------------------------------------------------------------
# 10. Generate Figure 1 panels
# -----------------------------------------------------------------------------

panel_A <- make_boxplot(
    df,
    expression_col,
    "DNMT3A expression\n(log2(norm_count + 1))",
    "A",
    stats_expression$overall$P_Value
)

panel_B <- make_boxplot(
    df,
    single_cpg_col,
    "cg21629895 methylation\n(beta value)",
    "B",
    stats_cg21629895$overall$P_Value
)

panel_C <- make_boxplot(
    df,
    promoter_mean_col,
    "11-CpG promoter mean methylation\n(beta value)",
    "C",
    stats_promoter_mean$overall$P_Value
)

panel_D <- make_scatter(
    df,
    single_cpg_col,
    expression_col,
    "cg21629895 methylation\n(beta value)",
    "DNMT3A expression\n(log2(norm_count + 1))",
    "D",
    cor_single
)

panel_E <- make_scatter(
    df,
    promoter_mean_col,
    expression_col,
    "11-CpG promoter mean methylation\n(beta value)",
    "DNMT3A expression\n(log2(norm_count + 1))",
    "E",
    cor_promoter_mean
)


# -----------------------------------------------------------------------------
# 11. Save individual panels
#
# These dimensions are chosen so that the combined final figure is 180 mm wide.
# The combined figure is the publication reference size.
# -----------------------------------------------------------------------------

save_panel(
    panel_A,
    "Fig1A",
    58,
    62
)

save_panel(
    panel_B,
    "Fig1B",
    58,
    62
)

save_panel(
    panel_C,
    "Fig1C",
    58,
    62
)

save_panel(
    panel_D,
    "Fig1D",
    88,
    66
)

save_panel(
    panel_E,
    "Fig1E",
    88,
    66
)


# -----------------------------------------------------------------------------
# 12. Assemble combined 180-mm Main Figure 1
# -----------------------------------------------------------------------------

top_row <- panel_A + panel_B + panel_C +
    plot_layout(
        widths = c(
            1,
            1,
            1
        )
    )

bottom_row <- panel_D + panel_E +
    plot_layout(
        widths = c(
            1,
            1
        )
    )

main_figure_1 <- top_row / bottom_row +
    plot_layout(
        heights = c(
            0.94,
            1.12
        ),
        guides = "collect"
    ) &
    theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center"
    )


combined_pdf <- file.path(
    output_dir,
    "Fig1_Combined_180mm_11Probe_Legacy_REFINED.pdf"
)

combined_png <- file.path(
    output_dir,
    "Fig1_Combined_180mm_11Probe_Legacy_REFINED_600dpi.png"
)


safe_save_pdf(
    main_figure_1,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

safe_save_png(
    main_figure_1,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)


# -----------------------------------------------------------------------------
# 13. Save statistics
# -----------------------------------------------------------------------------

overall_group_tests <- data.table::rbindlist(
    list(
        stats_expression$overall,
        stats_cg21629895$overall,
        stats_promoter_mean$overall
    ),
    fill = TRUE
)

data.table::fwrite(
    overall_group_tests,
    file.path(
        output_dir,
        "Fig1_Kruskal_Wallis.csv"
    )
)


pairwise_results <- data.table::rbindlist(
    list(
        cbind(
            Variable = "DNMT3A mRNA expression",
            stats_expression$pairwise
        ),
        cbind(
            Variable = "cg21629895 methylation",
            stats_cg21629895$pairwise
        ),
        cbind(
            Variable =
                "Mean methylation across 10 DNMT3A promoter-associated CpGs",
            stats_promoter_mean$pairwise
        )
    ),
    fill = TRUE
)

data.table::fwrite(
    pairwise_results,
    file.path(
        output_dir,
        "Fig1_Pairwise_Wilcoxon_BH.csv"
    )
)


group_summary_results <- data.table::rbindlist(
    list(
        cbind(
            Variable = "DNMT3A mRNA expression",
            stats_expression$summary
        ),
        cbind(
            Variable = "cg21629895 methylation",
            stats_cg21629895$summary
        ),
        cbind(
            Variable =
                "Mean methylation across 10 DNMT3A promoter-associated CpGs",
            stats_promoter_mean$summary
        )
    ),
    fill = TRUE
)

data.table::fwrite(
    group_summary_results,
    file.path(
        output_dir,
        "Fig1_Group_Summary.csv"
    )
)


correlation_results <- data.table::rbindlist(
    list(
        cor_single,
        cor_promoter_mean
    ),
    fill = TRUE
)

data.table::fwrite(
    correlation_results,
    file.path(
        output_dir,
        "Fig1_Pearson_Correlation.csv"
    )
)


# -----------------------------------------------------------------------------
# 14. QC summary
# -----------------------------------------------------------------------------

qc_summary <- data.frame(
    Metric = c(
        "Total rows",
        "Samples with molecular 3-group",
        "DNMT3A expression available",
        "cg21629895 available",
        "11-CpG promoter mean available",
        "Complete cg21629895-expression pairs",
        "Complete promoter mean-expression pairs",
        "Final combined figure width (mm)",
        "Final combined figure height (mm)",
        "Axis text size (pt)",
        "Axis title size (pt)",
        "Panel tag size (pt)"
    ),
    Value = c(
        nrow(df),
        sum(!is.na(df$Molecular_3Group_Plot)),
        sum(!is.na(df[[expression_col]])),
        sum(!is.na(df[[single_cpg_col]])),
        sum(!is.na(df[[promoter_mean_col]])),
        cor_single$N,
        cor_promoter_mean$N,
        FIGURE_WIDTH_MM,
        FIGURE_HEIGHT_MM,
        FONT_AXIS_TEXT,
        FONT_AXIS_TITLE,
        FONT_PANEL_TAG
    ),
    stringsAsFactors = FALSE
)

data.table::fwrite(
    qc_summary,
    file.path(
        output_dir,
        "Fig1_QC.csv"
    )
)


# -----------------------------------------------------------------------------
# 15. Final console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MAIN FIGURE 1 COMPLETED\n")
cat("============================================================\n")

cat("\nFinal publication dimensions:\n")
cat("Width:", FIGURE_WIDTH_MM, "mm\n")
cat("Height:", FIGURE_HEIGHT_MM, "mm\n")

cat("\nFinal font sizes:\n")
cat("Axis text:", FONT_AXIS_TEXT, "pt\n")
cat("Axis title:", FONT_AXIS_TITLE, "pt\n")
cat("Statistics:", FONT_STAT, "pt\n")
cat("Panel labels:", FONT_PANEL_TAG, "pt\n")

cat("\nPDF device:\n")
cat("grDevices::pdf (cairo_pdf is not used)\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nCombined figure:\n")
cat(combined_pdf, "\n")
cat(combined_png, "\n")

cat("\nKruskal-Wallis results:\n")
print(overall_group_tests)

cat("\nPearson correlation results:\n")
print(correlation_results)

cat("\n============================================================\n")
