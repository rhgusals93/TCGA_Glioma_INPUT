
# =============================================================================
# TCGA-LGG DNMT3A Supplementary Figure 2
# OS Kaplan-Meier analysis of promoter methylation and 10 individual CpG probes
#
# Panel A:
#   DNMT3A 11-CpG promoter mean methylation
#   -> Uses the EXISTING master-file median group:
#      DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group
#
# Panels B-K:
#   The 10 promoter CpGs other than cg21629895:
#     cg15998962
#     cg03463641
#     cg02208653
#     cg21708767
#     cg13344237
#     cg06112956
#     cg13076778
#     cg26470599
#     cg22731525
#     cg07150430
#
# For the 10 individual CpGs:
#   High/Low groups are generated using the median beta value of each CpG:
#      Low  = beta < median
#      High = beta >= median
#
# Survival endpoint:
#   Overall survival (OS) only
#   Source columns: OS_months, OS_status_raw
#
# Statistics:
#   - Kaplan-Meier
#   - Log-rank test
#   - Number-at-risk table
#   - Median cutoffs saved for every CpG
#
# Annotation:
#   The user selects the HM450 annotation CSV manually.
#   Plot titles include:
#      CpG ID
#      UCSC_RefGene_Group, Relation_to_Island
#
# Input and output:
#   1) Master CSV selected manually with file.choose()
#   2) HM450 annotation CSV selected manually with file.choose()
#   3) Output folder selected manually
#
# Outputs:
#   - Combined Supplementary Figure 2 PDF
#   - Combined 600-dpi PNG
#   - Individual panel PDFs/PNGs
#   - KM statistics CSV
#   - Cutoff/group QC CSV
#   - Excel workbook
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install/load packages
# -----------------------------------------------------------------------------

required_packages <- c(
    "data.table",
    "ggplot2",
    "patchwork",
    "survival",
    "cowplot",
    "grid",
    "openxlsx"
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
    library(survival)
    library(cowplot)
    library(grid)
    library(openxlsx)
})


# -----------------------------------------------------------------------------
# 1. Publication settings
# -----------------------------------------------------------------------------

set.seed(20260914)

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 215
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"

FONT_PANEL_TAG <- 7.0
FONT_TITLE <- 5.8
FONT_AXIS_TITLE <- 5.2
FONT_AXIS_TEXT <- 5.0
FONT_TABLE <- 4.8
FONT_STAT <- 4.9

INK <- "#202124"
MUTED <- "#6B6F76"

BINARY_COLORS <- c(
    "Low" = "#0072B2",
    "High" = "#D55E00"
)


# -----------------------------------------------------------------------------
# 2. Select master file manually
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE 2: PROMOTER/CpG OS SURVIVAL\n")
cat("============================================================\n")

cat("\n1. Select the completed TCGA-LGG DNMT3A master CSV file.\n")

MASTER_FILE <- file.choose()

if (!file.exists(MASTER_FILE)) {
    stop("Master-file selection failed.")
}

df <- data.table::fread(
    MASTER_FILE,
    check.names = FALSE,
    data.table = TRUE
)

cat("\nMaster file loaded.\n")
cat("Rows: ", nrow(df), "\n", sep = "")
cat("Columns: ", ncol(df), "\n", sep = "")


# -----------------------------------------------------------------------------
# 3. Select HM450 annotation file manually
# -----------------------------------------------------------------------------

cat(
    "\n2. Select the DNMT3A HM450 probe annotation CSV file.\n",
    "Required columns: Name, UCSC_RefGene_Group, Relation_to_Island\n",
    sep = ""
)

ANNOTATION_FILE <- file.choose()

if (!file.exists(ANNOTATION_FILE)) {
    stop("Annotation-file selection failed.")
}

annotation_df <- data.table::fread(
    ANNOTATION_FILE,
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


# -----------------------------------------------------------------------------
# 4. Select output folder manually
# -----------------------------------------------------------------------------

cat("\n3. Select the folder where Supplementary Figure 2 will be saved.\n")

if (.Platform$OS.type == "windows") {

    selected_output <- utils::choose.dir(
        default = Sys.getenv("USERPROFILE"),
        caption = "Choose output folder for Supplementary Figure 2"
    )

    if (
        is.na(selected_output) ||
        !nzchar(selected_output)
    ) {
        stop("Output-folder selection was cancelled.")
    }

} else {

    cat("Enter the full output-folder path and press Enter:\n")

    selected_output <- readline()

    if (
        !nzchar(selected_output) ||
        !dir.exists(selected_output)
    ) {
        stop("A valid output folder was not selected.")
    }
}

# Use the folder selected by the user DIRECTLY.
# Do not append another "DNMT3A_SuppFig2_Survival" folder.
OUTPUT_DIR <- selected_output

if (!dir.exists(OUTPUT_DIR)) {
    stop(
        paste0(
            "Selected output directory does not exist: ",
            OUTPUT_DIR
        )
    )
}

OUTPUT_DIR <- normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = TRUE
)

# Windows' PDF device can fail when the complete output path becomes too long.
# Detect this BEFORE plotting and let the user choose a shorter location.
if (
    .Platform$OS.type == "windows" &&
    nchar(OUTPUT_DIR) > 170
) {

    cat(
        "\nThe selected output path is too long for reliable PDF saving in Windows:\n",
        OUTPUT_DIR,
        "\n\nPlease select a SHORT folder, preferably directly on Desktop,\n",
        "for example: C:/Users/hko/Desktop/SFig2\n\n",
        sep = ""
    )

    selected_output_short <- utils::choose.dir(
        default = file.path(
            Sys.getenv("USERPROFILE"),
            "Desktop"
        ),
        caption = "Choose a SHORT output folder for Supplementary Figure 2"
    )

    if (
        is.na(selected_output_short) ||
        !nzchar(selected_output_short)
    ) {
        stop("Short output-folder selection was cancelled.")
    }

    OUTPUT_DIR <- normalizePath(
        selected_output_short,
        winslash = "/",
        mustWork = TRUE
    )
}

cat("\nFinal output directory:\n")
cat(OUTPUT_DIR, "\n")
cat("Path length: ", nchar(OUTPUT_DIR), " characters\n", sep = "")


# -----------------------------------------------------------------------------
# 5. Define variables
# -----------------------------------------------------------------------------

promoter_mean_col <- "DNMT3A_Promoter_Mean_Methylation_11Probe"
promoter_group_col <- "DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group"

main_figure_cpg <- "cg21629895"

supplementary_probes_10 <- c(
    "cg15998962",
    "cg03463641",
    "cg02208653",
    "cg21708767",
    "cg13344237",
    "cg06112956",
    "cg13076778",
    "cg26470599",
    "cg22731525",
    "cg07150430"
)

required_master_cols <- c(
    "Patient_ID",
    "OS_months",
    "OS_status_raw",
    promoter_mean_col,
    promoter_group_col,
    supplementary_probes_10
)

missing_master_cols <- setdiff(
    required_master_cols,
    names(df)
)

if (length(missing_master_cols) > 0) {
    stop(
        paste0(
            "Master file is missing required columns:\n",
            paste(
                missing_master_cols,
                collapse = "\n"
            )
        )
    )
}


# -----------------------------------------------------------------------------
# 6. Utility functions
# -----------------------------------------------------------------------------

safe_numeric <- function(x) {
    suppressWarnings(
        as.numeric(
            as.character(x)
        )
    )
}

parse_event_status <- function(x) {

    x <- toupper(
        trimws(
            as.character(x)
        )
    )

    data.table::fcase(
        is.na(x) | x == "", NA_integer_,
        grepl("^1:", x), 1L,
        grepl("^0:", x), 0L,
        x == "1", 1L,
        x == "0", 0L,
        default = NA_integer_
    )
}

pt_to_mm <- function(pt) {
    pt * 0.3527778
}

format_p <- function(p) {

    if (is.na(p)) {
        return("NA")
    }

    if (p < 0.001) {
        return("<0.001")
    }

    sprintf(
        "%.3f",
        p
    )
}

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


# -----------------------------------------------------------------------------
# 7. Prepare OS endpoint
# -----------------------------------------------------------------------------

df[
    ,
    OS_time := safe_numeric(
        OS_months
    )
]

df[
    ,
    OS_event := parse_event_status(
        OS_status_raw
    )
]

bad_os <- (
    !is.na(df$OS_time) &
        df$OS_time < 0
)

if (any(bad_os)) {
    df$OS_time[bad_os] <- NA_real_
    df$OS_event[bad_os] <- NA_integer_
}

cat("\nOS QC:\n")
cat(
    "Evaluable samples: ",
    sum(
        complete.cases(
            df$OS_time,
            df$OS_event
        )
    ),
    "\n",
    sep = ""
)
cat(
    "Events: ",
    sum(
        df$OS_event == 1,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 8. Prepare methylation variables
# -----------------------------------------------------------------------------

df[
    ,
    (promoter_mean_col) := safe_numeric(
        get(
            promoter_mean_col
        )
    )
]

for (probe in supplementary_probes_10) {

    df[
        ,
        (probe) := safe_numeric(
            get(
                probe
            )
        )
    ]
}


# -----------------------------------------------------------------------------
# 9. Promoter-mean High/Low group
#
# IMPORTANT:
# Use the EXISTING master-file median group exactly as provided.
# Do not recalculate this group.
# -----------------------------------------------------------------------------

df[
    ,
    Promoter11_Group_KM := factor(
        as.character(
            get(
                promoter_group_col
            )
        ),
        levels = c(
            "Low",
            "High"
        )
    )
]

promoter_median_for_report <- stats::median(
    df[[promoter_mean_col]],
    na.rm = TRUE
)


# -----------------------------------------------------------------------------
# 10. Create median High/Low groups for the 10 individual CpGs
# -----------------------------------------------------------------------------

probe_group_info <- data.table::rbindlist(
    lapply(
        supplementary_probes_10,
        function(probe) {

            cutoff <- stats::median(
                df[[probe]],
                na.rm = TRUE
            )

            group_col <- paste0(
                probe,
                "_Median_Group"
            )

            df[
                ,
                (group_col) := factor(
                    data.table::fcase(
                        is.na(
                            get(probe)
                        ),
                        NA_character_,

                        get(probe) < cutoff,
                        "Low",

                        get(probe) >= cutoff,
                        "High",

                        default = NA_character_
                    ),
                    levels = c(
                        "Low",
                        "High"
                    )
                )
            ]

            data.table::data.table(
                Variable = probe,
                Group_Source = "Calculated from individual CpG median",
                Median_Cutoff = cutoff,
                Low_N = sum(
                    df[[group_col]] == "Low",
                    na.rm = TRUE
                ),
                High_N = sum(
                    df[[group_col]] == "High",
                    na.rm = TRUE
                ),
                Missing_N = sum(
                    is.na(
                        df[[group_col]]
                    )
                )
            )
        }
    ),
    fill = TRUE
)

promoter_group_info <- data.table::data.table(
    Variable = "11-CpG promoter mean",
    Group_Source = promoter_group_col,
    Median_Cutoff = promoter_median_for_report,
    Low_N = sum(
        df$Promoter11_Group_KM == "Low",
        na.rm = TRUE
    ),
    High_N = sum(
        df$Promoter11_Group_KM == "High",
        na.rm = TRUE
    ),
    Missing_N = sum(
        is.na(
            df$Promoter11_Group_KM
        )
    )
)

group_qc <- data.table::rbindlist(
    list(
        promoter_group_info,
        probe_group_info
    ),
    fill = TRUE
)

cat("\nMedian-group QC:\n")
print(group_qc)


# -----------------------------------------------------------------------------
# 11. Prepare CpG annotation labels
# -----------------------------------------------------------------------------

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
            "The following probes were not found in the annotation file:\n",
            paste(
                missing_probe_annotations,
                collapse = "\n"
            )
        )
    )
}

probe_annotation[
    ,
    Plot_Title := paste0(
        CpG,
        "\n",
        UCSC_RefGene_Group,
        ", ",
        Relation_to_Island
    )
]


# -----------------------------------------------------------------------------
# 12. Define a COMMON OS x-axis for all panels
#
# A shared x-axis makes the 11 KM panels directly comparable.
# -----------------------------------------------------------------------------

os_complete <- df[
    !is.na(OS_time) &
        !is.na(OS_event)
]

COMMON_RISK_TIMES <- seq(
    0,
    ceiling(
        max(
            os_complete$OS_time,
            na.rm = TRUE
        ) / 48
    ) * 48,
    by = 48
)

COMMON_X_MAX <- max(
    COMMON_RISK_TIMES
) + 0.035 * max(
    COMMON_RISK_TIMES
)

cat("\nCommon risk-table times:\n")
print(COMMON_RISK_TIMES)


# -----------------------------------------------------------------------------
# 13. Publication-style Kaplan-Meier plot with risk table
#
# This is intentionally matched to the Main Figure 3 survival style:
#   - same colors
#   - same curve/censor widths
#   - same 48-month time points
#   - same risk-table geometry
#   - same font sizes
#   - same legend treatment
#
# The only Supplementary-Figure-specific change is that CpG titles can occupy
# two lines so the UCSC gene-region and CpG-island annotation can be shown.
# -----------------------------------------------------------------------------

make_km_with_risk_table <- function(
        data,
        time_col,
        event_col,
        group_col,
        panel_letter,
        panel_title,
        group_levels = c("Low", "High"),
        colors = BINARY_COLORS,
        show_confidence = TRUE
) {

    keep <- (
        !is.na(data[[time_col]]) &
        !is.na(data[[event_col]]) &
        !is.na(data[[group_col]])
    )

    km_data <- data.frame(
        Time = as.numeric(data[[time_col]][keep]),
        Event = as.numeric(data[[event_col]][keep]),
        KM_Group = factor(
            as.character(data[[group_col]][keep]),
            levels = group_levels
        ),
        stringsAsFactors = FALSE
    )

    km_data <- km_data[
        !is.na(km_data$KM_Group) &
        is.finite(km_data$Time) &
        km_data$Time >= 0 &
        km_data$Event %in% c(0, 1),
        ,
        drop = FALSE
    ]

    if (
        nrow(km_data) < 3 ||
        length(unique(km_data$KM_Group)) < 2
    ) {
        stop(
            paste0(
                "Insufficient evaluable samples for Kaplan-Meier analysis: ",
                group_col
            )
        )
    }

    fit <- survival::survfit(
        survival::Surv(Time, Event) ~ KM_Group,
        data = km_data,
        conf.type = "log"
    )

    surv_diff <- survival::survdiff(
        survival::Surv(Time, Event) ~ KM_Group,
        data = km_data
    )

    p_value <- stats::pchisq(
        surv_diff$chisq,
        df = length(surv_diff$n) - 1,
        lower.tail = FALSE
    )

    p_text <- paste0(
        "Log-rank P ",
        ifelse(
            p_value < 0.001,
            "< 0.001",
            paste0(
                "= ",
                sprintf("%.3f", p_value)
            )
        )
    )

    fit_summary <- summary(
        fit,
        censored = TRUE
    )

    strata_text <- sub(
        "^KM_Group=",
        "",
        as.character(fit_summary$strata)
    )

    curve_df <- data.frame(
        Time = fit_summary$time,
        Survival = fit_summary$surv,
        Lower = fit_summary$lower,
        Upper = fit_summary$upper,
        N_Censor = fit_summary$n.censor,
        Group = factor(
            strata_text,
            levels = group_levels
        ),
        stringsAsFactors = FALSE
    )

    zero_df <- data.frame(
        Time = 0,
        Survival = 1,
        Lower = 1,
        Upper = 1,
        N_Censor = 0,
        Group = factor(
            group_levels,
            levels = group_levels
        )
    )

    curve_df <- rbind(
        zero_df,
        curve_df
    )

    curve_df <- curve_df[
        order(
            curve_df$Group,
            curve_df$Time
        ),
        ,
        drop = FALSE
    ]

    censor_df <- curve_df[
        curve_df$N_Censor > 0,
        ,
        drop = FALSE
    ]

    group_counts <- table(
        factor(
            km_data$KM_Group,
            levels = group_levels
        )
    )

    final_legend_labels <- setNames(
        paste0(
            group_levels,
            " (n = ",
            as.integer(group_counts),
            ")"
        ),
        group_levels
    )

    # Compact common time points for the 3-column layout.
    # Four ticks preserve readability while keeping the full 0-240 month range.
    risk_times <- c(
        0,
        96,
        192,
        240
    )

    x_display_max <- 246

    risk_summary <- summary(
        fit,
        times = risk_times,
        extend = TRUE
    )

    risk_group <- sub(
        "^KM_Group=",
        "",
        as.character(risk_summary$strata)
    )

    risk_df <- data.frame(
        Time = risk_summary$time,
        N_Risk = risk_summary$n.risk,
        Group = factor(
            risk_group,
            levels = group_levels
        ),
        stringsAsFactors = FALSE
    )

    confidence_layer <- if (isTRUE(show_confidence)) {
        ggplot2::geom_ribbon(
            ggplot2::aes(
                ymin = pmax(Lower, 0),
                ymax = pmin(Upper, 1)
            ),
            alpha = 0.055,
            linewidth = 0,
            show.legend = FALSE
        )
    } else {
        NULL
    }

    survival_plot <- ggplot2::ggplot(
        curve_df,
        ggplot2::aes(
            x = Time,
            y = Survival,
            color = Group,
            fill = Group,
            group = Group
        )
    ) +
        confidence_layer +
        ggplot2::geom_step(
            linewidth = 0.64,
            direction = "hv"
        ) +
        ggplot2::geom_point(
            data = censor_df,
            ggplot2::aes(
                x = Time,
                y = Survival,
                color = Group
            ),
            inherit.aes = FALSE,
            shape = 3,
            size = 0.72,
            stroke = 0.32
        ) +
        ggplot2::scale_color_manual(
            values = colors[group_levels],
            breaks = group_levels,
            labels = unname(final_legend_labels[group_levels]),
            name = NULL,
            drop = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values = colors[group_levels],
            breaks = group_levels,
            guide = "none",
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(0, x_display_max),
            breaks = risk_times,
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1.02),
            breaks = c(0, 0.25, 0.50, 0.75, 1),
            expand = c(0, 0)
        ) +
        ggplot2::annotate(
            "label",
            x = 0.03 * x_display_max,
            y = 0.10,
            label = p_text,
            hjust = 0,
            size = pt_to_mm(5.0),
            fontface = "bold",
            family = FONT_FAMILY,
            color = INK,
            fill = "white",
            label.padding = grid::unit(
                0.10,
                "lines"
            )
        ) +
        ggplot2::labs(
            title = panel_title,
            x = "Overall survival (months)",
            y = "Survival probability"
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(
                nrow = 1,
                byrow = TRUE
            )
        ) +
        ggplot2::theme_classic(
            base_size = 7,
            base_family = FONT_FAMILY
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = FONT_TITLE,
                hjust = 0,
                color = INK,
                lineheight = 0.90,
                margin = ggplot2::margin(
                    t = 0.4,
                    b = 1.0
                )
            ),
            axis.title = ggplot2::element_text(
                size = FONT_AXIS_TITLE,
                face = "bold",
                color = INK
            ),
            axis.text = ggplot2::element_text(
                size = FONT_AXIS_TEXT,
                color = INK
            ),
            axis.line = ggplot2::element_line(
                linewidth = 0.38,
                color = INK
            ),
            axis.ticks = ggplot2::element_line(
                linewidth = 0.32,
                color = INK
            ),
            legend.position = "top",
            legend.justification = "center",
            legend.direction = "horizontal",
            legend.text = ggplot2::element_text(
                size = FONT_TABLE,
                color = INK
            ),
            legend.key.width = grid::unit(
                3.0,
                "mm"
            ),
            legend.key.height = grid::unit(
                2.0,
                "mm"
            ),
            plot.margin = ggplot2::margin(
                t = 0.6,
                r = 3.0,
                b = 0.3,
                l = 3.0
            )
        )

    risk_plot <- ggplot2::ggplot(
        risk_df,
        ggplot2::aes(
            x = Time,
            y = Group,
            label = N_Risk,
            color = Group
        )
    ) +
        ggplot2::geom_text(
            size = pt_to_mm(FONT_TABLE),
            show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
            values = colors[group_levels],
            guide = "none",
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(
                -0.12 * max(risk_times),
                x_display_max
            ),
            breaks = risk_times,
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_discrete(
            limits = rev(group_levels),
            labels = rev(group_levels)
        ) +
        ggplot2::labs(
            title = "Number at risk",
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_classic(
            base_size = 6,
            base_family = FONT_FAMILY
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = FONT_TABLE,
                hjust = 0,
                color = INK,
                margin = ggplot2::margin(
                    b = 1
                )
            ),
            axis.text.x = ggplot2::element_text(
                size = FONT_TABLE,
                color = INK
            ),
            axis.text.y = ggplot2::element_text(
                size = 5.0,
                color = INK,
                margin = ggplot2::margin(
                    r = 2.0
                )
            ),
            axis.ticks = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(
                linewidth = 0.28,
                color = "#A7ABB1"
            ),
            plot.margin = ggplot2::margin(
                t = 0,
                r = 3.0,
                b = 0.5,
                l = 4.0
            )
        )

    combined_plot <- cowplot::plot_grid(
        survival_plot,
        risk_plot,
        ncol = 1,
        rel_heights = c(
            0.81,
            0.19
        ),
        align = "v",
        axis = "lr"
    )

    combined_plot <- cowplot::ggdraw(
        combined_plot
    ) +
        cowplot::draw_label(
            panel_letter,
            x = 0.003,
            y = 0.997,
            fontface = "bold",
            size = FONT_PANEL_TAG,
            color = INK,
            hjust = 0,
            vjust = 1
        )

    list(
        plot = combined_plot,
        p_value = p_value,
        n = nrow(km_data),
        events = sum(
            km_data$Event == 1,
            na.rm = TRUE
        ),
        fit = fit,
        risk_table = risk_df,
        group_counts = group_counts
    )
}


# -----------------------------------------------------------------------------
# 14. Build promoter-mean panel A
# -----------------------------------------------------------------------------

panel_results <- list()

panel_results[["PromoterMean"]] <- make_km_with_risk_table(
    data = df,
    time_col = "OS_time",
    event_col = "OS_event",
    group_col = "Promoter11_Group_KM",
    panel_letter = "A",
    panel_title = "11-CpG promoter mean methylation",
    show_confidence = TRUE
)


# -----------------------------------------------------------------------------
# 15. Build the 10 individual CpG panels B-K
# -----------------------------------------------------------------------------

probe_panel_letters <- LETTERS[
    2:11
]

for (i in seq_along(supplementary_probes_10)) {

    probe <- supplementary_probes_10[i]

    group_col <- paste0(
        probe,
        "_Median_Group"
    )

    title_row <- probe_annotation[
        CpG == probe
    ]

    panel_results[[probe]] <- make_km_with_risk_table(
        data = df,
        time_col = "OS_time",
        event_col = "OS_event",
        group_col = group_col,
        panel_letter = probe_panel_letters[i],
        panel_title = title_row$Plot_Title[1],
        show_confidence = TRUE
    )
}


# -----------------------------------------------------------------------------
# 16. Assemble Supplementary Figure 2
#
# FINAL LAYOUT:
#   180 mm x 215 mm
#   3 columns x 4 rows
#
# This is a compromise between readability and total figure height:
#   - each panel is ~60 mm wide
#   - only four rows, so panels retain adequate vertical height
#   - all rows have exactly the same height
#   - all inter-row gaps are identical
#
# Row 1: A B C
# Row 2: D E F
# Row 3: G H I
# Row 4: J K [blank]
# -----------------------------------------------------------------------------

plot_list <- c(
    list(
        panel_results[["PromoterMean"]]$plot
    ),
    lapply(
        supplementary_probes_10,
        function(probe) {
            panel_results[[probe]]$plot
        }
    )
)

row_1 <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
    patchwork::plot_layout(
        widths = c(
            1,
            1,
            1
        )
    )

row_2 <- plot_list[[4]] + plot_list[[5]] + plot_list[[6]] +
    patchwork::plot_layout(
        widths = c(
            1,
            1,
            1
        )
    )

row_3 <- plot_list[[7]] + plot_list[[8]] + plot_list[[9]] +
    patchwork::plot_layout(
        widths = c(
            1,
            1,
            1
        )
    )

row_4 <- plot_list[[10]] + plot_list[[11]] + patchwork::plot_spacer() +
    patchwork::plot_layout(
        widths = c(
            1,
            1,
            1
        )
    )

vertical_gap <- patchwork::plot_spacer()

supplementary_figure_2 <- row_1 / vertical_gap /
                          row_2 / vertical_gap /
                          row_3 / vertical_gap /
                          row_4 +
    patchwork::plot_layout(
        heights = c(
            1,
            0.035,
            1,
            0.035,
            1,
            0.035,
            1
        )
    )


# -----------------------------------------------------------------------------
# 17. Save helpers
#
# Windows-safe:
#   - verify output folder exists
#   - use short filenames
#   - always close graphics device on error
# -----------------------------------------------------------------------------

save_pdf <- function(
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

    filename <- normalizePath(
        out_dir,
        winslash = "/",
        mustWork = TRUE
    ) |>
        file.path(
            basename(filename)
        )

    if (
        .Platform$OS.type == "windows" &&
        nchar(filename) > 240
    ) {
        stop(
            paste0(
                "PDF output path is still too long (",
                nchar(filename),
                " characters).\n",
                "Please rerun and choose a shorter output folder, such as:\n",
                "C:/Users/hko/Desktop/SFig2\n\n",
                "Current path:\n",
                filename
            )
        )
    }

    pdf_open <- FALSE

    tryCatch(
        {

            grDevices::pdf(
                file = filename,
                width = width_mm / 25.4,
                height = height_mm / 25.4,
                family = "sans",
                useDingbats = FALSE
            )

            pdf_open <- TRUE

            print(
                plot_object
            )

            grDevices::dev.off()
            pdf_open <- FALSE

        },
        error = function(e) {

            if (pdf_open && grDevices::dev.cur() > 1) {
                try(
                    grDevices::dev.off(),
                    silent = TRUE
                )
            }

            stop(
                paste0(
                    "PDF save failed.\n\n",
                    "Path:\n",
                    filename,
                    "\n\nError:\n",
                    conditionMessage(e)
                )
            )
        }
    )
}

save_png <- function(
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

    filename <- normalizePath(
        out_dir,
        winslash = "/",
        mustWork = TRUE
    ) |>
        file.path(
            basename(filename)
        )

    if (
        .Platform$OS.type == "windows" &&
        nchar(filename) > 240
    ) {
        stop(
            paste0(
                "PNG output path is too long.\n",
                "Please choose a shorter output folder.\n",
                filename
            )
        )
    }

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
# 18. Save combined Supplementary Figure 2
# -----------------------------------------------------------------------------

combined_pdf <- file.path(
    OUTPUT_DIR,
    "SuppFig2_KM_3col_FINAL.pdf"
)

combined_png <- file.path(
    OUTPUT_DIR,
    "SuppFig2_KM_3col_FINAL_600dpi.png"
)

save_pdf(
    supplementary_figure_2,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

save_png(
    supplementary_figure_2,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)


# -----------------------------------------------------------------------------
# 19. Save individual panels
# -----------------------------------------------------------------------------

individual_dir <- file.path(
    OUTPUT_DIR,
    "Panels"
)

dir.create(
    individual_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# Promoter mean.
save_pdf(
    panel_results[["PromoterMean"]]$plot,
    file.path(
        individual_dir,
        "A_PromoterMean.pdf"
    ),
    90,
    65
)

save_png(
    panel_results[["PromoterMean"]]$plot,
    file.path(
        individual_dir,
        "A_PromoterMean_600dpi.png"
    ),
    90,
    65
)

# Individual CpGs.
for (i in seq_along(supplementary_probes_10)) {

    probe <- supplementary_probes_10[i]
    letter <- probe_panel_letters[i]

    save_pdf(
        panel_results[[probe]]$plot,
        file.path(
            individual_dir,
            paste0(
                "",
                letter,
                "_",
                probe,
                ".pdf"
            )
        ),
        88,
        62
    )

    save_png(
        panel_results[[probe]]$plot,
        file.path(
            individual_dir,
            paste0(
                "",
                letter,
                "_",
                probe,
                "_600dpi.png"
            )
        ),
        88,
        62
    )
}


# -----------------------------------------------------------------------------
# 20. Build KM statistics table
# -----------------------------------------------------------------------------

km_statistics <- data.table::rbindlist(
    c(
        list(
            data.table::data.table(
                Panel = "A",
                Variable = "11-CpG promoter mean methylation",
                CpG = NA_character_,
                UCSC_RefGene_Group = NA_character_,
                Relation_to_Island = NA_character_,
                N = panel_results[["PromoterMean"]]$n,
                Events = panel_results[["PromoterMean"]]$events,
                Low_N = as.integer(
                    panel_results[["PromoterMean"]]$group_counts["Low"]
                ),
                High_N = as.integer(
                    panel_results[["PromoterMean"]]$group_counts["High"]
                ),
                LogRank_P = panel_results[["PromoterMean"]]$p_value
            )
        ),
        lapply(
            seq_along(
                supplementary_probes_10
            ),
            function(i) {

                probe <- supplementary_probes_10[i]

                ann <- probe_annotation[
                    CpG == probe
                ]

                data.table::data.table(
                    Panel = probe_panel_letters[i],
                    Variable = "Individual CpG methylation",
                    CpG = probe,
                    UCSC_RefGene_Group = ann$UCSC_RefGene_Group[1],
                    Relation_to_Island = ann$Relation_to_Island[1],
                    N = panel_results[[probe]]$n,
                    Events = panel_results[[probe]]$events,
                    Low_N = as.integer(
                        panel_results[[probe]]$group_counts["Low"]
                    ),
                    High_N = as.integer(
                        panel_results[[probe]]$group_counts["High"]
                    ),
                    LogRank_P = panel_results[[probe]]$p_value
                )
            }
        )
    ),
    fill = TRUE
)

# BH correction across the 10 individual CpG log-rank tests.
km_statistics[
    !is.na(CpG),
    BH_FDR_10CpGs := stats::p.adjust(
        LogRank_P,
        method = "BH"
    )
]

km_statistics[
    is.na(CpG),
    BH_FDR_10CpGs := NA_real_
]

data.table::fwrite(
    km_statistics,
    file.path(
        OUTPUT_DIR,
        "SuppFig2_KM_Statistics.csv"
    ),
    na = ""
)


# -----------------------------------------------------------------------------
# 21. Save group/cutoff QC
# -----------------------------------------------------------------------------

group_qc_with_annotation <- merge(
    group_qc,
    probe_annotation[
        ,
        .(
            Variable = CpG,
            UCSC_RefGene_Group,
            Relation_to_Island
        )
    ],
    by = "Variable",
    all.x = TRUE,
    sort = FALSE
)

data.table::fwrite(
    group_qc_with_annotation,
    file.path(
        OUTPUT_DIR,
        "SuppFig2_MedianCutoff_QC.csv"
    ),
    na = ""
)


# -----------------------------------------------------------------------------
# 22. Save Excel workbook
# -----------------------------------------------------------------------------

excel_file <- file.path(
    OUTPUT_DIR,
    "SuppFig2_Survival_Results.xlsx"
)

wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    textDecoration = "bold",
    fgFill = "#E7E7E7",
    border = "Bottom",
    borderColour = "#000000",
    halign = "center"
)

openxlsx::addWorksheet(
    wb,
    "KM_Statistics",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "KM_Statistics",
    km_statistics,
    headerStyle = header_style
)

openxlsx::freezePane(
    wb,
    "KM_Statistics",
    firstRow = TRUE
)

openxlsx::setColWidths(
    wb,
    "KM_Statistics",
    cols = 1:ncol(km_statistics),
    widths = "auto"
)

openxlsx::addWorksheet(
    wb,
    "Median_Cutoff_QC",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Median_Cutoff_QC",
    group_qc_with_annotation,
    headerStyle = header_style
)

openxlsx::freezePane(
    wb,
    "Median_Cutoff_QC",
    firstRow = TRUE
)

openxlsx::setColWidths(
    wb,
    "Median_Cutoff_QC",
    cols = 1:ncol(group_qc_with_annotation),
    widths = "auto"
)

openxlsx::saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
)


# -----------------------------------------------------------------------------
# 23. Final console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE 2 COMPLETE\n")
cat("============================================================\n")

cat("\nPanel A:\n")
cat("  11-CpG promoter mean methylation\n")
cat("  Group source: EXISTING master-file median group\n")

cat("\nPanels B-K:\n")
cat("  Remaining 10 promoter CpGs\n")
cat("  Group definition: Low < probe median; High >= probe median\n")

cat("\nMain-Figure CpG intentionally excluded:\n")
cat("  ", main_figure_cpg, "\n", sep = "")


cat("\nFigure geometry:\n")
cat("  Combined: 180 mm x 215 mm\n")
cat("  Layout: 3 columns x 4 rows\n")
cat("  Uniform row heights and vertical spacing\n")
cat("  Four common time points: 0, 96, 192, 240 months\n")
cat("  Main Figure 3 survival styling retained\n")

cat("\nSurvival endpoint:\n")
cat("  Overall survival (OS)\n")

cat("\nMultiple-testing correction:\n")
cat("  BH-FDR across the 10 individual CpG log-rank tests\n")

cat("\nCombined figure:\n")
cat("  ", combined_pdf, "\n", sep = "")
cat("  ", combined_png, "\n", sep = "")

cat("\nExcel workbook:\n")
cat("  ", excel_file, "\n", sep = "")

cat("\nKM statistics:\n")
print(
    km_statistics
)

cat("\n============================================================\n")
