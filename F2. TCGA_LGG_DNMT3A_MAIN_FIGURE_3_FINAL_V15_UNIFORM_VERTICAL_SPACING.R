# =============================================================================
# TCGA-LGG DNMT3A Main Figure 3
# Prognostic Significance
# FINAL V12 - Panel G removed; Panel H risk-table labels compacted to match A/B axis geometry
#
# Panels
# ------
# A. Integrated molecular states scatter plot
# B. Combined DNMT3A expression / cg21629895 OS Kaplan-Meier
# C. DNMT3A expression High vs Low OS Kaplan-Meier
# D. cg21629895 methylation Low vs High OS Kaplan-Meier
# E. Univariate OS Cox forest plot
# F. Multivariable OS Cox forest plot
# G. DNMT3A expression time-dependent OS ROC
# H. cg21629895 methylation time-dependent OS ROC
#
# Cox factors requested
# ---------------------
# Univariate:
#   Age (>45)
#   DNMT3A mRNA (High)
#   cg21629895 methylation (Low)
#   11-CpG promoter methylation (Low)
#   Gender (Male)
#   CDKN2A/B homozygous deletion
#   MGMT promoter status (Unmethylated)
#   EGFR amplification
#   DNMT3A expression (continuous, per SD)
#   cg21629895 methylation (continuous, per SD)
#
# Multivariable:
#   Age (>45)
#   DNMT3A mRNA (High)
#   cg21629895 methylation (Low)
#   Gender (Male)
#   CDKN2A/B homozygous deletion
#   MGMT promoter status (Unmethylated)
#
# Survival source columns
# -----------------------
# OS_months, OS_status_raw
# PFS_months, PFS_status_raw
# DSS_months, DSS_status_raw
#
# Cox regression is performed for OS ONLY.
#
# Current biomarker variables
# ---------------------------
# DNMT3A_expression
# DNMT3A_expression_Median_Group
# cg21629895
# cg21629895_Median_Group
# DNMT3A_Promoter_Mean_Methylation_11Probe
#
# Notes
# -----
# - Promoter methylation uses the current 11-CpG promoter mean.
# - Promoter High/Low is defined internally using the median 11-CpG mean.
# - Continuous expression and methylation are reported per 1-SD increase.
# - Time-dependent ROC is evaluated at 12, 24, and 36 months.
# - For cg21629895 ROC only, the marker sign is reversed so that lower
#   methylation corresponds to a higher risk-oriented marker value.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install/load packages
# -----------------------------------------------------------------------------

cran_packages <- c(
    "data.table",
    "ggplot2",
    "patchwork",
    "survival",
    "cowplot",
    "timeROC",
    "grid"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

library(data.table)
library(ggplot2)
library(patchwork)
library(survival)
library(cowplot)
library(timeROC)
library(grid)

cat("\n============================================================\n")
cat("RUNNING: FIGURE 3 FINAL V12 NO G; H AXIS MATCHES A/B\n")
cat("ROC objects: roc_expression and roc_cpg\n")
cat("============================================================\n\n")


# -----------------------------------------------------------------------------
# 1. Publication settings
# -----------------------------------------------------------------------------

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 245
FIGURE_DPI <- 600

FONT_FAMILY <- "sans"
FONT_PANEL_TAG <- 7.0
FONT_TITLE <- 6.6
FONT_AXIS_TITLE <- 5.8
FONT_AXIS_TEXT <- 5.3
FONT_TABLE <- 5.0
FONT_STAT <- 5.2

INK <- "#202124"
MUTED <- "#6B6F76"

# Panel A: red / orange / green / blue.
COMBINED_COLORS <- c(
    "Expression Low / Methylation High" = "#D55E00",
    "Expression Low / Methylation Low" = "#E69F00",
    "Expression High / Methylation High" = "#009E73",
    "Expression High / Methylation Low" = "#0072B2"
)

# Panels D and E use the same Low/High palette.
BINARY_COLORS <- c(
    "Low" = "#0072B2",
    "High" = "#D55E00"
)

ROC_COLORS <- c(
    "12 months" = "#3E8E8E",
    "24 months" = "#D07A43",
    "36 months" = "#574A7D"
)

COX_UNI_COLOR <- "#2E7D32"
COX_MULTI_COLOR <- "#17365D"


# -----------------------------------------------------------------------------
# 2. Input file
# -----------------------------------------------------------------------------

cat("\nSelect the completed TCGA-LGG DNMT3A master CSV file.\n")
input_file <- file.choose()

df <- data.table::fread(
    input_file,
    check.names = FALSE
)

cat("\nInput file loaded successfully.\n")
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n")


# -----------------------------------------------------------------------------
# 3. Output directory
# -----------------------------------------------------------------------------

choose_output_dir <- function() {

    selected <- NULL

    if (.Platform$OS.type == "windows") {
        selected <- tryCatch(
            utils::choose.dir(
                default = Sys.getenv("USERPROFILE"),
                caption = "Select a SHORT output folder for DNMT3A Figure 3"
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
        desktop <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
        if (!dir.exists(desktop)) {
            desktop <- getwd()
        }
        selected <- file.path(desktop, "DNMT3A_Fig3")
    }

    if (!dir.exists(selected)) {
        dir.create(selected, recursive = TRUE, showWarnings = FALSE)
    }

    normalizePath(selected, winslash = "/", mustWork = TRUE)
}

OUTPUT_DIR <- choose_output_dir()

cat("\nAll Figure 3 outputs will be saved to:\n")
cat(OUTPUT_DIR, "\n")


# -----------------------------------------------------------------------------
# 4. Required columns
# -----------------------------------------------------------------------------

required_columns <- c(
    "Patient_ID",
    "OS_months",
    "OS_status_raw",
    "PFS_months",
    "PFS_status_raw",
    "DSS_months",
    "DSS_status_raw",
    "DNMT3A_expression",
    "DNMT3A_expression_Median_Group",
    "cg21629895",
    "cg21629895_Median_Group",
    "DNMT3A_Promoter_Mean_Methylation_11Probe",
    "DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group",
    "Age..years.at.diagnosis.",
    "Gender",
    "EGFR_Amplification",
    "WHO2021_CDKN2AB_HomDel",
    "MGMT.promoter.status"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "Required columns are missing:\n",
            paste(missing_columns, collapse = ", ")
        )
    )
}


# -----------------------------------------------------------------------------
# 5. Utility functions
# -----------------------------------------------------------------------------

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

parse_event_status <- function(x) {

    x <- toupper(trimws(as.character(x)))

    data.table::fcase(
        is.na(x) | x == "", NA_integer_,
        grepl("^1:", x), 1L,
        grepl("^0:", x), 0L,
        x == "1", 1L,
        x == "0", 0L,
        default = NA_integer_
    )
}

safe_zscore <- function(x) {

    x <- safe_numeric(x)
    x_sd <- stats::sd(x, na.rm = TRUE)

    if (is.na(x_sd) || x_sd == 0) {
        return(rep(NA_real_, length(x)))
    }

    (x - mean(x, na.rm = TRUE)) / x_sd
}

pt_to_mm <- function(pt) {
    pt * 0.3527778
}

format_p <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
}


# -----------------------------------------------------------------------------
# 6. Prepare survival endpoints
# -----------------------------------------------------------------------------

df[, OS_time := safe_numeric(OS_months)]
df[, OS_event := parse_event_status(OS_status_raw)]
df[, PFS_time := safe_numeric(PFS_months)]
df[, PFS_event := parse_event_status(PFS_status_raw)]
df[, DSS_time := safe_numeric(DSS_months)]
df[, DSS_event := parse_event_status(DSS_status_raw)]

for (endpoint in c("OS", "PFS", "DSS")) {

    time_col <- paste0(endpoint, "_time")
    event_col <- paste0(endpoint, "_event")

    bad <- !is.na(df[[time_col]]) & df[[time_col]] < 0

    if (any(bad)) {
        df[[time_col]][bad] <- NA_real_
        df[[event_col]][bad] <- NA_integer_
    }
}


# -----------------------------------------------------------------------------
# 7. Prepare biomarker groups and Cox factors
# -----------------------------------------------------------------------------

df[, DNMT3A_expression := safe_numeric(DNMT3A_expression)]
df[, cg21629895 := safe_numeric(cg21629895)]
df[, DNMT3A_Promoter_Mean_Methylation_11Probe :=
       safe_numeric(DNMT3A_Promoter_Mean_Methylation_11Probe)]

# Median value used only for QC/reporting below.
# The categorical promoter High/Low groups continue to use the precomputed
# master-file column and are NOT recalculated here.
promoter_median <- stats::median(
    df$DNMT3A_Promoter_Mean_Methylation_11Probe,
    na.rm = TRUE
)

df[, Age_numeric := safe_numeric(`Age..years.at.diagnosis.`)]

# DNMT3A High vs Low; Low is the reference.
df[, DNMT3A_Group := factor(
    as.character(DNMT3A_expression_Median_Group),
    levels = c("Low", "High")
)]

# cg21629895 Cox coding: High is the reference so the HR represents Low vs High.
df[, cg21629895_Group_Cox := factor(
    as.character(cg21629895_Median_Group),
    levels = c("High", "Low")
)]

# cg21629895 KM display order uses Low then High for visual consistency.
df[, cg21629895_Group_KM := factor(
    as.character(cg21629895_Median_Group),
    levels = c("Low", "High")
)]

# Use the precomputed 11-CpG promoter median group from the master file.
# High is the reference so the Cox HR represents Low vs High.
df[, Promoter11_Group_Cox := factor(
    as.character(
        DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group
    ),
    levels = c("High", "Low")
)]

# Age >45 vs <=45.
df[, Age45_Group := factor(
    data.table::fcase(
        is.na(Age_numeric), NA_character_,
        Age_numeric > 45, ">45",
        Age_numeric <= 45, "<=45",
        default = NA_character_
    ),
    levels = c("<=45", ">45")
)]

# Gender.
df[, Sex := factor(
    data.table::fcase(
        tolower(trimws(as.character(Gender))) %in% c("female", "f"), "Female",
        tolower(trimws(as.character(Gender))) %in% c("male", "m"), "Male",
        default = NA_character_
    ),
    levels = c("Female", "Male")
)]

# CDKN2A/B homozygous deletion.
df[, CDKN2AB_Deletion := factor(
    data.table::fcase(
        safe_numeric(WHO2021_CDKN2AB_HomDel) == 0, "No",
        safe_numeric(WHO2021_CDKN2AB_HomDel) == 1, "Yes",
        default = NA_character_
    ),
    levels = c("No", "Yes")
)]

# EGFR amplification.
df[, EGFR_Amp := factor(
    data.table::fcase(
        safe_numeric(EGFR_Amplification) == 0, "No",
        safe_numeric(EGFR_Amplification) == 1, "Yes",
        default = NA_character_
    ),
    levels = c("No", "Yes")
)]

# MGMT promoter status.
df[, MGMT_Status := factor(
    trimws(as.character(MGMT.promoter.status)),
    levels = c("Methylated", "Unmethylated")
)]

# Continuous markers for univariate Cox and ROC.
df[, DNMT3A_expression_z := safe_zscore(DNMT3A_expression)]
df[, cg21629895_z := safe_zscore(cg21629895)]

# Risk-oriented continuous methylation variables for Cox regression.
# A 1-SD increase in these variables corresponds to a 1-SD DECREASE in the
# original methylation beta value. This makes the continuous Cox direction
# consistent with the categorical "Low vs High" comparisons.
df[, cg21629895_low_z := -safe_zscore(cg21629895)]
df[, Promoter11_low_z := -safe_zscore(
    DNMT3A_Promoter_Mean_Methylation_11Probe
)]

# Combined four-group biomarker variable.
df[, Combined_Biomarker_Group := data.table::fcase(
    DNMT3A_Group == "Low" & cg21629895_Group_KM == "High",
    "Expression Low / Methylation High",
    DNMT3A_Group == "Low" & cg21629895_Group_KM == "Low",
    "Expression Low / Methylation Low",
    DNMT3A_Group == "High" & cg21629895_Group_KM == "High",
    "Expression High / Methylation High",
    DNMT3A_Group == "High" & cg21629895_Group_KM == "Low",
    "Expression High / Methylation Low",
    default = NA_character_
)]

combined_levels <- c(
    "Expression Low / Methylation High",
    "Expression Low / Methylation Low",
    "Expression High / Methylation High",
    "Expression High / Methylation Low"
)

df[, Combined_Biomarker_Group := factor(
    Combined_Biomarker_Group,
    levels = combined_levels
)]


# -----------------------------------------------------------------------------
# 7B. Master-file median-group QC
# -----------------------------------------------------------------------------

master_group_qc <- data.table::data.table(
    Group = c(
        "DNMT3A expression median group",
        "cg21629895 median group",
        "11-CpG promoter mean median group"
    ),
    Source_Column = c(
        "DNMT3A_expression_Median_Group",
        "cg21629895_Median_Group",
        "DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group"
    ),
    Low_N = c(
        sum(df$DNMT3A_Group == "Low", na.rm = TRUE),
        sum(df$cg21629895_Group_Cox == "Low", na.rm = TRUE),
        sum(df$Promoter11_Group_Cox == "Low", na.rm = TRUE)
    ),
    High_N = c(
        sum(df$DNMT3A_Group == "High", na.rm = TRUE),
        sum(df$cg21629895_Group_Cox == "High", na.rm = TRUE),
        sum(df$Promoter11_Group_Cox == "High", na.rm = TRUE)
    ),
    Missing_N = c(
        sum(is.na(df$DNMT3A_Group)),
        sum(is.na(df$cg21629895_Group_Cox)),
        sum(is.na(df$Promoter11_Group_Cox))
    )
)

cat("\nMaster-file median-group QC:\n")
print(master_group_qc)


# -----------------------------------------------------------------------------
# 8. Survival QC
# -----------------------------------------------------------------------------

survival_qc <- data.table(
    Endpoint = c("OS", "PFS", "DSS"),
    N_complete = c(
        sum(complete.cases(df$OS_time, df$OS_event)),
        sum(complete.cases(df$PFS_time, df$PFS_event)),
        sum(complete.cases(df$DSS_time, df$DSS_event))
    ),
    N_events = c(
        sum(df$OS_event == 1, na.rm = TRUE),
        sum(df$PFS_event == 1, na.rm = TRUE),
        sum(df$DSS_event == 1, na.rm = TRUE)
    )
)

cat("\nSurvival QC:\n")
print(survival_qc)


# -----------------------------------------------------------------------------
# 9. Publication-style Kaplan-Meier plot with risk table
# -----------------------------------------------------------------------------

make_km_with_risk_table <- function(
        data,
        time_col,
        event_col,
        group_col,
        group_levels,
        colors,
        legend_labels,
        panel_letter,
        legend_nrow = 1,
        panel_title = NULL,
        show_confidence = TRUE,
        risk_table_labels = legend_labels
) {

    # Build a plain analysis data.frame using explicit standard names.
    # The Kaplan-Meier plot and risk table are constructed manually from
    # survival::survfit() output. This intentionally avoids survminer::ggsurvplot()
    # because that function can re-evaluate local object names and trigger
    # scoping errors such as "object 'km_data' not found" inside user functions.

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

    # -------------------------------------------------------------------------
    # Convert survfit output to plotting data.
    # -------------------------------------------------------------------------

    fit_summary <- summary(
        fit,
        censored = TRUE
    )

    strata_text <- as.character(
        fit_summary$strata
    )

    strata_text <- sub(
        "^KM_Group=",
        "",
        strata_text
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

    # Add time zero for each group so the curves and ribbons begin at survival 1.
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

    # -------------------------------------------------------------------------
    # Group counts and legend labels.
    # -------------------------------------------------------------------------

    group_counts <- table(
        factor(
            km_data$KM_Group,
            levels = group_levels
        )
    )

    final_legend_labels <- setNames(
        paste0(
            legend_labels,
            " (n = ",
            as.integer(
                group_counts
            ),
            ")"
        ),
        group_levels
    )

    # -------------------------------------------------------------------------
    # Risk table.
    # -------------------------------------------------------------------------

    x_max_raw <- max(
        km_data$Time,
        na.rm = TRUE
    )

    # Use 48-month intervals, while always including time zero.
    risk_times <- seq(
        0,
        ceiling(
            x_max_raw / 48
        ) * 48,
        by = 48
    )

    # Add a small right-side margin beyond the final labeled time point.
    # This prevents the 240-month tick and number-at-risk values from being
    # clipped when the plot is placed inside a multi-panel figure.
    x_display_max <- max(risk_times) + 0.035 * max(risk_times)

    risk_summary <- summary(
        fit,
        times = risk_times,
        extend = TRUE
    )

    risk_group <- sub(
        "^KM_Group=",
        "",
        as.character(
            risk_summary$strata
        )
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

    # -------------------------------------------------------------------------
    # Main survival plot.
    # -------------------------------------------------------------------------

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
            values = colors[
                group_levels
            ],
            breaks = group_levels,
            labels = unname(
                final_legend_labels[
                    group_levels
                ]
            ),
            name = NULL,
            drop = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values = colors[
                group_levels
            ],
            breaks = group_levels,
            guide = "none",
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(
                0,
                x_display_max
            ),
            breaks = risk_times,
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::scale_y_continuous(
            limits = c(
                0,
                1.02
            ),
            breaks = c(
                0,
                0.25,
                0.50,
                0.75,
                1
            ),
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::annotate(
            "label",
            x = 0.03 * x_display_max,
            y = 0.10,
            label = p_text,
            hjust = 0,
            size = pt_to_mm(
                5.0
            ),
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
                nrow = legend_nrow,
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
                margin = ggplot2::margin(
                    b = 2
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
            legend.justification = "left",
            legend.direction = "horizontal",
            legend.text = ggplot2::element_text(
                size = FONT_TABLE,
                color = INK
            ),
            legend.key.width = grid::unit(
                4.4,
                "mm"
            ),
            legend.key.height = grid::unit(
                2.0,
                "mm"
            ),
            plot.margin = ggplot2::margin(
                t = 1.5,
                r = 5,
                b = 0.5,
                l = 5
            )
        )

    # -------------------------------------------------------------------------
    # Manual number-at-risk table.
    # -------------------------------------------------------------------------

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
            size = pt_to_mm(
                FONT_TABLE
            ),
            show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
            values = colors[
                group_levels
            ],
            guide = "none",
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(
                -0.12 * max(risk_times),
                x_display_max
            ),
            breaks = risk_times,
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::scale_y_discrete(
            limits = rev(
                group_levels
            ),
            labels = rev(
                risk_table_labels
            )
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
                margin = ggplot2::margin(r = 2.0)
            ),
            axis.ticks = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(
                linewidth = 0.28,
                color = "#A7ABB1"
            ),
            plot.margin = ggplot2::margin(
                t = 0,
                r = 5,
                b = 1,
                l = 7
            )
        )

    combined_plot <- cowplot::plot_grid(
        survival_plot,
        risk_plot,
        ncol = 1,
        rel_heights = c(
            0.78,
            0.22
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
        n = nrow(
            km_data
        ),
        events = sum(
            km_data$Event == 1,
            na.rm = TRUE
        ),
        fit = fit,
        risk_table = risk_df
    )
}


# -----------------------------------------------------------------------------
# 10. Main OS Kaplan-Meier panels
# -----------------------------------------------------------------------------

km_A <- make_km_with_risk_table(
    data = df,
    time_col = "OS_time",
    event_col = "OS_event",
    group_col = "Combined_Biomarker_Group",
    group_levels = combined_levels,
    colors = COMBINED_COLORS,
    legend_labels = c(
        "RNA low / Meth high",
        "RNA low / Meth low",
        "RNA high / Meth high",
        "RNA high / Meth low"
    ),
    panel_letter = "H",
    legend_nrow = 2,
    panel_title = "Combined DNMT3A expression and cg21629895 methylation",
    show_confidence = FALSE,
    # Compact row labels prevent the four-group risk table from consuming
    # extra left-side width. The full group names remain in the legend.
    risk_table_labels = c(
        "L/H",
        "L/L",
        "H/H",
        "H/L"
    )
)

km_D <- make_km_with_risk_table(
    data = df,
    time_col = "OS_time",
    event_col = "OS_event",
    group_col = "DNMT3A_Group",
    group_levels = c("Low", "High"),
    colors = BINARY_COLORS,
    legend_labels = c("Low", "High"),
    panel_letter = "A",
    # Use the same two-row legend allocation as Panel H so the internal
    # plotting region has the same height across A, B, and H.
    legend_nrow = 2,
    panel_title = "DNMT3A expression",
    show_confidence = TRUE
)

km_E <- make_km_with_risk_table(
    data = df,
    time_col = "OS_time",
    event_col = "OS_event",
    group_col = "cg21629895_Group_KM",
    group_levels = c("Low", "High"),
    colors = BINARY_COLORS,
    legend_labels = c("Low", "High"),
    panel_letter = "B",
    # Use the same two-row legend allocation as Panel H so the internal
    # plotting region has the same height across A, B, and H.
    legend_nrow = 2,
    panel_title = "cg21629895 methylation",
    show_confidence = TRUE
)



# -----------------------------------------------------------------------------
# 11. Cox helpers
# -----------------------------------------------------------------------------

safe_coxph <- function(formula, data) {
    tryCatch(
        survival::coxph(
            formula,
            data = data,
            ties = "efron",
            x = TRUE,
            y = TRUE,
            model = TRUE
        ),
        error = function(e) NULL
    )
}

extract_cox <- function(fit, source_variable, display_label) {

    if (is.null(fit)) {
        return(data.table())
    }

    s <- summary(fit)
    coef_df <- as.data.frame(s$coefficients)
    ci_df <- as.data.frame(s$conf.int)

    data.table(
        Source_Variable = source_variable,
        Term = rownames(coef_df),
        Display_Label = display_label,
        HR = ci_df[["exp(coef)"]],
        CI_Lower = ci_df[["lower .95"]],
        CI_Upper = ci_df[["upper .95"]],
        P_Value = coef_df[["Pr(>|z|)"]],
        N = fit$n,
        Events = fit$nevent
    )
}


# -----------------------------------------------------------------------------
# 12. Univariate OS Cox analysis
# -----------------------------------------------------------------------------

univariate_specs <- list(
    list(
        formula = survival::Surv(OS_time, OS_event) ~ Age45_Group,
        variable = "Age45_Group",
        label = "Age (>45)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ DNMT3A_Group,
        variable = "DNMT3A_Group",
        label = "DNMT3A mRNA (High)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ cg21629895_Group_Cox,
        variable = "cg21629895_Group_Cox",
        label = "cg21629895 (Low)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ Promoter11_Group_Cox,
        variable = "Promoter11_Group_Cox",
        label = "11-CpG promoter (Low)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ Sex,
        variable = "Sex",
        label = "Gender (Male)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ CDKN2AB_Deletion,
        variable = "CDKN2AB_Deletion",
        label = "CDKN2A/B deletion"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ MGMT_Status,
        variable = "MGMT_Status",
        label = "MGMT (Unmethylated)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ EGFR_Amp,
        variable = "EGFR_Amp",
        label = "EGFR amplification"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ DNMT3A_expression_z,
        variable = "DNMT3A_expression_z",
        label = "DNMT3A expression (per SD)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ cg21629895_low_z,
        variable = "cg21629895_low_z",
        label = "cg21629895 (per SD lower)"
    ),
    list(
        formula = survival::Surv(OS_time, OS_event) ~ Promoter11_low_z,
        variable = "Promoter11_low_z",
        label = "11-CpG mean (per SD lower)"
    )
)

univariate_results <- data.table::rbindlist(
    lapply(
        univariate_specs,
        function(spec) {
            fit <- safe_coxph(spec$formula, df)
            extract_cox(fit, spec$variable, spec$label)
        }
    ),
    fill = TRUE
)


# -----------------------------------------------------------------------------
# 13. Multivariable OS Cox analysis
# -----------------------------------------------------------------------------

multivariable_formula <- survival::Surv(OS_time, OS_event) ~
    Age45_Group +
    DNMT3A_Group +
    cg21629895_Group_Cox +
    Sex +
    CDKN2AB_Deletion +
    MGMT_Status

multivariable_fit <- safe_coxph(multivariable_formula, df)

if (is.null(multivariable_fit)) {
    stop("The requested multivariable Cox model could not be fitted.")
}

multi_summary <- summary(multivariable_fit)
multi_coef <- as.data.frame(multi_summary$coefficients)
multi_ci <- as.data.frame(multi_summary$conf.int)

multi_label_map <- c(
    "Age45_Group>45" = "Age (>45)",
    "DNMT3A_GroupHigh" = "DNMT3A mRNA (High)",
    "cg21629895_Group_CoxLow" = "cg21629895 (Low)",
    "SexMale" = "Gender (Male)",
    "CDKN2AB_DeletionYes" = "CDKN2A/B deletion",
    "MGMT_StatusUnmethylated" = "MGMT (Unmethylated)"
)

multivariable_results <- data.table(
    Term = rownames(multi_coef),
    HR = multi_ci[["exp(coef)"]],
    CI_Lower = multi_ci[["lower .95"]],
    CI_Upper = multi_ci[["upper .95"]],
    P_Value = multi_coef[["Pr(>|z|)"]],
    N = multivariable_fit$n,
    Events = multivariable_fit$nevent
)

multivariable_results[, Display_Label := unname(multi_label_map[Term])]
multivariable_results[is.na(Display_Label), Display_Label := Term]


# -----------------------------------------------------------------------------
# 14. Cox forest plot in requested table + forest format
# -----------------------------------------------------------------------------

make_publication_cox_forest <- function(
        result_df,
        panel_letter,
        title_text,
        forest_color,
        desired_labels
) {

    plot_df <- data.table::copy(result_df)
    plot_df <- plot_df[
        !is.na(HR) &
            !is.na(CI_Lower) &
            !is.na(CI_Upper)
    ]

    plot_df[, Display_Order := match(Display_Label, desired_labels)]
    data.table::setorder(plot_df, Display_Order)

    # Use slightly expanded vertical spacing between Cox rows.
    row_step <- 1.12
    plot_df[, y := rev(seq_len(.N)) * row_step]

    plot_df[, HR_Text := sprintf("%.2f", HR)]
    plot_df[, CI_Text := sprintf("%.2f-%.2f", CI_Lower, CI_Upper)]
    plot_df[, P_Text := vapply(P_Value, format_p, character(1))]

    n_terms <- nrow(plot_df)
    y_top <- max(plot_df$y) + 0.90
    y_bottom <- 0.40

    # -------------------------------------------------------------------------
    # Left table
    # -------------------------------------------------------------------------
    table_panel <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(y = y)
    ) +
        ggplot2::geom_hline(
            yintercept = plot_df$y - row_step / 2,
            linewidth = 0.10,
            color = "#ECEDEF"
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                x = 0.02,
                label = Display_Label
            ),
            hjust = 0,
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                x = 2.08,
                label = HR_Text
            ),
            hjust = 0.5,
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                x = 2.86,
                label = CI_Text
            ),
            hjust = 0.5,
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                x = 3.92,
                label = P_Text
            ),
            hjust = 0.5,
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::annotate(
            "text",
            x = 0.02,
            y = y_top,
            label = "Variables",
            hjust = 0,
            fontface = "bold",
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::annotate(
            "text",
            x = 2.08,
            y = y_top,
            label = "HR",
            hjust = 0.5,
            fontface = "bold",
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::annotate(
            "text",
            x = 2.86,
            y = y_top,
            label = "95% CI",
            hjust = 0.5,
            fontface = "bold",
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::annotate(
            "text",
            x = 3.92,
            y = y_top,
            label = "p-value",
            hjust = 0.5,
            fontface = "bold",
            size = pt_to_mm(5.00),
            color = INK
        ) +
        ggplot2::scale_y_continuous(
            limits = c(
                y_bottom,
                y_top + 0.35
            ),
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::xlim(
            0,
            4.30
        ) +
        ggplot2::coord_cartesian(
            clip = "off"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                t = 2.5,
                r = 1.2,
                b = 1,
                l = 1.5
            )
        )

    # -------------------------------------------------------------------------
    # Forest core
    # -------------------------------------------------------------------------
    x_lower <- max(
        min(
            plot_df$CI_Lower,
            na.rm = TRUE
        ) * 0.85,
        0.25
    )

    x_upper <- min(
        max(
            plot_df$CI_Upper,
            na.rm = TRUE
        ) * 1.08,
        16
    )

    if (x_lower >= 1) {
        x_lower <- 0.5
    }

    if (x_upper <= 1) {
        x_upper <- 2
    }

    forest_panel <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = HR,
            y = y
        )
    ) +
        ggplot2::geom_hline(
            yintercept = plot_df$y - row_step / 2,
            linewidth = 0.10,
            color = "#ECEDEF"
        ) +
        ggplot2::geom_vline(
            xintercept = 1,
            linewidth = 0.30,
            linetype = "22",
            color = "#A8ADB3"
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = pmax(
                    CI_Lower,
                    x_lower
                ),
                xend = pmin(
                    CI_Upper,
                    x_upper
                ),
                yend = y
            ),
            linewidth = 0.48,
            color = forest_color,
            lineend = "round"
        ) +
        ggplot2::geom_point(
            shape = 15,
            size = 1.18,
            color = forest_color
        ) +
        ggplot2::scale_x_log10(
            limits = c(
                x_lower,
                x_upper
            ),
            breaks = c(
                0.3,
                1,
                3,
                10
            )
        ) +
        ggplot2::scale_y_continuous(
            limits = c(
                y_bottom,
                y_top + 0.35
            ),
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::labs(
            x = "Hazard ratio (95% CI)",
            y = NULL
        ) +
        ggplot2::theme_classic(
            base_size = 7,
            base_family = FONT_FAMILY
        ) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(
                size = 5.0,
                color = INK,
                margin = ggplot2::margin(
                    t = 1.5
                )
            ),
            axis.text.x = ggplot2::element_text(
                size = 5.0,
                color = INK
            ),
            axis.line.x = ggplot2::element_line(
                linewidth = 0.30,
                color = INK
            ),
            plot.margin = ggplot2::margin(
                t = 2.5,
                r = 4,
                b = 1,
                l = 1
            )
        )

    # Keep a generous table area while leaving enough room for the HR axis.
    combined <- table_panel + forest_panel +
        patchwork::plot_layout(
            widths = c(
                2.18,
                0.82
            )
        )

    cowplot::ggdraw(
        combined
    ) +
        cowplot::draw_label(
            title_text,
            x = 0.50,
            y = 0.992,
            hjust = 0.5,
            vjust = 1,
            fontface = "bold",
            size = 6.15,
            color = INK
        ) +
        cowplot::draw_label(
            panel_letter,
            x = 0.003,
            y = 0.997,
            hjust = 0,
            vjust = 1,
            fontface = "bold",
            size = FONT_PANEL_TAG,
            color = INK
        )
}


univariate_order <- c(
    "Age (>45)",
    "DNMT3A mRNA (High)",
    "DNMT3A expression (per SD)",
    "cg21629895 (Low)",
    "cg21629895 (per SD lower)",
    "11-CpG promoter (Low)",
    "11-CpG mean (per SD lower)",
    "Gender (Male)",
    "CDKN2A/B deletion",
    "MGMT (Unmethylated)",
    "EGFR amplification"
)

multivariable_order <- c(
    "Age (>45)",
    "DNMT3A mRNA (High)",
    "cg21629895 (Low)",
    "Gender (Male)",
    "CDKN2A/B deletion",
    "MGMT (Unmethylated)"
)

panel_E <- make_publication_cox_forest(
    result_df = univariate_results,
    panel_letter = "E",
    title_text = "Univariate Analysis",
    forest_color = COX_UNI_COLOR,
    desired_labels = univariate_order
)

panel_F <- make_publication_cox_forest(
    result_df = multivariable_results,
    panel_letter = "F",
    title_text = "Multivariate Analysis",
    forest_color = COX_MULTI_COLOR,
    desired_labels = multivariable_order
)


# -----------------------------------------------------------------------------
# 15. Time-dependent ROC helper
# -----------------------------------------------------------------------------

make_time_roc <- function(time, event, marker, panel_letter, subtitle_text) {

    roc_input <- data.table(
        Time = safe_numeric(time),
        Event = safe_numeric(event),
        Marker = safe_numeric(marker)
    )

    roc_input <- roc_input[
        complete.cases(Time, Event, Marker) &
        Time >= 0 &
        Event %in% c(0, 1)
    ]

    times_eval <- c(12, 24, 36)

    roc_fit <- timeROC::timeROC(
        T = roc_input$Time,
        delta = roc_input$Event,
        marker = roc_input$Marker,
        cause = 1,
        weighting = "marginal",
        times = times_eval,
        ROC = TRUE,
        iid = TRUE
    )

    plot_data <- data.table::rbindlist(
        lapply(
            seq_along(times_eval),
            function(i) {
                data.table(
                    FPR = roc_fit$FP[, i],
                    TPR = roc_fit$TP[, i],
                    Timepoint = paste0(times_eval[i], " months")
                )
            }
        )
    )

    plot_data[, Timepoint := factor(
        Timepoint,
        levels = names(ROC_COLORS)
    )]

    auc_annotation <- data.table(
        x = 0.96,
        y = c(0.31, 0.22, 0.13),
        Timepoint = factor(
            names(ROC_COLORS),
            levels = names(ROC_COLORS)
        ),
        Label = paste0(
            c("12 mo", "24 mo", "36 mo"),
            "  AUC ",
            sprintf("%.2f", roc_fit$AUC)
        )
    )

    p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = FPR,
            y = TPR,
            color = Timepoint,
            group = Timepoint
        )
    ) +
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            linetype = "22",
            linewidth = 0.28,
            color = "#A4A7AB"
        ) +
        ggplot2::geom_step(
            linewidth = 0.60,
            direction = "hv"
        ) +
        ggplot2::geom_text(
            data = auc_annotation,
            ggplot2::aes(
                x = x,
                y = y,
                label = Label,
                color = Timepoint
            ),
            inherit.aes = FALSE,
            hjust = 1,
            size = pt_to_mm(5.1),
            fontface = "bold"
        ) +
        ggplot2::scale_color_manual(
            values = ROC_COLORS,
            guide = "none"
        ) +
        ggplot2::scale_x_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.25),
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.25),
            expand = c(0, 0)
        ) +
        ggplot2::labs(
            subtitle = subtitle_text,
            x = "1 - Specificity",
            y = "Sensitivity"
        ) +
        ggplot2::theme_classic(base_size = 7) +
        ggplot2::theme(
            plot.subtitle = ggplot2::element_text(
                size = FONT_TITLE,
                face = "bold",
                hjust = 0,
                color = INK,
                margin = ggplot2::margin(b = 2)
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
                linewidth = 0.30,
                color = INK
            ),
            plot.margin = ggplot2::margin(2, 4, 2, 3)
        )

    wrapped <- cowplot::ggdraw(p) +
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
        plot = wrapped,
        fit = roc_fit,
        n = nrow(roc_input),
        auc = roc_fit$AUC
    )
}

roc_expression <- make_time_roc(
    time = df$OS_time,
    event = df$OS_event,
    marker = df$DNMT3A_expression_z,
    panel_letter = "C",
    subtitle_text = "DNMT3A expression"
)

roc_cpg <- make_time_roc(
    time = df$OS_time,
    event = df$OS_event,
    marker = -df$cg21629895_z,
    panel_letter = "D",
    subtitle_text = "cg21629895 methylation"
)


# -----------------------------------------------------------------------------
# 16. Assemble final seven-panel figure
#
# Layout:
#   A. DNMT3A expression OS Kaplan-Meier
#   B. cg21629895 methylation OS Kaplan-Meier
#   C. DNMT3A expression time-dependent ROC
#   D. cg21629895 methylation time-dependent ROC
#   E. Univariate OS Cox
#   F. Multivariable OS Cox
#   G. Removed
#   H. Combined four-group OS Kaplan-Meier
#
# The former integrated molecular states scatter panel has been removed.
# The combined four-group OS Kaplan-Meier panel remains labeled H to preserve
# the existing panel lettering used in the manuscript.
# -----------------------------------------------------------------------------

panel_A <- km_D$plot
panel_B <- km_E$plot
panel_C <- roc_expression$plot
panel_D <- roc_cpg$plot
# panel_E and panel_F are the Cox plots created above.
panel_H <- km_A$plot

row_1 <- panel_A + panel_B +
    patchwork::plot_layout(
        widths = c(
            1,
            1
        )
    )

row_2 <- panel_C + panel_D +
    patchwork::plot_layout(
        widths = c(
            1,
            1
        )
    )

row_3 <- panel_E + panel_F +
    patchwork::plot_layout(
        widths = c(
            1,
            1
        )
    )

# Panel H uses the same half-page-width column geometry as Panels A and B.
row_4 <- panel_H + patchwork::plot_spacer() +
    patchwork::plot_layout(
        widths = c(
            1,
            1
        )
    )

# Add a small, uniform vertical gutter between figure rows.
# The spacer rows are equal-height, so A/B-C/D, C/D-E/F, and E/F-H
# have visually consistent vertical separation. Panel contents are unchanged.
vertical_gap <- patchwork::plot_spacer()

final_figure_3 <- row_1 / vertical_gap /
                  row_2 / vertical_gap /
                  row_3 / vertical_gap /
                  row_4 +
    patchwork::plot_layout(
        heights = c(
            1.08,
            0.055,
            0.88,
            0.055,
            1.22,
            0.055,
            1.08
        )
    )


# -----------------------------------------------------------------------------
# 17. Save helpers
# -----------------------------------------------------------------------------

save_pdf <- function(plot_object, filename, width_mm, height_mm) {

    grDevices::pdf(
        file = filename,
        width = width_mm / 25.4,
        height = height_mm / 25.4,
        family = "sans",
        useDingbats = FALSE
    )

    print(plot_object)
    grDevices::dev.off()
}

save_png <- function(plot_object, filename, width_mm, height_mm) {

    ggplot2::ggsave(
        filename = filename,
        plot = plot_object,
        width = width_mm,
        height = height_mm,
        units = "mm",
        dpi = FIGURE_DPI,
        bg = "white"
    )
}


# -----------------------------------------------------------------------------
# 18. Save combined figure
# -----------------------------------------------------------------------------

combined_pdf <- file.path(
    OUTPUT_DIR,
    "Fig3_Combined_180mm_Prognostic_Significance_FINAL_V12_NO_G_H_AXIS_MATCH_AB.pdf"
)

combined_png <- file.path(
    OUTPUT_DIR,
    "Fig3_Combined_180mm_Prognostic_Significance_FINAL_V12_NO_G_H_AXIS_MATCH_AB_600dpi.png"
)

save_pdf(
    final_figure_3,
    combined_pdf,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)

save_png(
    final_figure_3,
    combined_png,
    FIGURE_WIDTH_MM,
    FIGURE_HEIGHT_MM
)


# -----------------------------------------------------------------------------
# 19. Save individual panels
# -----------------------------------------------------------------------------

save_pdf(
    panel_A,
    file.path(
        OUTPUT_DIR,
        "Fig3C_DNMT3A_TimeROC.pdf"
    ),
    90,
    65
)

save_pdf(
    panel_B,
    file.path(
        OUTPUT_DIR,
        "Fig3D_cg21629895_TimeROC.pdf"
    ),
    90,
    65
)

save_pdf(
    panel_C,
    file.path(
        OUTPUT_DIR,
        "Fig3A_DNMT3A_OS_KM.pdf"
    ),
    90,
    65
)

save_pdf(
    panel_D,
    file.path(
        OUTPUT_DIR,
        "Fig3B_cg21629895_OS_KM.pdf"
    ),
    90,
    65
)

save_pdf(
    panel_E,
    file.path(
        OUTPUT_DIR,
        "Fig3E_Univariate_OS_Cox.pdf"
    ),
    90,
    66
)

save_pdf(
    panel_F,
    file.path(
        OUTPUT_DIR,
        "Fig3F_Multivariable_OS_Cox.pdf"
    ),
    90,
    66
)

save_pdf(
    panel_H,
    file.path(
        OUTPUT_DIR,
        "Fig3H_Combined_OS_KM.pdf"
    ),
    90,
    65
)


# -----------------------------------------------------------------------------
# 20. Save Cox tables
# -----------------------------------------------------------------------------

data.table::fwrite(
    univariate_results,
    file.path(OUTPUT_DIR, "Fig3_Univariate_OS_Cox_Results.csv"),
    na = ""
)

data.table::fwrite(
    multivariable_results,
    file.path(OUTPUT_DIR, "Fig3_Multivariable_OS_Cox_Results.csv"),
    na = ""
)


# -----------------------------------------------------------------------------
# 21. Cox proportional-hazards assumption
# -----------------------------------------------------------------------------

multivariable_zph <- survival::cox.zph(multivariable_fit)

multivariable_zph_table <- data.table::as.data.table(
    as.data.frame(multivariable_zph$table),
    keep.rownames = "Term"
)

data.table::fwrite(
    multivariable_zph_table,
    file.path(OUTPUT_DIR, "Fig3_Multivariable_Cox_PH_Assumption.csv"),
    na = ""
)



# -----------------------------------------------------------------------------
# ROC object integrity check
# -----------------------------------------------------------------------------

if (!exists("roc_expression", inherits = FALSE)) {
    stop("Internal error: roc_expression was not created.")
}

if (!exists("roc_cpg", inherits = FALSE)) {
    stop("Internal error: roc_cpg was not created.")
}

if (length(roc_expression$auc) != 3 || length(roc_cpg$auc) != 3) {
    stop("Internal error: expected exactly three AUC values for each ROC analysis.")
}

# -----------------------------------------------------------------------------
# 22. Save ROC statistics
# -----------------------------------------------------------------------------

roc_results <- data.table(
    Biomarker = rep(
        c("DNMT3A expression", "cg21629895 methylation"),
        each = 3
    ),
    Time_months = rep(c(12, 24, 36), 2),
    AUC = c(roc_expression$auc, roc_cpg$auc),
    N = c(
        rep(roc_expression$n, 3),
        rep(roc_cpg$n, 3)
    )
)

data.table::fwrite(
    roc_results,
    file.path(OUTPUT_DIR, "Fig3_TimeDependent_ROC_Results.csv"),
    na = ""
)


# -----------------------------------------------------------------------------
# 23. Save KM statistics and promoter cutoff
# -----------------------------------------------------------------------------

km_statistics <- data.table(
    Panel = c("A", "D", "E"),
    Analysis = c(
        "Combined expression/methylation",
        "DNMT3A expression",
        "cg21629895 methylation"
    ),
    N = c(km_A$n, km_D$n, km_E$n),
    Events = c(km_A$events, km_D$events, km_E$events),
    LogRank_P = c(km_A$p_value, km_D$p_value, km_E$p_value)
)

data.table::fwrite(
    km_statistics,
    file.path(OUTPUT_DIR, "Fig3_OS_KM_Statistics.csv"),
    na = ""
)

promoter_cutoff_qc <- data.table(
    Variable = "DNMT3A 11-CpG promoter mean methylation",
    Median_Cutoff = promoter_median,
    Definition = "Low < median; High >= median"
)

data.table::fwrite(
    promoter_cutoff_qc,
    file.path(OUTPUT_DIR, "Fig3_11CpG_Promoter_Median_Cutoff.csv"),
    na = ""
)


# -----------------------------------------------------------------------------
# 24. Save survival QC and overall QC
# -----------------------------------------------------------------------------

data.table::fwrite(
    survival_qc,
    file.path(OUTPUT_DIR, "Fig3_Survival_Endpoint_QC.csv"),
    na = ""
)

qc <- data.table(
    Metric = c(
        "Master rows",
        "OS evaluable samples",
        "OS events",
        "Age cutoff",
        "Promoter methylation definition",
        "Univariate Cox variables",
        "Multivariable Cox variables",
        "Multivariable Cox N",
        "Multivariable Cox events",
        "ROC time points",
        "Figure width mm",
        "Figure height mm"
    ),
    Value = c(
        nrow(df),
        sum(complete.cases(df$OS_time, df$OS_event)),
        sum(df$OS_event == 1, na.rm = TRUE),
        ">45 vs <=45 years",
        "11-CpG mean; Low < median, High >= median",
        paste(univariate_order, collapse = " | "),
        paste(multivariable_order, collapse = " | "),
        multivariable_fit$n,
        multivariable_fit$nevent,
        "12, 24, 36 months",
        FIGURE_WIDTH_MM,
        FIGURE_HEIGHT_MM
    )
)

data.table::fwrite(
    qc,
    file.path(OUTPUT_DIR, "Fig3_QC.csv"),
    na = ""
)


data.table::fwrite(
    master_group_qc,
    file.path(
        OUTPUT_DIR,
        "Fig3_Master_Median_Group_QC.csv"
    ),
    na = ""
)


# -----------------------------------------------------------------------------
# 25. Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MAIN FIGURE 3 FINAL V12 COMPLETED\n")
cat("============================================================\n")

cat("\nMaster-file group definitions:\n")
cat("DNMT3A expression High/Low: DNMT3A_expression_Median_Group\n")
cat("cg21629895 High/Low: cg21629895_Median_Group\n")
cat("11-CpG promoter High/Low: DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group\n")

cat("\nFinal Cox orientation notes:\n")
cat("DNMT3A expression continuous HR: per 1-SD increase in expression\n")
cat("cg21629895 continuous HR: per 1-SD decrease in methylation\n")
cat("11-CpG promoter continuous HR: per 1-SD decrease in methylation\n")
cat("Panel H confidence shading: OFF\n")
cat("Panels A/B/H KM geometry: matched 90 x 65 mm, matched row height, matched two-row legend allocation; H risk labels compacted to preserve axis width\n")

cat("\nPanel structure:\n")
cat("A. DNMT3A expression OS KM with risk table\n")
cat("B. cg21629895 methylation OS KM with risk table\n")
cat("C. DNMT3A expression time-dependent ROC\n")
cat("D. cg21629895 methylation time-dependent ROC\n")
cat("E. Univariate OS Cox\n")
cat("F. Multivariable OS Cox\n")
cat("G. Removed\n")
cat("H. Combined OS KM with risk table\n")

cat("\nUnivariate Cox factors:\n")
cat(paste(univariate_order, collapse = "\n"), "\n")

cat("\nMultivariable Cox factors:\n")
cat(paste(multivariable_order, collapse = "\n"), "\n")

cat("\n11-CpG promoter median cutoff:\n")
cat(promoter_median, "\n")

cat("\nMultivariable model:\n")
cat("N =", multivariable_fit$n, "\n")
cat("Events =", multivariable_fit$nevent, "\n")

cat("\nCombined figure:\n")
cat(combined_pdf, "\n")
cat(combined_png, "\n")

cat("\n============================================================\n")
