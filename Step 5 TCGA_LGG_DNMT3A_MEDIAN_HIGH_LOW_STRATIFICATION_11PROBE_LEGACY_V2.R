# =============================================================================
# TCGA-LGG DNMT3A Median-Based High/Low Stratification
#
# This script:
# 1. Selects the current TCGA-LGG DNMT3A master CSV.
# 2. Identifies DNMT3A mRNA expression, all cg methylation probes, and the
#    11-probe promoter-associated CpG mean.
# 3. Calculates the median separately for each continuous variable using
#    available values only.
# 4. Adds a High/Low group column immediately after each corresponding variable.
# 5. Saves the updated master file and QC files.
#
# Stratification rule:
#   Low  = value < median
#   High = value >= median
#   NA   = missing original value
#
# NOTE:
# Values exactly equal to the median are assigned to the High group.
#
# Current analysis context:
#   - TCGA-LGG legacy HiSeqV2 DNMT3A expression
#   - TCGA-LGG legacy HumanMethylation450 beta values
#   - DNMT3A 11-probe promoter-associated methylation mean
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install/load package
# -----------------------------------------------------------------------------

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}

library(data.table)


# -----------------------------------------------------------------------------
# 1. Select input master file
# -----------------------------------------------------------------------------

cat("\nSelect the current TCGA-LGG DNMT3A master CSV file containing the 11-probe promoter mean.\n")

input_file <- file.choose()

master <- data.table::fread(
    input_file,
    check.names = FALSE
)

cat("\nInput file loaded successfully.\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")


# -----------------------------------------------------------------------------
# 2. Define variables to stratify
# -----------------------------------------------------------------------------

expression_col <- "DNMT3A_expression"

methylation_cols <- grep(
    "^cg[0-9]+$",
    names(master),
    value = TRUE
)

promoter_mean_col <- "DNMT3A_Promoter_Mean_Methylation_11Probe"

variables_to_stratify <- c(
    expression_col,
    methylation_cols
)

if (promoter_mean_col %in% names(master)) {
    variables_to_stratify <- c(
        variables_to_stratify,
        promoter_mean_col
    )
}


# -----------------------------------------------------------------------------
# 3. Confirm required variables
# -----------------------------------------------------------------------------

if (!(expression_col %in% names(master))) {
    stop(
        paste0(
            "Required expression column was not found: ",
            expression_col
        )
    )
}

if (length(methylation_cols) == 0) {
    stop(
        "No methylation probe columns matching ^cg[0-9]+$ were found."
    )
}

cat("\nVariables detected for median stratification:\n")
cat("DNMT3A expression:", expression_col, "\n")
cat("Methylation probes:", length(methylation_cols), "\n")

if (promoter_mean_col %in% variables_to_stratify) {
    cat("11-probe promoter mean: Yes\n")
} else {
    cat("11-probe promoter mean: No\n")
}

cat("Total variables to stratify:", length(variables_to_stratify), "\n")


# -----------------------------------------------------------------------------
# 4. Convert variables to numeric and calculate median-based groups
# -----------------------------------------------------------------------------

qc_list <- vector(
    "list",
    length(variables_to_stratify)
)

names(qc_list) <- variables_to_stratify

new_group_columns <- list()

for (variable in variables_to_stratify) {

    values <- suppressWarnings(
        as.numeric(master[[variable]])
    )

    master[[variable]] <- values

    n_available <- sum(
        !is.na(values)
    )

    n_missing <- sum(
        is.na(values)
    )

    if (n_available == 0) {

        median_value <- NA_real_

        group <- rep(
            NA_character_,
            length(values)
        )

        n_low <- 0
        n_high <- 0
        n_equal_median <- 0

    } else {

        median_value <- median(
            values,
            na.rm = TRUE
        )

        group <- ifelse(
            is.na(values),
            NA_character_,
            ifelse(
                values < median_value,
                "Low",
                "High"
            )
        )

        n_low <- sum(
            group == "Low",
            na.rm = TRUE
        )

        n_high <- sum(
            group == "High",
            na.rm = TRUE
        )

        n_equal_median <- sum(
            values == median_value,
            na.rm = TRUE
        )
    }

    group_col <- paste0(
        variable,
        "_Median_Group"
    )

    new_group_columns[[variable]] <- list(
        group_name = group_col,
        group_values = group
    )

    qc_list[[variable]] <- data.frame(
        Variable = variable,
        Variable_Type = ifelse(
            variable == expression_col,
            "mRNA expression",
            ifelse(
                variable == promoter_mean_col,
                "Mean methylation across 11 promoter-associated CpGs",
                "Individual methylation probe"
            )
        ),
        Median = median_value,
        N_Total = length(values),
        N_Available = n_available,
        N_Missing = n_missing,
        N_Low = n_low,
        N_High = n_high,
        N_Equal_to_Median = n_equal_median,
        Low_Definition = "value < median",
        High_Definition = "value >= median",
        stringsAsFactors = FALSE
    )
}


# -----------------------------------------------------------------------------
# 5. Insert each High/Low group immediately after its source variable
# -----------------------------------------------------------------------------

original_names <- names(master)

for (variable in rev(variables_to_stratify)) {

    group_col <- new_group_columns[[variable]]$group_name
    group_values <- new_group_columns[[variable]]$group_values

    source_position <- match(
        variable,
        names(master)
    )

    master[
        ,
        (group_col) := group_values
    ]

    current_names <- names(master)

    current_names <- current_names[
        current_names != group_col
    ]

    source_position <- match(
        variable,
        current_names
    )

    reordered_names <- append(
        current_names,
        group_col,
        after = source_position
    )

    data.table::setcolorder(
        master,
        reordered_names
    )
}


# -----------------------------------------------------------------------------
# 6. Combine QC table
# -----------------------------------------------------------------------------

qc_table <- data.table::rbindlist(
    qc_list,
    fill = TRUE
)


# -----------------------------------------------------------------------------
# 7. Additional QC for balance and missingness
# -----------------------------------------------------------------------------

qc_table[
    ,
    Low_Percent_of_Available := round(
        N_Low / N_Available * 100,
        2
    )
]

qc_table[
    ,
    High_Percent_of_Available := round(
        N_High / N_Available * 100,
        2
    )
]

qc_table[
    N_Available == 0,
    `:=`(
        Low_Percent_of_Available = NA_real_,
        High_Percent_of_Available = NA_real_
    )
]

qc_table[
    ,
    Group_Count_Check := ifelse(
        N_Low + N_High == N_Available,
        "PASS",
        "FAIL"
    )
]


# -----------------------------------------------------------------------------
# 8. Overall summary QC
# -----------------------------------------------------------------------------

summary_qc <- data.frame(
    Metric = c(
        "Total samples",
        "Total variables stratified",
        "mRNA expression variables stratified",
        "Individual methylation probes stratified",
        "Promoter mean variables stratified",
        "Variables with no available values",
        "Variables passing Low + High = available count check",
        "Variables failing Low + High = available count check"
    ),
    Value = c(
        nrow(master),
        nrow(qc_table),
        sum(
            qc_table$Variable_Type == "mRNA expression"
        ),
        sum(
            qc_table$Variable_Type == "Individual methylation probe"
        ),
        sum(
            qc_table$Variable_Type ==
                "Mean methylation across 11 promoter-associated CpGs"
        ),
        sum(
            qc_table$N_Available == 0
        ),
        sum(
            qc_table$Group_Count_Check == "PASS"
        ),
        sum(
            qc_table$Group_Count_Check == "FAIL"
        )
    ),
    stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 9. Select a short output directory and save outputs
#
# Deep OneDrive paths can exceed the Windows path length handled by fwrite().
# Therefore, the user is asked to choose a short output folder.
# If folder selection is cancelled, a short folder on the Desktop is used.
# -----------------------------------------------------------------------------

choose_output_dir <- function() {

    selected_dir <- NULL

    if (.Platform$OS.type == "windows") {

        selected_dir <- tryCatch(
            utils::choose.dir(
                default = path.expand("~/Desktop"),
                caption = "Select a SHORT output folder for DNMT3A results"
            ),
            error = function(e) NULL
        )
    }

    if (
        is.null(selected_dir) ||
        is.na(selected_dir) ||
        !nzchar(selected_dir)
    ) {

        desktop_dir <- file.path(
            path.expand("~"),
            "Desktop"
        )

        if (!dir.exists(desktop_dir)) {
            desktop_dir <- tempdir()
        }

        selected_dir <- file.path(
            desktop_dir,
            "DNMT3A_11Probe_Output"
        )
    }

    if (!dir.exists(selected_dir)) {
        dir.create(
            selected_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    selected_dir <- normalizePath(
        selected_dir,
        winslash = "/",
        mustWork = TRUE
    )

    return(selected_dir)
}


output_dir <- choose_output_dir()

cat("\nOutput directory selected:\n")
cat(output_dir, "\n")

master_output <- file.path(
    output_dir,
    "DNMT3A_MASTER_11PROBE_HIGH_LOW.csv"
)

qc_output <- file.path(
    output_dir,
    "DNMT3A_11PROBE_HIGH_LOW_QC.csv"
)

summary_output <- file.path(
    output_dir,
    "DNMT3A_11PROBE_HIGH_LOW_SUMMARY.csv"
)


safe_fwrite <- function(
    x,
    output_file
) {

    tryCatch(
        {
            data.table::fwrite(
                x,
                output_file,
                na = ""
            )
        },
        error = function(e) {

            stop(
                paste0(
                    "\nFailed to save file:\n",
                    output_file,
                    "\n\nReason:\n",
                    conditionMessage(e),
                    "\n\nPlease rerun the script and select a shorter local ",
                    "output folder, such as C:/DNMT3A_Output or Desktop."
                ),
                call. = FALSE
            )
        }
    )
}


safe_fwrite(
    master,
    master_output
)

safe_fwrite(
    qc_table,
    qc_output
)

safe_fwrite(
    summary_qc,
    summary_output
)


# -----------------------------------------------------------------------------
# 10. Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MEDIAN-BASED HIGH/LOW STRATIFICATION COMPLETED\n")
cat("============================================================\n")

cat("\nRule used:\n")
cat("Low  = value < median\n")
cat("High = value >= median\n")
cat("Missing original value = NA group\n")

cat("\nSummary QC:\n")
print(
    summary_qc
)

cat("\nDNMT3A expression QC:\n")
print(
    qc_table[
        Variable == expression_col
    ]
)

if (promoter_mean_col %in% qc_table$Variable) {
    cat("\n11-probe promoter mean QC:\n")
    print(
        qc_table[
            Variable == promoter_mean_col
        ]
    )
}

cat("\nMethylation probes with missing values:\n")
print(
    qc_table[
        Variable_Type == "Individual methylation probe" &
        N_Missing > 0,
        .(
            Variable,
            Median,
            N_Available,
            N_Missing,
            N_Low,
            N_High,
            N_Equal_to_Median
        )
    ]
)

cat("\nOutput master:\n")
cat(master_output, "\n")

cat("\nVariable-level QC:\n")
cat(qc_output, "\n")

cat("\nSummary QC:\n")
cat(summary_output, "\n")
