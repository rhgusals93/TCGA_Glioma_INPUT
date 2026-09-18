# =============================================================================
# Add DNMT3A 11-Probe Promoter Mean Methylation + QC
#
# Current analysis context:
#   - TCGA-LGG legacy HumanMethylation450 beta values
#   - hg19/GRCh37-compatible methylation annotation
#   - Current integrated/classified TCGA-LGG master CSV
#
# Promoter definition:
# The DNMT3A promoter-associated set contains 11 HumanMethylation450 probes
# annotated to TSS200 and/or TSS1500 regions.
#
# All 11 probes are included in the promoter methylation score:
#
#   cg15998962
#   cg03463641
#   cg02208653
#   cg21708767
#   cg13344237
#   cg06112956
#   cg13076778
#   cg21629895
#   cg26470599
#   cg22731525
#   cg07150430
#
# The promoter methylation score is calculated as the row-wise mean beta value
# across all available values among these 11 probes.
#
# Output:
#   1. Updated master CSV with the 11-probe promoter mean
#   2. Sample-level QC CSV
#   3. Probe-level QC CSV
#   4. Summary QC CSV
#   5. Beta-value range QC CSV
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install/load packages
# -----------------------------------------------------------------------------

required_packages <- c(
    "data.table",
    "dplyr"
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
    install.packages(missing_packages)
}

library(data.table)
library(dplyr)


# -----------------------------------------------------------------------------
# 1. Select the current TCGA-LGG master CSV
# -----------------------------------------------------------------------------

cat("\nSelect the current TCGA-LGG DNMT3A master CSV file.\n")

input_file <- file.choose()

master <- data.table::fread(
    input_file,
    check.names = FALSE
)

cat("\nInput file loaded successfully.\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")


# -----------------------------------------------------------------------------
# 2. Define the 11 DNMT3A promoter-associated CpG probes
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# 3. Confirm that all 11 promoter probes are present
# -----------------------------------------------------------------------------

missing_probes <- setdiff(
    promoter_probes_11,
    names(master)
)

if (length(missing_probes) > 0) {
    stop(
        paste(
            "The following DNMT3A promoter-associated probes are missing:",
            paste(
                missing_probes,
                collapse = ", "
            )
        )
    )
}

cat("\nAll 11 DNMT3A promoter-associated CpG probes are present.\n")


# -----------------------------------------------------------------------------
# 4. Convert promoter probe columns to numeric
# -----------------------------------------------------------------------------

for (probe in promoter_probes_11) {

    master[[probe]] <- suppressWarnings(
        as.numeric(master[[probe]])
    )
}


# -----------------------------------------------------------------------------
# 5. Probe-level QC
# -----------------------------------------------------------------------------

probe_qc <- lapply(
    promoter_probes_11,
    function(probe) {

        values <- master[[probe]]
        non_missing_values <- values[!is.na(values)]

        data.frame(
            Probe = probe,
            Included_in_11Probe_Mean = "Yes",
            N_Total = length(values),
            N_Available = sum(!is.na(values)),
            N_Missing = sum(is.na(values)),
            Completeness_Percent = round(
                mean(!is.na(values)) * 100,
                2
            ),
            Mean_Beta = if (
                length(non_missing_values) > 0
            ) {
                mean(non_missing_values)
            } else {
                NA_real_
            },
            Median_Beta = if (
                length(non_missing_values) > 0
            ) {
                median(non_missing_values)
            } else {
                NA_real_
            },
            Min_Beta = if (
                length(non_missing_values) > 0
            ) {
                min(non_missing_values)
            } else {
                NA_real_
            },
            Max_Beta = if (
                length(non_missing_values) > 0
            ) {
                max(non_missing_values)
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }
)

probe_qc <- dplyr::bind_rows(
    probe_qc
)


# -----------------------------------------------------------------------------
# 6. Calculate sample-level 11-probe promoter methylation mean
#
# The mean is calculated using all non-missing values among the 11 probes.
# Therefore, a sample can still receive a promoter mean if one or more probe
# values are missing. The number of available probes is retained for QC.
# -----------------------------------------------------------------------------

promoter_matrix <- as.matrix(
    master[
        ,
        ..promoter_probes_11
    ]
)

storage.mode(
    promoter_matrix
) <- "numeric"

n_available_per_sample <- rowSums(
    !is.na(promoter_matrix)
)

promoter_mean <- rowMeans(
    promoter_matrix,
    na.rm = TRUE
)

promoter_mean[
    n_available_per_sample == 0
] <- NA_real_


# -----------------------------------------------------------------------------
# 7. Add the new promoter-mean variables after the methylation probe block
# -----------------------------------------------------------------------------

all_cpg_columns <- grep(
    "^cg[0-9]+$",
    names(master),
    value = TRUE
)

if (length(all_cpg_columns) == 0) {
    stop(
        "No methylation CpG columns were detected in the master file."
    )
}

last_cpg_column <- tail(
    all_cpg_columns,
    1
)

master$DNMT3A_Promoter_Mean_Methylation_11Probe <- promoter_mean

master$DNMT3A_Promoter_Probes_Available_N <- n_available_per_sample

master$DNMT3A_Promoter_Probes_Missing_N <- (
    length(promoter_probes_11) -
    n_available_per_sample
)

master <- master %>%
    dplyr::relocate(
        DNMT3A_Promoter_Mean_Methylation_11Probe,
        DNMT3A_Promoter_Probes_Available_N,
        DNMT3A_Promoter_Probes_Missing_N,
        .after = dplyr::all_of(last_cpg_column)
    )


# -----------------------------------------------------------------------------
# 8. Sample-level QC
# -----------------------------------------------------------------------------

sample_id_columns <- intersect(
    c(
        "Patient_ID",
        "Sample_Core",
        "Sample_ID"
    ),
    names(master)
)

sample_qc <- master %>%
    dplyr::select(
        dplyr::all_of(sample_id_columns),
        dplyr::all_of(promoter_probes_11),
        DNMT3A_Promoter_Mean_Methylation_11Probe,
        DNMT3A_Promoter_Probes_Available_N,
        DNMT3A_Promoter_Probes_Missing_N
    ) %>%
    dplyr::mutate(
        Promoter_Mean_QC = dplyr::case_when(
            DNMT3A_Promoter_Probes_Available_N == 11 ~
                "Complete: 11/11 probes",

            DNMT3A_Promoter_Probes_Available_N >= 9 ~
                "Partial: 9-10/11 probes",

            DNMT3A_Promoter_Probes_Available_N > 0 ~
                "Low completeness: <9/11 probes",

            TRUE ~
                "No evaluable promoter probes"
        )
    )


# -----------------------------------------------------------------------------
# 9. Summary QC
# -----------------------------------------------------------------------------

mean_values <- master$DNMT3A_Promoter_Mean_Methylation_11Probe

summary_qc <- data.frame(
    Metric = c(
        "Total samples",
        "Promoter probes defined",
        "Promoter probes used for mean",
        "Samples with promoter mean",
        "Samples missing promoter mean",
        "Samples with 11/11 probes",
        "Samples with 10/11 probes",
        "Samples with 9/11 probes",
        "Samples with fewer than 9/11 probes",
        "Mean of promoter mean beta values",
        "Median of promoter mean beta values",
        "Minimum promoter mean beta value",
        "Maximum promoter mean beta value"
    ),

    Value = c(
        nrow(master),
        length(promoter_probes_11),
        length(promoter_probes_11),
        sum(!is.na(mean_values)),
        sum(is.na(mean_values)),
        sum(n_available_per_sample == 11),
        sum(n_available_per_sample == 10),
        sum(n_available_per_sample == 9),
        sum(n_available_per_sample < 9),
        ifelse(
            all(is.na(mean_values)),
            NA,
            mean(
                mean_values,
                na.rm = TRUE
            )
        ),
        ifelse(
            all(is.na(mean_values)),
            NA,
            median(
                mean_values,
                na.rm = TRUE
            )
        ),
        ifelse(
            all(is.na(mean_values)),
            NA,
            min(
                mean_values,
                na.rm = TRUE
            )
        ),
        ifelse(
            all(is.na(mean_values)),
            NA,
            max(
                mean_values,
                na.rm = TRUE
            )
        )
    ),

    stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 10. Validate beta-value range
# -----------------------------------------------------------------------------

out_of_range_probe_values <- 0

for (probe in promoter_probes_11) {

    values <- master[[probe]]

    out_of_range_probe_values <- (
        out_of_range_probe_values +
        sum(
            values < 0 |
            values > 1,
            na.rm = TRUE
        )
    )
}

out_of_range_mean_values <- sum(
    mean_values < 0 |
    mean_values > 1,
    na.rm = TRUE
)

range_qc <- data.frame(
    QC_Item = c(
        "Individual promoter probe beta values outside 0-1",
        "11-probe promoter mean values outside 0-1"
    ),
    N = c(
        out_of_range_probe_values,
        out_of_range_mean_values
    ),
    stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 11. Save output files
# -----------------------------------------------------------------------------

output_dir <- dirname(
    normalizePath(
        input_file,
        winslash = "/",
        mustWork = TRUE
    )
)

master_output <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_MASTER_WITH_11PROBE_PROMOTER_MEAN.csv"
)

sample_qc_output <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_11PROBE_PROMOTER_MEAN_SAMPLE_QC.csv"
)

probe_qc_output <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_11PROBE_PROMOTER_MEAN_PROBE_QC.csv"
)

summary_qc_output <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_11PROBE_PROMOTER_MEAN_SUMMARY_QC.csv"
)

range_qc_output <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_11PROBE_PROMOTER_MEAN_RANGE_QC.csv"
)

data.table::fwrite(
    master,
    master_output,
    na = ""
)

data.table::fwrite(
    sample_qc,
    sample_qc_output,
    na = ""
)

data.table::fwrite(
    probe_qc,
    probe_qc_output,
    na = ""
)

data.table::fwrite(
    summary_qc,
    summary_qc_output,
    na = ""
)

data.table::fwrite(
    range_qc,
    range_qc_output,
    na = ""
)


# -----------------------------------------------------------------------------
# 12. Console QC report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("DNMT3A 11-PROBE PROMOTER MEAN COMPLETED\n")
cat("============================================================\n")

cat("\n11 probes used:\n")
print(
    promoter_probes_11
)

cat("\nProbe-level completeness:\n")
print(
    probe_qc[
        ,
        c(
            "Probe",
            "Included_in_11Probe_Mean",
            "N_Available",
            "N_Missing",
            "Completeness_Percent"
        )
    ]
)

cat("\nNumber of evaluable promoter probes per sample:\n")
print(
    table(
        n_available_per_sample,
        useNA = "ifany"
    )
)

cat("\nPromoter mean summary:\n")
print(
    summary(
        mean_values
    )
)

cat("\nBeta-value range QC:\n")
print(
    range_qc
)

cat("\nUpdated master file:\n")
cat(
    master_output,
    "\n"
)

cat("\nSample-level QC file:\n")
cat(
    sample_qc_output,
    "\n"
)

cat("\nProbe-level QC file:\n")
cat(
    probe_qc_output,
    "\n"
)

cat("\nSummary QC file:\n")
cat(
    summary_qc_output,
    "\n"
)

cat("\nRange QC file:\n")
cat(
    range_qc_output,
    "\n"
)

cat(
    "\nMain downstream methylation variable:\n",
    "DNMT3A_Promoter_Mean_Methylation_11Probe\n",
    sep = ""
)
