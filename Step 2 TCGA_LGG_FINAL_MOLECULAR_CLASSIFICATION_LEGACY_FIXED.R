# =============================================================================
# TCGA-LGG Molecular Classification
#
# This script adds TWO parallel classifications:
#
# A. Main three molecular groups:
#    1. IDH-mutant/1p19q-codeleted
#    2. IDH-mutant/1p19q-non-codeleted
#    3. IDH-wildtype
#
# B. Strict WHO 2021-aligned classification:
#    1. Oligodendroglioma, IDH-mutant and 1p/19q-codeleted
#    2. Astrocytoma, IDH-mutant
#    3. Glioblastoma, IDH-wildtype
#    4. IDH-wildtype, no GBM-defining molecular feature detected
#    5. IDH-wildtype, incomplete molecular information
#
# For IDH-wildtype tumors, strict molecular GBM classification requires at
# least one available GBM-defining molecular feature:
#    - TERT promoter mutation
#    - EGFR amplification
#    - Combined chromosome 7 gain / chromosome 10 loss
#
# The main downstream analysis variable is Molecular_3Group.
# Strict_WHO2021_Classification is retained as a parallel annotation.
#
# Current legacy multi-omics master:
#   - Sample_Core (15-character TCGA sample barcode) is the primary sample key.
#   - Sample_ID may be unavailable/NA and is therefore not required.
# =============================================================================

required_packages <- c("data.table", "dplyr")

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

library(data.table)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Select input file
# -----------------------------------------------------------------------------

cat("\nSelect the TCGA-LGG integrated master CSV file.\n")
input_file <- file.choose()

master <- fread(input_file)

# -----------------------------------------------------------------------------
# 2. Confirm required columns
# -----------------------------------------------------------------------------

required_columns <- c(
    "Patient_ID",
    "Sample_Core",
    "IDH.codel.subtype",
    "TERT.promoter.status",
    "Chr.7.gain.Chr.10.loss",
    "EGFR_Amplification"
)

missing_columns <- setdiff(required_columns, names(master))

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing required columns:",
            paste(missing_columns, collapse = ", ")
        )
    )
}

# -----------------------------------------------------------------------------
# 3. Standardize source variables
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        IDH_codel_source = trimws(as.character(`IDH.codel.subtype`)),
        TERT_source = trimws(as.character(`TERT.promoter.status`)),
        Chr7_10_source = trimws(as.character(`Chr.7.gain.Chr.10.loss`)),
        EGFR_amp_source = suppressWarnings(as.numeric(`EGFR_Amplification`))
    )

# -----------------------------------------------------------------------------
# 4. Main three molecular groups
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        Molecular_3Group = case_when(
            IDH_codel_source == "IDHmut-codel" ~
                "IDH-mutant/1p19q-codeleted",

            IDH_codel_source == "IDHmut-non-codel" ~
                "IDH-mutant/1p19q-non-codeleted",

            IDH_codel_source == "IDHwt" ~
                "IDH-wildtype",

            TRUE ~ NA_character_
        ),

        Molecular_3Group_Factor = factor(
            Molecular_3Group,
            levels = c(
                "IDH-mutant/1p19q-codeleted",
                "IDH-mutant/1p19q-non-codeleted",
                "IDH-wildtype"
            )
        )
    )

# -----------------------------------------------------------------------------
# 5. WHO 2021 molecular feature indicators
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        WHO2021_TERT_Mutant = case_when(
            TERT_source == "Mutant" ~ 1L,
            TERT_source == "WT" ~ 0L,
            TRUE ~ NA_integer_
        ),

        WHO2021_EGFR_Amplified = case_when(
            EGFR_amp_source == 1 ~ 1L,
            EGFR_amp_source == 0 ~ 0L,
            TRUE ~ NA_integer_
        ),

        WHO2021_Chr7Gain_Chr10Loss = case_when(
            Chr7_10_source == "Gain chr 7 & loss chr 10" ~ 1L,
            Chr7_10_source == "No combined CNA" ~ 0L,
            TRUE ~ NA_integer_
        )
    )

# -----------------------------------------------------------------------------
# 6. Determine GBM-defining molecular feature status in IDH-wildtype tumors
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        WHO2021_GBM_Molecular_Feature = case_when(
            Molecular_3Group == "IDH-wildtype" &
            (
                WHO2021_TERT_Mutant == 1 |
                WHO2021_EGFR_Amplified == 1 |
                WHO2021_Chr7Gain_Chr10Loss == 1
            ) ~ 1L,

            Molecular_3Group == "IDH-wildtype" &
            WHO2021_TERT_Mutant == 0 &
            WHO2021_EGFR_Amplified == 0 &
            WHO2021_Chr7Gain_Chr10Loss == 0 ~ 0L,

            Molecular_3Group == "IDH-wildtype" ~ NA_integer_,

            TRUE ~ NA_integer_
        ),

        WHO2021_GBM_Criteria = case_when(
            Molecular_3Group != "IDH-wildtype" ~ NA_character_,

            WHO2021_TERT_Mutant == 1 &
            WHO2021_EGFR_Amplified == 1 &
            WHO2021_Chr7Gain_Chr10Loss == 1 ~
                "TERT promoter mutation; EGFR amplification; +7/-10",

            WHO2021_TERT_Mutant == 1 &
            WHO2021_EGFR_Amplified == 1 ~
                "TERT promoter mutation; EGFR amplification",

            WHO2021_TERT_Mutant == 1 &
            WHO2021_Chr7Gain_Chr10Loss == 1 ~
                "TERT promoter mutation; +7/-10",

            WHO2021_EGFR_Amplified == 1 &
            WHO2021_Chr7Gain_Chr10Loss == 1 ~
                "EGFR amplification; +7/-10",

            WHO2021_TERT_Mutant == 1 ~
                "TERT promoter mutation",

            WHO2021_EGFR_Amplified == 1 ~
                "EGFR amplification",

            WHO2021_Chr7Gain_Chr10Loss == 1 ~
                "+7/-10",

            WHO2021_TERT_Mutant == 0 &
            WHO2021_EGFR_Amplified == 0 &
            WHO2021_Chr7Gain_Chr10Loss == 0 ~
                "No GBM-defining molecular feature detected",

            TRUE ~
                "Incomplete molecular information"
        )
    )

# -----------------------------------------------------------------------------
# 7. Strict WHO 2021-aligned classification
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        Strict_WHO2021_Classification = case_when(
            IDH_codel_source == "IDHmut-codel" ~
                "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted",

            IDH_codel_source == "IDHmut-non-codel" ~
                "Astrocytoma, IDH-mutant",

            IDH_codel_source == "IDHwt" &
            WHO2021_GBM_Molecular_Feature == 1 ~
                "Glioblastoma, IDH-wildtype",

            IDH_codel_source == "IDHwt" &
            WHO2021_GBM_Molecular_Feature == 0 ~
                "IDH-wildtype, no GBM-defining molecular feature detected",

            IDH_codel_source == "IDHwt" &
            is.na(WHO2021_GBM_Molecular_Feature) ~
                "IDH-wildtype, incomplete molecular information",

            TRUE ~ NA_character_
        )
    )

# -----------------------------------------------------------------------------
# 8. Strict WHO 2021 three-group variable
#
# Only cases assignable to one of the three adult-type diffuse glioma
# categories using the available molecular information are retained.
# Unresolved IDH-wildtype cases are NA.
# -----------------------------------------------------------------------------

master <- master %>%
    dplyr::mutate(
        Strict_WHO2021_3Group = case_when(
            Strict_WHO2021_Classification ==
                "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted" ~
                "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted",

            Strict_WHO2021_Classification ==
                "Astrocytoma, IDH-mutant" ~
                "Astrocytoma, IDH-mutant",

            Strict_WHO2021_Classification ==
                "Glioblastoma, IDH-wildtype" ~
                "Glioblastoma, IDH-wildtype",

            TRUE ~ NA_character_
        ),

        Strict_WHO2021_3Group_Factor = factor(
            Strict_WHO2021_3Group,
            levels = c(
                "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted",
                "Astrocytoma, IDH-mutant",
                "Glioblastoma, IDH-wildtype"
            )
        )
    )

# -----------------------------------------------------------------------------
# 9. Optional CDKN2A/B homozygous deletion annotation
# -----------------------------------------------------------------------------

if (
    "CDKN2A_HomDel" %in% names(master) &&
    "CDKN2B_HomDel" %in% names(master)
) {

    master <- master %>%
        mutate(
            WHO2021_CDKN2AB_HomDel = case_when(
                `CDKN2A_HomDel` == 1 | `CDKN2B_HomDel` == 1 ~ 1L,
                `CDKN2A_HomDel` == 0 & `CDKN2B_HomDel` == 0 ~ 0L,
                TRUE ~ NA_integer_
            )
        )
}

# -----------------------------------------------------------------------------
# 10. QC summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MOLECULAR CLASSIFICATION SUMMARY\n")
cat("============================================================\n")

cat("\nSource IDH.codel.subtype:\n")
print(table(master$IDH_codel_source, useNA = "ifany"))

cat("\nMain three molecular groups:\n")
print(table(master$Molecular_3Group, useNA = "ifany"))

cat("\nMain three-group classification completeness:\n")
cat(
    "Classified:",
    sum(!is.na(master$Molecular_3Group)),
    "/",
    nrow(master),
    "\n"
)

cat("\nStrict WHO 2021-aligned classification:\n")
print(table(master$Strict_WHO2021_Classification, useNA = "ifany"))

cat("\nStrict WHO 2021 three-group variable:\n")
print(table(master$Strict_WHO2021_3Group, useNA = "ifany"))

cat("\nIDH-wildtype molecular GBM criteria:\n")
print(
    master %>%
        dplyr::filter(Molecular_3Group == "IDH-wildtype") %>%
        dplyr::count(WHO2021_GBM_Criteria, name = "N") %>%
        dplyr::arrange(dplyr::desc(N))
)

# -----------------------------------------------------------------------------
# 11. Cross-check table
# -----------------------------------------------------------------------------

qc_columns <- c(
    "Patient_ID",
    "Sample_Core",
    "IDH.codel.subtype",
    "Molecular_3Group",
    "TERT.promoter.status",
    "Chr.7.gain.Chr.10.loss",
    "EGFR_Amplification",
    "WHO2021_TERT_Mutant",
    "WHO2021_EGFR_Amplified",
    "WHO2021_Chr7Gain_Chr10Loss",
    "WHO2021_GBM_Molecular_Feature",
    "WHO2021_GBM_Criteria",
    "Strict_WHO2021_Classification",
    "Strict_WHO2021_3Group"
)

if ("Sample_ID" %in% names(master)) {
    qc_columns <- append(
        qc_columns,
        "Sample_ID",
        after = match("Sample_Core", qc_columns)
    )
}

if ("IDH.status" %in% names(master)) {
    qc_columns <- append(
        qc_columns,
        "IDH.status",
        after = match("IDH.codel.subtype", qc_columns)
    )
}

if ("X1p.19q.codeletion" %in% names(master)) {
    insert_after_pos <- if ("IDH.status" %in% qc_columns) {
        match("IDH.status", qc_columns)
    } else {
        match("IDH.codel.subtype", qc_columns)
    }

    qc_columns <- append(
        qc_columns,
        "X1p.19q.codeletion",
        after = insert_after_pos
    )
}

classification_qc <- master %>%
    dplyr::select(dplyr::any_of(qc_columns))

# -----------------------------------------------------------------------------
# 12. Save outputs
# -----------------------------------------------------------------------------

output_dir <- dirname(input_file)

output_master <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_MASTER_LEGACY_FINAL_MOLECULAR_CLASSIFICATION.csv"
)

output_qc <- file.path(
    output_dir,
    "TCGA_LGG_LEGACY_FINAL_MOLECULAR_CLASSIFICATION_QC.csv"
)

fwrite(master, output_master)
fwrite(classification_qc, output_qc)

# -----------------------------------------------------------------------------
# 13. Final report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CLASSIFICATION COMPLETED\n")
cat("============================================================\n")

cat("\nFinal master file:\n")
cat(output_master, "\n")

cat("\nClassification QC file:\n")
cat(output_qc, "\n")

cat(
    "\nPrimary downstream variable: Molecular_3Group\n",
    "  IDH-mutant/1p19q-codeleted\n",
    "  IDH-mutant/1p19q-non-codeleted\n",
    "  IDH-wildtype\n",
    sep = ""
)

cat(
    "\nParallel strict annotation: Strict_WHO2021_Classification\n",
    "Strict three-group variable: Strict_WHO2021_3Group\n",
    sep = ""
)

cat(
    "\nImportant: Molecular_3Group retains all IDH-wildtype tumors as one ",
    "molecular group. Strict_WHO2021_3Group includes only IDH-wildtype tumors ",
    "with at least one available molecular GBM-defining feature.\n",
    sep = ""
)
