# =============================================================================
# TCGA-LGG DNMT3A Multi-Omics Integration
#
# Integrates:
#   1. TCGA-LGG legacy HiSeqV2 RNA-seq from UCSC Xena
#      - annotated expression matrix
#      - values are already log2-transformed normalized expression
#   2. TCGA-LGG legacy HumanMethylation450 DNMT3A methylation
#   3. TCGAbiolinks LGG subtype annotation
#   4. PanCanAtlas survival endpoints
#   5. TCGA-LGG GISTIC2 thresholded CNA for CDKN2A, CDKN2B, and EGFR
#
# Important barcode rule:
#   - Legacy HiSeqV2 RNA sample columns are 15-character sample barcodes
#     such as TCGA-XX-XXXX-01 and therefore do NOT contain a vial character.
#   - RNA and methylation are matched using Sample_Core (first 15 characters).
#   - The 16-character methylation barcode is retained as Sample_ID when present.
#   - Patient_ID = first 12 characters.
#   - Primary tumor sample type = 01.
#
# No expression or methylation transformation is performed in this script.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c("data.table", "dplyr", "stringr")

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

library(data.table)
library(dplyr)
library(stringr)


# -----------------------------------------------------------------------------
# 1. Select input files
# -----------------------------------------------------------------------------

cat("\nSelect the annotated TCGA-LGG legacy HiSeqV2 RNA-seq CSV file.\n")
rna_file <- file.choose()

cat("\nSelect the TCGA-LGG legacy HumanMethylation450 DNMT3A methylation TSV file.\n")
methylation_file <- file.choose()

cat("\nSelect the TCGAbiolinks TCGA-LGG subtype annotation CSV file.\n")
subtype_file <- file.choose()

cat("\nSelect the PanCanAtlas TCGA-LGG clinical/survival TSV file.\n")
survival_file <- file.choose()

cat("\nSelect the TCGA-LGG GISTIC2 thresholded-by-genes CNA file.\n")
cna_file <- file.choose()


# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

normalize_tcga <- function(x) {
    x <- as.character(x)
    x <- toupper(trimws(x))
    x[x == ""] <- NA_character_
    x
}

patient_id <- function(x) {
    x <- normalize_tcga(x)
    ifelse(!is.na(x) & nchar(x) >= 12, substr(x, 1, 12), NA_character_)
}

sample_id_16 <- function(x) {
    x <- normalize_tcga(x)
    ifelse(!is.na(x) & nchar(x) >= 16, substr(x, 1, 16), NA_character_)
}

sample_core_15 <- function(x) {
    x <- normalize_tcga(x)
    ifelse(!is.na(x) & nchar(x) >= 15, substr(x, 1, 15), NA_character_)
}

sample_type_code <- function(x) {
    x <- normalize_tcga(x)
    ifelse(!is.na(x) & nchar(x) >= 15, substr(x, 14, 15), NA_character_)
}


# -----------------------------------------------------------------------------
# 3. RNA-seq: extract DNMT3A expression and define primary tumor cohort
# -----------------------------------------------------------------------------

cat("\nReading RNA-seq header...\n")

rna_header <- names(fread(rna_file, nrows = 0))
rna_sample_cols <- rna_header[grepl("^TCGA-", rna_header, ignore.case = TRUE)]

if (!"Gene_Symbol" %in% rna_header) {
    stop("RNA-seq file must contain a 'Gene_Symbol' column.")
}

if (length(rna_sample_cols) == 0) {
    stop("No TCGA sample columns were detected in the RNA-seq file.")
}

cat("Extracting DNMT3A expression...\n")

rna_dnmt3a <- fread(
    rna_file,
    select = c("Gene_Symbol", rna_sample_cols)
)

rna_dnmt3a <- rna_dnmt3a[
    toupper(trimws(as.character(Gene_Symbol))) == "DNMT3A"
]

if (nrow(rna_dnmt3a) != 1) {
    stop(
        paste0(
            "Expected exactly one DNMT3A row in RNA-seq, but found ",
            nrow(rna_dnmt3a),
            "."
        )
    )
}

rna_long <- data.table(
    Raw_RNA_Sample = rna_sample_cols,
    DNMT3A_expression = as.numeric(unlist(rna_dnmt3a[1, ..rna_sample_cols]))
)

rna_long[, Raw_RNA_Sample := normalize_tcga(Raw_RNA_Sample)]
rna_long[, Sample_Core := sample_core_15(Raw_RNA_Sample)]
rna_long[, Patient_ID := patient_id(Raw_RNA_Sample)]
rna_long[, Sample_Type := sample_type_code(Raw_RNA_Sample)]

# Legacy HiSeqV2 uses 15-character sample barcodes, so Sample_Core is the
# correct RNA sample-level key. Do not require a 16-character Sample_ID here.
rna_primary <- rna_long[
    Sample_Type == "01" & !is.na(Sample_Core)
]

if (anyDuplicated(rna_primary$Sample_Core) > 0) {
    stop("Duplicate primary RNA Sample_Core values were detected.")
}

if (anyDuplicated(rna_primary$Patient_ID) > 0) {
    warning("More than one primary RNA sample exists for at least one Patient_ID.")
}

cat(
    "Legacy HiSeqV2 RNA-seq: ",
    nrow(rna_long),
    " total samples; ",
    nrow(rna_primary),
    " primary tumor samples.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 4. DNA methylation: retain all DNMT3A CpGs
# -----------------------------------------------------------------------------

cat("\nReading DNMT3A methylation data...\n")

meth <- fread(methylation_file)

if (!"sample" %in% names(meth)) {
    stop("Methylation file must contain a 'sample' column.")
}

meth[, Raw_Methylation_Sample := normalize_tcga(sample)]
meth[, Sample_ID := sample_id_16(Raw_Methylation_Sample)]
meth[, Sample_Core := sample_core_15(Raw_Methylation_Sample)]
meth[, Patient_ID := patient_id(Raw_Methylation_Sample)]
meth[, Sample_Type := sample_type_code(Raw_Methylation_Sample)]

meth_primary <- meth[
    Sample_Type == "01" & !is.na(Sample_Core)
]

if (anyDuplicated(meth_primary$Sample_Core) > 0) {
    duplicate_meth <- meth_primary[
        duplicated(Sample_Core) | duplicated(Sample_Core, fromLast = TRUE),
        .(Raw_Methylation_Sample, Sample_ID, Sample_Core, Patient_ID)
    ]
    print(duplicate_meth)
    stop(
        "Duplicate primary methylation Sample_Core values were detected. ",
        "The RNA-to-methylation merge would be ambiguous."
    )
}

meth_meta_cols <- c(
    "sample",
    "Raw_Methylation_Sample",
    "Sample_ID",
    "Sample_Core",
    "Patient_ID",
    "Sample_Type"
)

meth_cpg_cols <- setdiff(names(meth_primary), meth_meta_cols)

if (length(meth_cpg_cols) == 0) {
    stop("No methylation CpG columns were detected.")
}

if (!"cg21629895" %in% meth_cpg_cols) {
    warning("cg21629895 was not detected in the methylation file.")
}

# Convert CpG columns to numeric without transforming beta values.
for (cpg in meth_cpg_cols) {
    meth_primary[[cpg]] <- suppressWarnings(as.numeric(meth_primary[[cpg]]))
}

# Basic beta-value QC.
beta_values <- unlist(meth_primary[, ..meth_cpg_cols], use.names = FALSE)
beta_values <- beta_values[!is.na(beta_values)]

if (length(beta_values) > 0) {
    beta_out_of_range <- sum(beta_values < 0 | beta_values > 1)

    if (beta_out_of_range > 0) {
        warning(
            beta_out_of_range,
            " non-missing methylation values fall outside the expected beta-value range [0, 1]."
        )
    }
}

# Match RNA and methylation by the common 15-character sample barcode.
# Retain the vial-level methylation Sample_ID as the final Sample_ID.
meth_for_merge <- meth_primary %>%
    as.data.frame() %>%
    dplyr::select(
        Sample_Core,
        Sample_ID,
        all_of(meth_cpg_cols)
    )

cat(
    "Legacy HumanMethylation450: ",
    nrow(meth),
    " total samples; ",
    nrow(meth_primary),
    " primary tumor samples; ",
    length(meth_cpg_cols),
    " CpG columns retained.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 5. TCGAbiolinks subtype annotation
# -----------------------------------------------------------------------------

cat("\nReading TCGAbiolinks subtype annotation...\n")

subtype <- fread(subtype_file)

if (!"patient" %in% names(subtype)) {
    stop("TCGAbiolinks subtype file must contain a 'patient' column.")
}

subtype[, Patient_ID := patient_id(patient)]

if (anyDuplicated(subtype$Patient_ID) > 0) {
    stop("Duplicate Patient_ID values were detected in the TCGAbiolinks subtype file.")
}

subtype_for_merge <- subtype %>%
    as.data.frame() %>%
    dplyr::select(Patient_ID, dplyr::everything())

cat(
    "TCGAbiolinks subtype annotation: ",
    nrow(subtype_for_merge),
    " patients; ",
    ncol(subtype_for_merge) - 1,
    " original annotation columns retained.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 6. PanCanAtlas survival endpoints
# -----------------------------------------------------------------------------

cat("\nReading PanCanAtlas survival data...\n")

surv <- fread(survival_file)

required_survival_cols <- c(
    "Patient ID",
    "Sample ID",
    "Overall Survival (Months)",
    "Overall Survival Status",
    "Progress Free Survival (Months)",
    "Progression Free Status",
    "Months of disease-specific survival",
    "Disease-specific Survival status"
)

missing_survival_cols <- setdiff(required_survival_cols, names(surv))

if (length(missing_survival_cols) > 0) {
    stop(
        paste(
            "Missing required PanCanAtlas survival columns:",
            paste(missing_survival_cols, collapse = ", ")
        )
    )
}

surv_for_merge <- surv %>%
    as.data.frame() %>%
    dplyr::transmute(
        Patient_ID = patient_id(`Patient ID`),
        PanCan_Sample_ID = normalize_tcga(`Sample ID`),
        OS_months = suppressWarnings(as.numeric(`Overall Survival (Months)`)),
        OS_status_raw = as.character(`Overall Survival Status`),
        PFS_months = suppressWarnings(as.numeric(`Progress Free Survival (Months)`)),
        PFS_status_raw = as.character(`Progression Free Status`),
        DSS_months = suppressWarnings(as.numeric(`Months of disease-specific survival`)),
        DSS_status_raw = as.character(`Disease-specific Survival status`)
    ) %>%
    dplyr::mutate(
        OS_status = case_when(
            str_detect(OS_status_raw, "^1:") ~ 1L,
            str_detect(OS_status_raw, "^0:") ~ 0L,
            TRUE ~ NA_integer_
        ),
        PFS_status = case_when(
            str_detect(PFS_status_raw, "^1:") ~ 1L,
            str_detect(PFS_status_raw, "^0:") ~ 0L,
            TRUE ~ NA_integer_
        ),
        DSS_status = case_when(
            str_detect(DSS_status_raw, "^1:") ~ 1L,
            str_detect(DSS_status_raw, "^0:") ~ 0L,
            TRUE ~ NA_integer_
        )
    )

if (anyDuplicated(surv_for_merge$Patient_ID) > 0) {
    stop("Duplicate Patient_ID values were detected in the PanCanAtlas survival file.")
}

cat(
    "PanCanAtlas survival: ",
    nrow(surv_for_merge),
    " patients.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 7. GISTIC2 CNA: extract CDKN2A, CDKN2B, and EGFR
# -----------------------------------------------------------------------------

cat("\nReading GISTIC2 CNA data...\n")

cna <- fread(cna_file)

if (!"Gene Symbol" %in% names(cna)) {
    stop("GISTIC2 CNA file must contain a 'Gene Symbol' column.")
}

cna_genes <- c("CDKN2A", "CDKN2B", "EGFR")
cna_selected <- cna[`Gene Symbol` %in% cna_genes]

missing_cna_genes <- setdiff(cna_genes, cna_selected$`Gene Symbol`)

if (length(missing_cna_genes) > 0) {
    stop(
        paste(
            "The following CNA genes were not found:",
            paste(missing_cna_genes, collapse = ", ")
        )
    )
}

if (anyDuplicated(cna_selected$`Gene Symbol`) > 0) {
    stop("Duplicate gene rows were detected for CDKN2A/CDKN2B/EGFR in the CNA file.")
}

cna_long <- melt(
    cna_selected,
    id.vars = "Gene Symbol",
    variable.name = "Raw_CNA_Sample",
    value.name = "GISTIC"
)

cna_long[, Sample_Core := sample_core_15(Raw_CNA_Sample)]
cna_long[, Patient_ID := patient_id(Raw_CNA_Sample)]

# Remove non-TCGA/non-sample columns if present.
cna_long <- cna_long[
    !is.na(Sample_Core) & grepl("^TCGA-", normalize_tcga(Raw_CNA_Sample))
]

cna_wide <- dcast(
    cna_long,
    Sample_Core + Patient_ID ~ `Gene Symbol`,
    value.var = "GISTIC"
)

setnames(
    cna_wide,
    old = c("CDKN2A", "CDKN2B", "EGFR"),
    new = c("CDKN2A_GISTIC", "CDKN2B_GISTIC", "EGFR_GISTIC")
)

cna_wide[, CDKN2A_GISTIC := suppressWarnings(as.numeric(CDKN2A_GISTIC))]
cna_wide[, CDKN2B_GISTIC := suppressWarnings(as.numeric(CDKN2B_GISTIC))]
cna_wide[, EGFR_GISTIC := suppressWarnings(as.numeric(EGFR_GISTIC))]

cna_wide[, CDKN2A_HomDel := fifelse(
    is.na(CDKN2A_GISTIC), NA_integer_,
    fifelse(CDKN2A_GISTIC == -2, 1L, 0L)
)]

cna_wide[, CDKN2B_HomDel := fifelse(
    is.na(CDKN2B_GISTIC), NA_integer_,
    fifelse(CDKN2B_GISTIC == -2, 1L, 0L)
)]

cna_wide[, EGFR_Amplification := fifelse(
    is.na(EGFR_GISTIC), NA_integer_,
    fifelse(EGFR_GISTIC == 2, 1L, 0L)
)]

cna_for_merge <- cna_wide %>%
    as.data.frame() %>%
    dplyr::select(
        Sample_Core,
        CDKN2A_GISTIC,
        CDKN2B_GISTIC,
        EGFR_GISTIC,
        CDKN2A_HomDel,
        CDKN2B_HomDel,
        EGFR_Amplification
    )

if (anyDuplicated(cna_for_merge$Sample_Core) > 0) {
    stop("Duplicate Sample_Core values were detected in the GISTIC2 CNA data.")
}

cat(
    "GISTIC2 CNA: ",
    nrow(cna_for_merge),
    " samples with CDKN2A/CDKN2B/EGFR values.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 8. Build the primary-tumor master table
# -----------------------------------------------------------------------------

cat("\nBuilding the integrated master table...\n")

master <- rna_primary %>%
    as.data.frame() %>%
    dplyr::select(
        Patient_ID,
        Sample_Core,
        DNMT3A_expression
    ) %>%
    dplyr::left_join(
        meth_for_merge,
        by = "Sample_Core"
    ) %>%
    dplyr::left_join(
        subtype_for_merge,
        by = "Patient_ID"
    ) %>%
    dplyr::left_join(
        surv_for_merge,
        by = "Patient_ID"
    ) %>%
    dplyr::left_join(
        cna_for_merge,
        by = "Sample_Core"
    ) %>%
    dplyr::relocate(
        Patient_ID,
        Sample_ID,
        Sample_Core,
        DNMT3A_expression
    )

if (nrow(master) != nrow(rna_primary)) {
    stop(
        paste0(
            "Integration changed the row count from ",
            nrow(rna_primary),
            " to ",
            nrow(master),
            ". Check duplicate merge keys."
        )
    )
}

if (anyDuplicated(master$Sample_Core) > 0) {
    stop("Duplicate Sample_Core values were created in the integrated master table.")
}

if (anyDuplicated(master$Patient_ID) > 0) {
    warning(
        "More than one primary sample exists for at least one Patient_ID. ",
        "Review the duplicate-patient QC file before patient-level analyses."
    )
}


# -----------------------------------------------------------------------------
# 9. QC summary
# -----------------------------------------------------------------------------
# Explicit dplyr:: prefixes are used below to prevent namespace conflicts with
# packages such as matrixStats/Bioconductor that may mask dplyr::count().

cat("\nGenerating integration QC...\n")

rna_meth_overlap <- sum(
    rna_primary$Sample_Core %in% meth_primary$Sample_Core
)

qc_summary <- data.frame(
    Metric = c(
        "RNA total samples",
        "RNA primary samples",
        "RNA primary unique Patient_ID",
        "RNA primary unique Sample_Core",
        "Methylation total samples",
        "Methylation primary samples",
        "RNA-methylation primary Sample_Core overlap",
        "TCGAbiolinks patients",
        "TCGAbiolinks overlap with master",
        "PanCanAtlas survival patients",
        "PanCanAtlas overlap with master",
        "GISTIC2 selected-gene samples",
        "GISTIC2 Sample_Core overlap with master",
        "Final master rows",
        "Final unique Sample_Core",
        "Final unique Patient_ID",
        "Final non-missing vial-level Sample_ID",
        "DNMT3A expression non-missing",
        "cg21629895 non-missing",
        "OS non-missing",
        "PFS non-missing",
        "DSS non-missing",
        "CDKN2A GISTIC non-missing",
        "CDKN2B GISTIC non-missing",
        "EGFR GISTIC non-missing"
    ),
    N = c(
        nrow(rna_long),
        nrow(rna_primary),
        dplyr::n_distinct(rna_primary$Patient_ID),
        dplyr::n_distinct(rna_primary$Sample_Core),
        nrow(meth),
        nrow(meth_primary),
        rna_meth_overlap,
        nrow(subtype_for_merge),
        sum(unique(master$Patient_ID) %in% subtype_for_merge$Patient_ID),
        nrow(surv_for_merge),
        sum(unique(master$Patient_ID) %in% surv_for_merge$Patient_ID),
        nrow(cna_for_merge),
        sum(master$Sample_Core %in% cna_for_merge$Sample_Core),
        nrow(master),
        dplyr::n_distinct(master$Sample_Core),
        dplyr::n_distinct(master$Patient_ID),
        sum(!is.na(master$Sample_ID)),
        sum(!is.na(master$DNMT3A_expression)),
        if ("cg21629895" %in% names(master)) sum(!is.na(master$cg21629895)) else NA_integer_,
        sum(!is.na(master$OS_months)),
        sum(!is.na(master$PFS_months)),
        sum(!is.na(master$DSS_months)),
        sum(!is.na(master$CDKN2A_GISTIC)),
        sum(!is.na(master$CDKN2B_GISTIC)),
        sum(!is.na(master$EGFR_GISTIC))
    ),
    stringsAsFactors = FALSE
)

duplicate_patients <- master %>%
    dplyr::count(Patient_ID, name = "N_primary_samples") %>%
    dplyr::filter(N_primary_samples > 1)

unmatched_summary <- data.frame(
    Patient_ID = master$Patient_ID,
    Sample_ID = master$Sample_ID,
    Sample_Core = master$Sample_Core,
    Missing_Methylation = !master$Sample_Core %in% meth_primary$Sample_Core,
    Missing_TCGAbiolinks = !master$Patient_ID %in% subtype_for_merge$Patient_ID,
    Missing_PanCanAtlas_Survival = !master$Patient_ID %in% surv_for_merge$Patient_ID,
    Missing_GISTIC2_CNA = !master$Sample_Core %in% cna_for_merge$Sample_Core,
    stringsAsFactors = FALSE
)

unmatched_summary <- unmatched_summary %>%
    dplyr::filter(
        Missing_Methylation |
        Missing_TCGAbiolinks |
        Missing_PanCanAtlas_Survival |
        Missing_GISTIC2_CNA
    )

cpg_missingness <- data.frame(
    CpG = meth_cpg_cols,
    N_non_missing = vapply(
        meth_cpg_cols,
        function(x) sum(!is.na(master[[x]])),
        FUN.VALUE = numeric(1)
    ),
    N_missing = vapply(
        meth_cpg_cols,
        function(x) sum(is.na(master[[x]])),
        FUN.VALUE = numeric(1)
    ),
    stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 10. Save outputs
# -----------------------------------------------------------------------------

output_dir <- dirname(normalizePath(rna_file))

master_file <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_PRIMARY_MULTIOMICS_MASTER_LEGACY.csv"
)

qc_file <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_INTEGRATION_QC_LEGACY.csv"
)

unmatched_file <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_UNMATCHED_SAMPLES_LEGACY.csv"
)

duplicate_file <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_DUPLICATE_PATIENTS_LEGACY.csv"
)

cpg_qc_file <- file.path(
    output_dir,
    "TCGA_LGG_DNMT3A_CPG_MISSINGNESS_QC_LEGACY.csv"
)

fwrite(master, master_file)
fwrite(qc_summary, qc_file)
fwrite(unmatched_summary, unmatched_file)
fwrite(duplicate_patients, duplicate_file)
fwrite(cpg_missingness, cpg_qc_file)


# -----------------------------------------------------------------------------
# 11. Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("INTEGRATION COMPLETED\n")
cat("============================================================\n\n")

print(qc_summary)

cat("\nFinal master dimensions: ")
cat(nrow(master), " rows x ", ncol(master), " columns\n", sep = "")

cat("\nKey barcode rule used:\n")
cat("RNA HiSeqV2 was matched to methylation by 15-character Sample_Core.\n")
cat("The 16-character methylation barcode was retained as Sample_ID when available.\n")

cat("\nAll TCGAbiolinks annotation columns were retained.\n")
cat("All DNMT3A methylation CpG columns were retained.\n")
cat("DNMT3A expression was retained as a continuous value without transformation.\n")
cat("No expression or methylation stratification was performed.\n")

cat("\nOutput files:\n")
cat(master_file, "\n")
cat(qc_file, "\n")
cat(unmatched_file, "\n")
cat(duplicate_file, "\n")
cat(cpg_qc_file, "\n")

cat("\nReview the QC files before starting downstream analyses.\n")
