
# ==============================================================================
# TCGA-LGG DNMT3A MAIN TABLE
# BUILD DIRECTLY FROM:
#   1) Original DNMT3A master CSV
#   2) MC3 LGG mutation .txt.gz file
#
# IMPORTANT:
#   - Existing median-group assignments in the master file are used exactly
#     as provided. No median split is recalculated.
#   - Selected mutation genes are added directly from the MC3 file.
#   - IDH1 is intentionally excluded from the added mutation rows.
#
# Existing grouping variables used:
#   1) DNMT3A_expression_Median_Group
#   2) DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group
#   3) cg21629895_Median_Group
#
# Added mutation genes:
#   FUBP1, PTEN, NOTCH1, NF1, EGFR, TTN, CIC, ATRX, TP53
#
# Mutation definition:
#   Non-silent coding/splice variants:
#     Missense_Mutation
#     Nonsense_Mutation
#     Frame_Shift_Del
#     Frame_Shift_Ins
#     In_Frame_Del
#     In_Frame_Ins
#     Splice_Site
#     Translation_Start_Site
#     Nonstop_Mutation
#     large deletion
#
# Outputs:
#   - Master_with_Selected_Mutations.csv
#   - Table_1_DNMT3A_Clinical_Molecular_FINAL.docx
#   - Table_1_DNMT3A_Clinical_Molecular_FINAL.html
#   - Table_1_Analysis_Data_FINAL.csv
#   - Mutation_Merge_QC.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Package setup
# ------------------------------------------------------------------------------

packages <- c(
    "data.table",
    "dplyr",
    "gtsummary",
    "flextable",
    "officer",
    "openxlsx"
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
    library(dplyr)
    library(gtsummary)
    library(flextable)
    library(officer)
    library(openxlsx)
})


# ------------------------------------------------------------------------------
# 1. Select the TWO input files
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SELECT TWO INPUT FILES\n")
cat("============================================================\n")

cat("\n1. Select the ORIGINAL DNMT3A master CSV file.\n")
MASTER_FILE <- file.choose()

if (!file.exists(MASTER_FILE)) {
    stop(
        paste0(
            "Master file does not exist: ",
            MASTER_FILE
        )
    )
}

cat("\n2. Select the MC3 LGG mutation .txt.gz file.\n")
MUTATION_FILE <- file.choose()

if (!file.exists(MUTATION_FILE)) {
    stop(
        paste0(
            "Mutation file does not exist: ",
            MUTATION_FILE
        )
    )
}


# ------------------------------------------------------------------------------
# 2. Read the master file
# ------------------------------------------------------------------------------

master <- data.table::fread(
    MASTER_FILE,
    check.names = FALSE,
    data.table = TRUE
)

cat("\nMaster file loaded.\n")
cat("Rows: ", nrow(master), "\n", sep = "")
cat("Columns: ", ncol(master), "\n", sep = "")


# ------------------------------------------------------------------------------
# 3. Read the gzipped mutation file
#
# This version does NOT call an external gzip executable, so it works on
# Windows without installing gzip.
# ------------------------------------------------------------------------------

mutation_connection <- gzfile(
    MUTATION_FILE,
    open = "rt"
)

mutation <- utils::read.delim(
    mutation_connection,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

close(
    mutation_connection
)

mutation <- data.table::as.data.table(
    mutation
)

cat("\nMutation file loaded.\n")
cat("Rows: ", nrow(mutation), "\n", sep = "")

if ("sample" %in% names(mutation)) {
    cat(
        "Unique mutation-file samples: ",
        data.table::uniqueN(
            mutation$sample
        ),
        "\n",
        sep = ""
    )
}


# ------------------------------------------------------------------------------
# 4. Validate required columns
# ------------------------------------------------------------------------------

required_master_columns <- c(
    "Patient_ID",
    "Sample_Core",
    "DNMT3A_expression_Median_Group",
    "DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group",
    "cg21629895_Median_Group",
    "Strict_WHO2021_3Group",
    "Age..years.at.diagnosis.",
    "Gender",
    "CDKN2A_GISTIC",
    "CDKN2B_GISTIC",
    "EGFR_GISTIC",
    "TERT.promoter.status",
    "Chr.7.gain.Chr.10.loss",
    "MGMT.promoter.status",
    "OS_months",
    "OS_status"
)

missing_master_columns <- setdiff(
    required_master_columns,
    names(master)
)

if (length(missing_master_columns) > 0) {
    stop(
        paste0(
            "Master file is missing required columns:\n",
            paste(
                missing_master_columns,
                collapse = "\n"
            )
        )
    )
}


required_mutation_columns <- c(
    "sample",
    "gene",
    "effect"
)

missing_mutation_columns <- setdiff(
    required_mutation_columns,
    names(mutation)
)

if (length(missing_mutation_columns) > 0) {
    stop(
        paste0(
            "Mutation file is missing required columns:\n",
            paste(
                missing_mutation_columns,
                collapse = "\n"
            )
        )
    )
}


# ------------------------------------------------------------------------------
# 5. Confirm the existing median groups
# ------------------------------------------------------------------------------

group_columns <- c(
    "DNMT3A_expression_Median_Group",
    "DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group",
    "cg21629895_Median_Group"
)

cat("\nExisting median-group assignments in the master file:\n")

for (group_col in group_columns) {

    cat("\n", group_col, "\n", sep = "")

    print(
        table(
            master[[group_col]],
            useNA = "ifany"
        )
    )
}


# ------------------------------------------------------------------------------
# 6. Define selected mutation genes and mutation effects
# ------------------------------------------------------------------------------

SELECTED_GENES <- c(
    "FUBP1",
    "PTEN",
    "NOTCH1",
    "NF1",
    "EGFR",
    "TTN",
    "CIC",
    "ATRX",
    "TP53"
)

NON_SILENT_EFFECTS <- c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
    "large deletion"
)


# ------------------------------------------------------------------------------
# 7. Build sample-level mutation status
#
# Rules:
#   Mutant  = at least one selected non-silent mutation for that gene
#   WT      = sample is present in MC3 but has no selected non-silent mutation
#   Unknown = sample is absent from the MC3 mutation file entirely
# ------------------------------------------------------------------------------

all_mc3_samples <- unique(
    mutation$sample
)

mutation_selected <- mutation[
    gene %in% SELECTED_GENES &
        effect %in% NON_SILENT_EFFECTS,
    .(
        Mutant = TRUE
    ),
    by = .(
        sample,
        gene
    )
]


for (gene_name in SELECTED_GENES) {

    mutant_samples_gene <- mutation_selected[
        gene == gene_name,
        sample
    ]

    new_column <- paste0(
        gene_name,
        "_Mutation"
    )

    master[
        ,
        (new_column) := data.table::fcase(
            Sample_Core %in% mutant_samples_gene,
            "Mutant",

            Sample_Core %in% all_mc3_samples,
            "WT",

            default = "Unknown"
        )
    ]
}


# ------------------------------------------------------------------------------
# 8. Mutation merge QC
# ------------------------------------------------------------------------------

master_samples <- unique(
    master$Sample_Core
)

matched_samples <- intersect(
    master_samples,
    all_mc3_samples
)

unmatched_samples <- setdiff(
    master_samples,
    all_mc3_samples
)

cat("\nMutation merge QC:\n")
cat(
    "  Master samples: ",
    length(master_samples),
    "\n",
    sep = ""
)
cat(
    "  MC3 samples: ",
    length(all_mc3_samples),
    "\n",
    sep = ""
)
cat(
    "  Matched samples: ",
    length(matched_samples),
    "\n",
    sep = ""
)
cat(
    "  Unmatched master samples: ",
    length(unmatched_samples),
    "\n",
    sep = ""
)


mutation_counts <- data.table::rbindlist(
    lapply(
        SELECTED_GENES,
        function(gene_name) {

            col_name <- paste0(
                gene_name,
                "_Mutation"
            )

            data.table::data.table(
                Gene = gene_name,
                Mutant = sum(
                    master[[col_name]] == "Mutant",
                    na.rm = TRUE
                ),
                WT = sum(
                    master[[col_name]] == "WT",
                    na.rm = TRUE
                ),
                Unknown = sum(
                    master[[col_name]] == "Unknown",
                    na.rm = TRUE
                )
            )
        }
    )
)

print(
    mutation_counts
)


# ------------------------------------------------------------------------------
# 9. Select output folder
#
# The selected folder itself is used directly. No extra nested output folder
# is appended, which avoids Windows/OneDrive "No such file or directory" errors.
# ------------------------------------------------------------------------------

if (.Platform$OS.type == "windows") {

    OUTPUT_DIR <- utils::choose.dir(
        caption = "Select the folder where the final Table 1 files will be saved"
    )

    if (is.na(OUTPUT_DIR) || !nzchar(OUTPUT_DIR)) {
        stop("No output folder was selected.")
    }

} else {

    OUTPUT_DIR <- readline(
        prompt = "Enter the full path to the output folder: "
    )

    if (!nzchar(OUTPUT_DIR)) {
        stop("No output folder was provided.")
    }
}

OUTPUT_DIR <- normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = FALSE
)

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(
        OUTPUT_DIR,
        recursive = TRUE,
        showWarnings = FALSE
    )
}

if (!dir.exists(OUTPUT_DIR)) {
    stop(
        "The selected output directory could not be created or accessed:\n",
        OUTPUT_DIR
    )
}

# Test write permission before creating any output.
write_test_file <- file.path(
    OUTPUT_DIR,
    paste0(".__DNMT3A_TABLE_WRITE_TEST_", Sys.getpid(), ".tmp")
)

write_test_ok <- tryCatch(
    {
        writeLines("write test", write_test_file)
        file.exists(write_test_file)
    },
    error = function(e) FALSE
)

if (file.exists(write_test_file)) {
    unlink(write_test_file)
}

if (!isTRUE(write_test_ok)) {
    stop(
        "R cannot write to the selected output folder:\n",
        OUTPUT_DIR,
        "\nSelect another writable folder, preferably with a shorter path."
    )
}

cat("\nOutput directory confirmed:\n")
cat(OUTPUT_DIR, "\n")

# Check an output file before saving. This also gives a clear message when a
# previous Word/Excel output is still open and therefore locked by Windows.
assert_output_ready <- function(filename) {

    parent_dir <- dirname(filename)

    if (!dir.exists(parent_dir)) {
        dir.create(
            parent_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    if (!dir.exists(parent_dir)) {
        stop(
            "The parent output directory does not exist:\n",
            parent_dir
        )
    }

    if (file.exists(filename)) {

        removed <- tryCatch(
            suppressWarnings(file.remove(filename)),
            error = function(e) FALSE
        )

        if (!isTRUE(removed)) {
            stop(
                "The output file cannot be overwritten:\n",
                filename,
                "\nClose this file in Microsoft Word/Excel and run the script again."
            )
        }
    }

    invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 10. Save the merged master file
# ------------------------------------------------------------------------------

MERGED_MASTER_FILE <- file.path(
    OUTPUT_DIR,
    "Master_with_Selected_Mutations.csv"
)

data.table::fwrite(
    master,
    MERGED_MASTER_FILE
)


# ------------------------------------------------------------------------------
# 11. Prepare variables for the manuscript table
# ------------------------------------------------------------------------------

table_data <- as.data.frame(
    master,
    check.names = FALSE
) %>%
    mutate(

        # Existing median groups are used exactly as provided.
        DNMT3A_mRNA_Group = factor(
            DNMT3A_expression_Median_Group,
            levels = c(
                "High",
                "Low"
            )
        ),

        DNMT3A_Promoter_Methylation_Group = factor(
            DNMT3A_Promoter_Mean_Methylation_11Probe_Median_Group,
            levels = c(
                "High",
                "Low"
            )
        ),

        cg21629895_Group = factor(
            cg21629895_Median_Group,
            levels = c(
                "High",
                "Low"
            )
        ),

        WHO_2021_Classification = factor(
            Strict_WHO2021_3Group,
            levels = c(
                "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted",
                "Astrocytoma, IDH-mutant",
                "Glioblastoma, IDH-wildtype"
            )
        ),

        Age = suppressWarnings(
            as.numeric(
                Age..years.at.diagnosis.
            )
        ),

        Sex = factor(
            case_when(
                tolower(
                    as.character(
                        Gender
                    )
                ) == "female" ~ "Female",

                tolower(
                    as.character(
                        Gender
                    )
                ) == "male" ~ "Male",

                TRUE ~ NA_character_
            ),
            levels = c(
                "Female",
                "Male"
            )
        ),

        CDKN2A_CNA = factor(
            as.character(
                CDKN2A_GISTIC
            ),
            levels = c(
                "-2",
                "-1",
                "0",
                "1",
                "2"
            ),
            labels = c(
                "Deep deletion",
                "Shallow deletion",
                "Diploid",
                "Gain",
                "Amplification"
            )
        ),

        CDKN2B_CNA = factor(
            as.character(
                CDKN2B_GISTIC
            ),
            levels = c(
                "-2",
                "-1",
                "0",
                "1",
                "2"
            ),
            labels = c(
                "Deep deletion",
                "Shallow deletion",
                "Diploid",
                "Gain",
                "Amplification"
            )
        ),

        EGFR_CNA = factor(
            as.character(
                EGFR_GISTIC
            ),
            levels = c(
                "-2",
                "-1",
                "0",
                "1",
                "2"
            ),
            labels = c(
                "Deep deletion",
                "Shallow deletion",
                "Diploid",
                "Gain",
                "Amplification"
            )
        ),

        TERT_Promoter_Status = factor(
            as.character(
                TERT.promoter.status
            ),
            levels = c(
                "Mutant",
                "WT"
            )
        ),

        Chr7_10_Status = factor(
            as.character(
                Chr.7.gain.Chr.10.loss
            ),
            levels = c(
                "Gain chr 7 & loss chr 10",
                "No combined CNA"
            )
        ),

        MGMT_Promoter_Status = factor(
            as.character(
                MGMT.promoter.status
            ),
            levels = c(
                "Methylated",
                "Unmethylated"
            )
        ),

        OS_MONTHS = suppressWarnings(
            as.numeric(
                OS_months
            )
        ),

        OS_STATUS = factor(
            case_when(
                as.character(OS_status) %in% c(
                    "0",
                    "0.0",
                    "Alive",
                    "alive"
                ) ~ "Alive",

                as.character(OS_status) %in% c(
                    "1",
                    "1.0",
                    "Deceased",
                    "deceased",
                    "Dead",
                    "dead"
                ) ~ "Deceased",

                TRUE ~ NA_character_
            ),
            levels = c(
                "Alive",
                "Deceased"
            )
        )
    )


# ------------------------------------------------------------------------------
# 12. Create concise mutation variables
#
# Unknown MC3 status is converted to NA.
# Each gene is displayed as one row: Mutant n (%).
# ------------------------------------------------------------------------------

for (gene_name in SELECTED_GENES) {

    source_col <- paste0(
        gene_name,
        "_Mutation"
    )

    table_data[[gene_name]] <- factor(
        dplyr::case_when(
            table_data[[source_col]] == "Mutant" ~ "Mutant",
            table_data[[source_col]] == "WT" ~ "WT",
            TRUE ~ NA_character_
        ),
        levels = c(
            "Mutant",
            "WT"
        )
    )
}


# ------------------------------------------------------------------------------
# 13. Final table variables
# ------------------------------------------------------------------------------

analysis_variables <- c(
    # Clinical characteristics
    "WHO_2021_Classification",
    "Age",
    "Sex",

    # Copy-number alterations
    "CDKN2A_CNA",
    "CDKN2B_CNA",
    "EGFR_CNA",
    "Chr7_10_Status",

    # Molecular markers
    "TERT_Promoter_Status",
    "MGMT_Promoter_Status",

    # Somatic mutations
    SELECTED_GENES,

    # Survival
    "OS_MONTHS",
    "OS_STATUS"
)


# ------------------------------------------------------------------------------
# 13A. Sample sizes for table headers
#
# These N values are calculated independently for each DNMT3A grouping variable
# after excluding samples with missing group assignment.
# ------------------------------------------------------------------------------

get_group_n <- function(data, group_variable) {

    x <- as.character(data[[group_variable]])
    x <- x[!is.na(x)]

    c(
        High = sum(x == "High"),
        Low = sum(x == "Low"),
        Total = length(x)
    )
}

n_mrna <- get_group_n(
    table_data,
    "DNMT3A_mRNA_Group"
)

n_promoter <- get_group_n(
    table_data,
    "DNMT3A_Promoter_Methylation_Group"
)

n_cg21629895 <- get_group_n(
    table_data,
    "cg21629895_Group"
)

cat("\nTable-header sample sizes:\n")
cat(
    "  DNMT3A mRNA: High n = ", n_mrna["High"],
    ", Low n = ", n_mrna["Low"],
    ", Total N = ", n_mrna["Total"], "\n",
    sep = ""
)
cat(
    "  Promoter mean methylation: High n = ", n_promoter["High"],
    ", Low n = ", n_promoter["Low"],
    ", Total N = ", n_promoter["Total"], "\n",
    sep = ""
)
cat(
    "  cg21629895: High n = ", n_cg21629895["High"],
    ", Low n = ", n_cg21629895["Low"],
    ", Total N = ", n_cg21629895["Total"], "\n",
    sep = ""
)


variable_labels <- list(
    WHO_2021_Classification ~ "WHO 2021 classification",
    Age ~ "Age, years",
    Sex ~ "Sex",
    CDKN2A_CNA ~ "CDKN2A copy-number status",
    CDKN2B_CNA ~ "CDKN2B copy-number status",
    EGFR_CNA ~ "EGFR copy-number status",
    TERT_Promoter_Status ~ "TERT promoter status",
    Chr7_10_Status ~ "Chromosome 7 gain / chromosome 10 loss",
    MGMT_Promoter_Status ~ "MGMT promoter status",
    FUBP1 ~ "FUBP1 mutation",
    PTEN ~ "PTEN mutation",
    NOTCH1 ~ "NOTCH1 mutation",
    NF1 ~ "NF1 mutation",
    EGFR ~ "EGFR mutation",
    TTN ~ "TTN mutation",
    CIC ~ "CIC mutation",
    ATRX ~ "ATRX mutation",
    TP53 ~ "TP53 mutation",
    OS_MONTHS ~ "Overall survival, months",
    OS_STATUS ~ "Overall survival status"
)


# ------------------------------------------------------------------------------
# 14. Build one comparison table
# ------------------------------------------------------------------------------

build_comparison_table <- function(
        data,
        group_variable
) {

    data_for_table <- data[
        !is.na(
            data[[group_variable]]
        ),
        ,
        drop = FALSE
    ]

    data_for_table %>%
        select(
            all_of(
                c(
                    group_variable,
                    analysis_variables
                )
            )
        ) %>%

        gtsummary::tbl_summary(
            by = all_of(
                group_variable
            ),

            include = all_of(
                analysis_variables
            ),

            type = list(
                all_of(
                    SELECTED_GENES
                ) ~ "dichotomous"
            ),

            value = list(
                all_of(
                    SELECTED_GENES
                ) ~ "Mutant"
            ),

            statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})",
                all_categorical() ~ "{n} ({p}%)"
            ),

            digits = list(
                all_continuous() ~ 1,
                all_categorical() ~ c(
                    0,
                    1
                )
            ),

            label = variable_labels,

            missing = "no"
        ) %>%

        gtsummary::add_p(
            pvalue_fun = gtsummary::label_style_pvalue(
                digits = 3
            )
        ) %>%

        gtsummary::modify_header(
            label ~ "**Characteristic**",
            all_stat_cols() ~ "**{level}**<br>N = {n}",
            p.value ~ "**P value**"
        )
}


# ------------------------------------------------------------------------------
# 15. Build all three existing median-group comparisons
# ------------------------------------------------------------------------------

table_mrna <- build_comparison_table(
    table_data,
    "DNMT3A_mRNA_Group"
)

table_promoter <- build_comparison_table(
    table_data,
    "DNMT3A_Promoter_Methylation_Group"
)

table_cg21629895 <- build_comparison_table(
    table_data,
    "cg21629895_Group"
)


# ------------------------------------------------------------------------------
# 16. Nature-style formatting helper for the split Word tables
#
# Word outputs:
#   1) Split portrait version:
#        A. DNMT3A mRNA expression
#        B. DNMT3A promoter-region mean methylation
#        C. DNMT3A cg21629895 methylation
#
#   2) Integrated landscape version:
#        One 10-column table containing all three comparisons side by side.
#
# Excel output:
#   - Table1_Integrated
#   - Table1_mRNA
#   - Table1_Promoter
#   - Table1_cg21629895
#   - Master_with_mutations
#   - Mutation_QC
# ------------------------------------------------------------------------------

format_nature_flextable <- function(
        gtsummary_table,
        high_n,
        low_n
) {

    ft <- gtsummary::as_flex_table(
        gtsummary_table
    )

    # Clean header labels.
    ft <- flextable::compose(
        ft,
        j = 2,
        value = flextable::as_paragraph(
            flextable::as_chunk(
                paste0(
                    "High\n(n = ",
                    high_n,
                    ")"
                )
            )
        ),
        part = "header"
    )

    ft <- flextable::compose(
        ft,
        j = 3,
        value = flextable::as_paragraph(
            flextable::as_chunk(
                paste0(
                    "Low\n(n = ",
                    low_n,
                    ")"
                )
            )
        ),
        part = "header"
    )

    ft <- flextable::compose(
        ft,
        j = 4,
        value = flextable::as_paragraph(
            flextable::as_chunk(
                "P value"
            )
        ),
        part = "header"
    )

    # Minimal black-and-white journal styling.
    ft <- flextable::border_remove(
        ft
    )

    thin_black <- officer::fp_border(
        color = "#000000",
        width = 0.6
    )

    ft <- flextable::hline_top(
        ft,
        border = thin_black,
        part = "header"
    )

    ft <- flextable::hline_bottom(
        ft,
        border = thin_black,
        part = "header"
    )

    ft <- flextable::hline_bottom(
        ft,
        border = thin_black,
        part = "body"
    )

    ft <- ft %>%
        flextable::font(
            fontname = "Arial",
            part = "all"
        ) %>%
        flextable::fontsize(
            size = 8.5,
            part = "body"
        ) %>%
        flextable::fontsize(
            size = 8.5,
            part = "header"
        ) %>%
        flextable::bold(
            bold = TRUE,
            part = "header"
        ) %>%
        flextable::align(
            j = 1,
            align = "left",
            part = "all"
        ) %>%
        flextable::align(
            j = 2:4,
            align = "center",
            part = "all"
        ) %>%
        flextable::valign(
            valign = "center",
            part = "all"
        ) %>%
        flextable::padding(
            padding.top = 1.5,
            padding.bottom = 1.5,
            padding.left = 2.0,
            padding.right = 2.0,
            part = "all"
        ) %>%
        flextable::set_table_properties(
            layout = "fixed",
            width = 1
        )

    ft <- flextable::width(
        ft,
        j = 1,
        width = 3.45
    )

    ft <- flextable::width(
        ft,
        j = 2,
        width = 1.15
    )

    ft <- flextable::width(
        ft,
        j = 3,
        width = 1.15
    )

    ft <- flextable::width(
        ft,
        j = 4,
        width = 0.85
    )

    ft
}


# ------------------------------------------------------------------------------
# 17. Build the three split Nature-style Word subtables
# ------------------------------------------------------------------------------

format_sectioned_split_flextable <- function(
        df,
        high_n,
        low_n
) {

    section_labels_local <- c(
        "Clinical characteristics",
        "Copy-number alterations",
        "Molecular markers",
        "Somatic mutations",
        "Survival outcomes"
    )

    ft <- flextable::flextable(df)

    ft <- flextable::set_header_labels(
        ft,
        Characteristic = "Characteristic",
        High = paste0("High\n(n = ", high_n, ")"),
        Low = paste0("Low\n(n = ", low_n, ")"),
        `P value` = "P value"
    )

    section_rows <- which(
        df$Characteristic %in% section_labels_local
    )

    ft <- flextable::border_remove(ft)

    rule <- officer::fp_border(
        color = "#000000",
        width = 0.6
    )

    ft <- flextable::hline_top(
        ft,
        border = rule,
        part = "header"
    )

    ft <- flextable::hline_bottom(
        ft,
        border = rule,
        part = "header"
    )

    ft <- flextable::hline_bottom(
        ft,
        border = rule,
        part = "body"
    )

    ft <- ft %>%
        flextable::font(
            fontname = "Arial",
            part = "all"
        ) %>%
        flextable::fontsize(
            size = 7.0,
            part = "body"
        ) %>%
        flextable::fontsize(
            size = 7.2,
            part = "header"
        ) %>%
        flextable::bold(
            bold = TRUE,
            part = "header"
        ) %>%
        flextable::align(
            j = 1,
            align = "left",
            part = "all"
        ) %>%
        flextable::align(
            j = 2:4,
            align = "center",
            part = "all"
        ) %>%
        flextable::valign(
            valign = "center",
            part = "all"
        ) %>%
        flextable::padding(
            padding.top = 1.3,
            padding.bottom = 1.3,
            padding.left = 1.2,
            padding.right = 1.2,
            part = "all"
        ) %>%
        flextable::set_table_properties(
            layout = "fixed",
            width = 1
        )

    # Add subtle separators before each main variable row.
    # This improves within-section readability while avoiding a full grid.
    parent_labels <- c(
        "WHO 2021 classification",
        "Age, years",
        "Sex",
        "CDKN2A copy-number status",
        "CDKN2B copy-number status",
        "EGFR copy-number status",
        "Chromosome 7 gain / chromosome 10 loss",
        "TERT promoter status",
        "MGMT promoter status",
        "FUBP1 mutation",
        "PTEN mutation",
        "NOTCH1 mutation",
        "NF1 mutation",
        "EGFR mutation",
        "TTN mutation",
        "CIC mutation",
        "ATRX mutation",
        "TP53 mutation",
        "Overall survival, months",
        "Overall survival status"
    )

    parent_rows <- which(
        df$Characteristic %in% parent_labels
    )

    subtle_rule <- officer::fp_border(
        color = "#BFBFBF",
        width = 0.35
    )

    if (length(parent_rows) > 0) {
        for (parent_i in parent_rows) {
            line_i <- max(1, parent_i - 1)
            ft <- flextable::hline(
                ft,
                i = line_i,
                border = subtle_rule,
                part = "body"
            )
        }
    }

    if (length(section_rows) > 0) {

        # merge_at() cannot merge several non-consecutive rows in one call.
        # Merge each section divider independently across all four columns.
        for (section_i in section_rows) {
            ft <- flextable::merge_at(
                ft,
                i = section_i,
                j = 1:4,
                part = "body"
            )
        }

        ft <- flextable::bold(
            ft,
            i = section_rows,
            bold = TRUE,
            part = "body"
        )

        ft <- flextable::bg(
            ft,
            i = section_rows,
            bg = "#F2F2F2",
            part = "body"
        )

        ft <- flextable::align(
            ft,
            i = section_rows,
            j = 1,
            align = "left",
            part = "body"
        )

        ft <- flextable::padding(
            ft,
            i = section_rows,
            padding.top = 2.2,
            padding.bottom = 2.2,
            padding.left = 1.2,
            padding.right = 1.2,
            part = "body"
        )
    }

    ft <- flextable::width(
        ft,
        j = 1,
        width = 3.25
    )

    ft <- flextable::width(
        ft,
        j = 2:3,
        width = 1.35
    )

    ft <- flextable::width(
        ft,
        j = 4,
        width = 0.90
    )

    ft
}

# ------------------------------------------------------------------------------
# 18. Build clean data frames from the three gtsummary objects
#
# These are also used to create the integrated Excel and integrated Word table.
# ------------------------------------------------------------------------------

gtsummary_to_clean_df <- function(tbl) {

    out <- gtsummary::as_tibble(
        tbl,
        col_labels = FALSE
    )

    # Retain the first four summary columns only.
    out <- out[
        ,
        seq_len(
            min(
                4,
                ncol(out)
            )
        ),
        drop = FALSE
    ]

    names(out) <- c(
        "Characteristic",
        "High",
        "Low",
        "P value"
    )[seq_len(ncol(out))]

    as.data.frame(
        out,
        check.names = FALSE
    )
}


df_mrna <- gtsummary_to_clean_df(
    table_mrna
)

df_promoter <- gtsummary_to_clean_df(
    table_promoter
)

df_cg21629895 <- gtsummary_to_clean_df(
    table_cg21629895
)


# ------------------------------------------------------------------------------
# Global section labels used by both split and integrated flextable formatters
# ------------------------------------------------------------------------------

SECTION_LABELS <- c(
    "Clinical characteristics",
    "Copy-number alterations",
    "Molecular markers",
    "Somatic mutations",
    "Survival outcomes"
)

# ------------------------------------------------------------------------------
# 18A. Clean labels and add manuscript-style section rows
#
# Section rows replace the visually busy "__Variable__" appearance and organize
# the table into clinically meaningful blocks.
# ------------------------------------------------------------------------------

clean_characteristic_labels <- function(df) {

    df$Characteristic <- gsub(
        "^_+|_+$",
        "",
        df$Characteristic
    )

    df$Characteristic <- gsub(
        "\\*\\*",
        "",
        df$Characteristic
    )

    df
}

add_section_rows <- function(df) {

    section_map <- list(
        "Clinical characteristics" = "WHO 2021 classification",
        "Copy-number alterations" = "CDKN2A copy-number status",
        "Molecular markers" = "TERT promoter status",
        "Somatic mutations" = "FUBP1 mutation",
        "Survival outcomes" = "Overall survival, months"
    )

    out <- df[0, , drop = FALSE]

    for (i in seq_len(nrow(df))) {

        current_label <- trimws(df$Characteristic[i])

        for (section_name in names(section_map)) {

            if (identical(current_label, section_map[[section_name]])) {

                section_row <- as.list(
                    rep("", ncol(df))
                )

                names(section_row) <- names(df)
                section_row$Characteristic <- section_name

                out <- rbind(
                    out,
                    as.data.frame(
                        section_row,
                        check.names = FALSE,
                        stringsAsFactors = FALSE
                    )
                )
            }
        }

        out <- rbind(
            out,
            df[i, , drop = FALSE]
        )
    }

    rownames(out) <- NULL
    out
}

df_mrna <- clean_characteristic_labels(df_mrna)
df_promoter <- clean_characteristic_labels(df_promoter)
df_cg21629895 <- clean_characteristic_labels(df_cg21629895)

df_mrna <- add_section_rows(df_mrna)
df_promoter <- add_section_rows(df_promoter)
df_cg21629895 <- add_section_rows(df_cg21629895)


# ------------------------------------------------------------------------------
# 18B. Build the three split Nature-style Word subtables
#
# IMPORTANT:
# These flextables must be created only AFTER df_mrna, df_promoter, and
# df_cg21629895 have been generated, cleaned, and section rows have been added.
# ------------------------------------------------------------------------------

ft_mrna <- format_sectioned_split_flextable(
    df_mrna,
    high_n = n_mrna["High"],
    low_n = n_mrna["Low"]
)

ft_promoter <- format_sectioned_split_flextable(
    df_promoter,
    high_n = n_promoter["High"],
    low_n = n_promoter["Low"]
)

ft_cg21629895 <- format_sectioned_split_flextable(
    df_cg21629895,
    high_n = n_cg21629895["High"],
    low_n = n_cg21629895["Low"]
)



# ------------------------------------------------------------------------------
# 19. Build the integrated 10-column table
#
# The three gtsummary tables are generated from the same variable list, so their
# row structure should match. The code verifies this before combining them.
# ------------------------------------------------------------------------------

if (
    nrow(df_mrna) != nrow(df_promoter) ||
    nrow(df_mrna) != nrow(df_cg21629895)
) {
    stop(
        paste0(
            "The three summary tables do not contain the same number of rows. ",
            "Integrated output cannot be created safely."
        )
    )
}

integrated_df <- data.frame(
    Characteristic = df_mrna$Characteristic,

    `mRNA High` = df_mrna$High,
    `mRNA Low` = df_mrna$Low,
    `mRNA P value` = df_mrna$`P value`,

    `Promoter High` = df_promoter$High,
    `Promoter Low` = df_promoter$Low,
    `Promoter P value` = df_promoter$`P value`,

    `cg21629895 High` = df_cg21629895$High,
    `cg21629895 Low` = df_cg21629895$Low,
    `cg21629895 P value` = df_cg21629895$`P value`,

    check.names = FALSE,
    stringsAsFactors = FALSE
)

# IMPORTANT:
# flextable requires unique internal column names. Keep unique technical names
# here and replace only the displayed header labels after flextable creation.
names(integrated_df) <- c(
    "Characteristic",
    "mRNA_High",
    "mRNA_Low",
    "mRNA_P",
    "Promoter_High",
    "Promoter_Low",
    "Promoter_P",
    "cg21629895_High",
    "cg21629895_Low",
    "cg21629895_P"
)


# ------------------------------------------------------------------------------
# 20. Common table footnote
# ------------------------------------------------------------------------------

footnote_text <- paste0(
    "Data are presented as n (%) for categorical variables and median (Q1, Q3) ",
    "for continuous variables. DNMT3A mRNA, promoter-region mean methylation, ",
    "and cg21629895 High/Low groups were taken directly from the existing ",
    "master-file median-group variables and were not recalculated. Somatic ",
    "percentages were calculated among samples with available data for each ",
    "variable. Somatic mutation rows represent qualifying non-silent somatic ",
    "alterations. Mutation percentages were calculated among samples with ",
    "available MC3 mutation data; samples absent from the MC3 mutation file were ",
    "treated as missing rather than wild type. ",
    "GISTIC copy-number categories are displayed as deep deletion (-2), shallow ",
    "deletion (-1), diploid (0), gain (1), and amplification (2). ",
    "IDH1 was intentionally excluded from the added mutation rows. Statistical ",
    "tests: Pearson's chi-squared test, Fisher's exact test, or Wilcoxon rank-sum ",
    "test, as appropriate."
)


# ------------------------------------------------------------------------------
# 21. WORD OUTPUT 1: Split portrait Nature-style table
# ------------------------------------------------------------------------------

WORD_SPLIT_FILE <- file.path(
    OUTPUT_DIR,
    "Table_1_DNMT3A_NatureStyle_SPLIT_SECTIONED_CLEAN_FINAL_v10_ROW_DIVIDERS_FIXED.docx"
)

doc_split <- officer::read_docx()

doc_split <- officer::body_add_par(
    doc_split,
    "Table 1 | Clinical and molecular characteristics according to DNMT3A expression and methylation status",
    style = "Normal"
)

doc_split <- officer::body_add_par(
    doc_split,
    "",
    style = "Normal"
)

doc_split <- officer::body_add_par(
    doc_split,
    "A. DNMT3A mRNA expression",
    style = "Normal"
)

doc_split <- flextable::body_add_flextable(
    doc_split,
    ft_mrna,
    align = "left"
)

doc_split <- officer::body_add_par(
    doc_split,
    "",
    style = "Normal"
)

doc_split <- officer::body_add_par(
    doc_split,
    "B. DNMT3A promoter-region mean methylation",
    style = "Normal"
)

doc_split <- flextable::body_add_flextable(
    doc_split,
    ft_promoter,
    align = "left"
)

doc_split <- officer::body_add_par(
    doc_split,
    "",
    style = "Normal"
)

doc_split <- officer::body_add_par(
    doc_split,
    "C. DNMT3A cg21629895 site-specific methylation",
    style = "Normal"
)

doc_split <- flextable::body_add_flextable(
    doc_split,
    ft_cg21629895,
    align = "left"
)

doc_split <- officer::body_add_par(
    doc_split,
    "",
    style = "Normal"
)

doc_split <- officer::body_add_par(
    doc_split,
    footnote_text,
    style = "Normal"
)

assert_output_ready(WORD_SPLIT_FILE)

print(
    doc_split,
    target = WORD_SPLIT_FILE
)


# ------------------------------------------------------------------------------
# 22. WORD OUTPUT 2: Integrated landscape Nature-style table
#
# This version keeps all three comparisons together in one table. It is saved
# separately because 10 columns require landscape orientation and smaller text.
# ------------------------------------------------------------------------------

ft_integrated <- flextable::flextable(
    integrated_df
)

# Display publication-facing subgroup headers while preserving unique internal
# column keys required by flextable.
ft_integrated <- flextable::set_header_labels(
    ft_integrated,
    Characteristic = "Characteristic",
    mRNA_High = paste0("High\n(n = ", n_mrna["High"], ")"),
    mRNA_Low = paste0("Low\n(n = ", n_mrna["Low"], ")"),
    mRNA_P = "P value",
    Promoter_High = paste0("High\n(n = ", n_promoter["High"], ")"),
    Promoter_Low = paste0("Low\n(n = ", n_promoter["Low"], ")"),
    Promoter_P = "P value",
    cg21629895_High = paste0("High\n(n = ", n_cg21629895["High"], ")"),
    cg21629895_Low = paste0("Low\n(n = ", n_cg21629895["Low"], ")"),
    cg21629895_P = "P value"
)

# First header row: detailed column names.
ft_integrated <- flextable::border_remove(
    ft_integrated
)

thin_black <- officer::fp_border(
    color = "#000000",
    width = 0.6
)

ft_integrated <- flextable::hline_top(
    ft_integrated,
    border = thin_black,
    part = "header"
)

ft_integrated <- flextable::hline_bottom(
    ft_integrated,
    border = thin_black,
    part = "header"
)

ft_integrated <- flextable::hline_bottom(
    ft_integrated,
    border = thin_black,
    part = "body"
)

# Add grouped spanner header.
ft_integrated <- flextable::add_header_row(
    ft_integrated,
    values = c(
        "",
        paste0(
            "DNMT3A mRNA expression (N = ",
            n_mrna["Total"],
            ")"
        ),
        paste0(
            "DNMT3A promoter-region mean methylation (N = ",
            n_promoter["Total"],
            ")"
        ),
        paste0(
            "DNMT3A cg21629895 site-specific methylation (N = ",
            n_cg21629895["Total"],
            ")"
        )
    ),
    colwidths = c(
        1,
        3,
        3,
        3
    )
)

ft_integrated <- ft_integrated %>%
    flextable::font(
        fontname = "Arial",
        part = "all"
    ) %>%
    flextable::fontsize(
        size = 6.4,
        part = "body"
    ) %>%
    flextable::fontsize(
        size = 6.6,
        part = "header"
    ) %>%
    flextable::bold(
        bold = TRUE,
        part = "header"
    ) %>%
    flextable::align(
        j = 1,
        align = "left",
        part = "all"
    ) %>%
    flextable::align(
        j = 2:10,
        align = "center",
        part = "all"
    ) %>%
    flextable::valign(
        valign = "center",
        part = "all"
    ) %>%
    flextable::padding(
        padding.top = 1.0,
        padding.bottom = 1.0,
        padding.left = 1.0,
        padding.right = 1.0,
        part = "all"
    ) %>%
    flextable::set_table_properties(
        layout = "fixed",
        width = 1
    )

# Landscape page widths in inches.
ft_integrated <- flextable::width(
    ft_integrated,
    j = 1,
    width = 2.70
)

for (j in c(2, 3, 5, 6, 8, 9)) {
    ft_integrated <- flextable::width(
        ft_integrated,
        j = j,
        width = 0.88
    )
}

for (j in c(4, 7, 10)) {
    ft_integrated <- flextable::width(
        ft_integrated,
        j = j,
        width = 0.68
    )
}


# Add subtle separators before each main variable row.
integrated_parent_labels <- c(
    "WHO 2021 classification",
    "Age, years",
    "Sex",
    "CDKN2A copy-number status",
    "CDKN2B copy-number status",
    "EGFR copy-number status",
    "Chromosome 7 gain / chromosome 10 loss",
    "TERT promoter status",
    "MGMT promoter status",
    "FUBP1 mutation",
    "PTEN mutation",
    "NOTCH1 mutation",
    "NF1 mutation",
    "EGFR mutation",
    "TTN mutation",
    "CIC mutation",
    "ATRX mutation",
    "TP53 mutation",
    "Overall survival, months",
    "Overall survival status"
)

integrated_parent_rows <- which(
    integrated_df$Characteristic %in% integrated_parent_labels
)

integrated_subtle_rule <- officer::fp_border(
    color = "#BFBFBF",
    width = 0.35
)

if (length(integrated_parent_rows) > 0) {
    for (parent_i in integrated_parent_rows) {
        line_i <- max(1, parent_i - 1)
        ft_integrated <- flextable::hline(
            ft_integrated,
            i = line_i,
            border = integrated_subtle_rule,
            part = "body"
        )
    }
}

# Emphasize clean section dividers across the full integrated table.
integrated_section_rows <- which(
    integrated_df$Characteristic %in% SECTION_LABELS
)

if (length(integrated_section_rows) > 0) {

    # Section rows are non-consecutive; merge each divider independently.
    for (section_i in integrated_section_rows) {
        ft_integrated <- flextable::merge_at(
            ft_integrated,
            i = section_i,
            j = 1:10,
            part = "body"
        )
    }

    ft_integrated <- flextable::bold(
        ft_integrated,
        i = integrated_section_rows,
        bold = TRUE,
        part = "body"
    )

    ft_integrated <- flextable::bg(
        ft_integrated,
        i = integrated_section_rows,
        bg = "#F2F2F2",
        part = "body"
    )

    ft_integrated <- flextable::align(
        ft_integrated,
        i = integrated_section_rows,
        j = 1,
        align = "left",
        part = "body"
    )

    ft_integrated <- flextable::padding(
        ft_integrated,
        i = integrated_section_rows,
        padding.top = 2.0,
        padding.bottom = 2.0,
        padding.left = 1.0,
        padding.right = 1.0,
        part = "body"
    )
}

WORD_INTEGRATED_FILE <- file.path(
    OUTPUT_DIR,
    "Table_1_DNMT3A_NatureStyle_INTEGRATED_SECTIONED_CLEAN_FINAL_v10_ROW_DIVIDERS_FIXED.docx"
)

doc_integrated <- officer::read_docx()

doc_integrated <- officer::body_add_par(
    doc_integrated,
    "Table 1 | Clinical and molecular characteristics according to DNMT3A expression and methylation status",
    style = "Normal"
)

doc_integrated <- flextable::body_add_flextable(
    doc_integrated,
    ft_integrated,
    align = "center"
)

doc_integrated <- officer::body_add_par(
    doc_integrated,
    "",
    style = "Normal"
)

doc_integrated <- officer::body_add_par(
    doc_integrated,
    footnote_text,
    style = "Normal"
)

# Landscape orientation for the integrated table.
doc_integrated <- officer::body_end_section_landscape(
    doc_integrated
)

assert_output_ready(WORD_INTEGRATED_FILE)

print(
    doc_integrated,
    target = WORD_INTEGRATED_FILE
)


# ------------------------------------------------------------------------------
# 23. EXCEL OUTPUT
#
# Includes BOTH:
#   - one integrated summary sheet
#   - three separate comparison sheets
# ------------------------------------------------------------------------------

EXCEL_FILE <- file.path(
    OUTPUT_DIR,
    "Table_1_DNMT3A_NatureStyle_SECTIONED_CLEAN_FINAL_v10_ROW_DIVIDERS_FIXED.xlsx"
)

wb <- openxlsx::createWorkbook()


# ------------------------------------------------------------------------------
# 23A. Excel style definitions
# ------------------------------------------------------------------------------

style_title <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 11,
    textDecoration = "bold",
    halign = "left",
    valign = "center"
)

style_spanner <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = c(
        "top",
        "bottom"
    ),
    borderStyle = "thin",
    borderColour = "#000000"
)

style_header <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = c(
        "bottom"
    ),
    borderStyle = "thin",
    borderColour = "#000000"
)

style_body_left <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    halign = "left",
    valign = "center"
)

style_body_center <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    halign = "center",
    valign = "center"
)

style_note <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 8,
    fontColour = "#444444",
    wrapText = TRUE,
    valign = "top"
)


# ------------------------------------------------------------------------------
# 23B. Helper for split Excel summary sheets
# ------------------------------------------------------------------------------

write_nature_sheet <- function(
        wb,
        sheet_name,
        title_text,
        df,
        high_n,
        low_n,
        total_n
) {

    df <- as.data.frame(
        df,
        check.names = FALSE
    )

    names(df) <- c(
        "Characteristic",
        paste0("High (n = ", high_n, ")"),
        paste0("Low (n = ", low_n, ")"),
        "P value"
    )

    title_text <- paste0(
        title_text,
        " (N = ",
        total_n,
        ")"
    )

    openxlsx::addWorksheet(
        wb,
        sheet_name,
        gridLines = FALSE
    )

    openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x = title_text,
        startRow = 1,
        startCol = 1
    )

    openxlsx::mergeCells(
        wb,
        sheet = sheet_name,
        cols = 1:4,
        rows = 1
    )

    openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = style_title,
        rows = 1,
        cols = 1:4,
        gridExpand = TRUE
    )

    openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x = df,
        startRow = 3,
        startCol = 1,
        headerStyle = style_header,
        borders = "none"
    )

    last_row <- 3 + nrow(df)

    openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = style_body_left,
        rows = 4:last_row,
        cols = 1,
        gridExpand = TRUE
    )

    openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = style_body_center,
        rows = 4:last_row,
        cols = 2:4,
        gridExpand = TRUE
    )

    note_row <- last_row + 2

    openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x = footnote_text,
        startRow = note_row,
        startCol = 1
    )

    openxlsx::mergeCells(
        wb,
        sheet = sheet_name,
        cols = 1:4,
        rows = note_row
    )

    openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = style_note,
        rows = note_row,
        cols = 1:4,
        gridExpand = TRUE
    )

    openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = 1,
        widths = 46
    )

    openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = 2:3,
        widths = 18
    )

    openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = 4,
        widths = 12
    )

    openxlsx::freezePane(
        wb,
        sheet = sheet_name,
        firstActiveRow = 4
    )
}


# ------------------------------------------------------------------------------
# 23C. Integrated Excel summary sheet
# ------------------------------------------------------------------------------

openxlsx::addWorksheet(
    wb,
    "Table1_Integrated",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    "Table 1. Clinical and molecular characteristics according to DNMT3A expression and methylation status",
    startRow = 1,
    startCol = 1
)

openxlsx::mergeCells(
    wb,
    "Table1_Integrated",
    cols = 1:10,
    rows = 1
)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_title,
    rows = 1,
    cols = 1:10,
    gridExpand = TRUE
)

# Grouped headers.
openxlsx::writeData(
    wb,
    "Table1_Integrated",
    x = "Characteristic",
    startRow = 3,
    startCol = 1
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    x = paste0(
        "DNMT3A mRNA expression (N = ",
        n_mrna["Total"],
        ")"
    ),
    startRow = 3,
    startCol = 2
)

openxlsx::mergeCells(
    wb,
    "Table1_Integrated",
    cols = 2:4,
    rows = 3
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    x = paste0(
        "DNMT3A promoter-region mean methylation (N = ",
        n_promoter["Total"],
        ")"
    ),
    startRow = 3,
    startCol = 5
)

openxlsx::mergeCells(
    wb,
    "Table1_Integrated",
    cols = 5:7,
    rows = 3
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    x = paste0(
        "DNMT3A cg21629895 site-specific methylation (N = ",
        n_cg21629895["Total"],
        ")"
    ),
    startRow = 3,
    startCol = 8
)

openxlsx::mergeCells(
    wb,
    "Table1_Integrated",
    cols = 8:10,
    rows = 3
)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_spanner,
    rows = 3,
    cols = 1:10,
    gridExpand = TRUE
)

integrated_subheaders <- data.frame(
    Characteristic = "",
    `mRNA High` = paste0("High\n(n = ", n_mrna["High"], ")"),
    `mRNA Low` = paste0("Low\n(n = ", n_mrna["Low"], ")"),
    `mRNA P value` = "P value",
    `Promoter High` = paste0("High\n(n = ", n_promoter["High"], ")"),
    `Promoter Low` = paste0("Low\n(n = ", n_promoter["Low"], ")"),
    `Promoter P value` = "P value",
    `cg21629895 High` = paste0("High\n(n = ", n_cg21629895["High"], ")"),
    `cg21629895 Low` = paste0("Low\n(n = ", n_cg21629895["Low"], ")"),
    `cg21629895 P value` = "P value",
    check.names = FALSE
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    integrated_subheaders,
    startRow = 4,
    startCol = 1,
    colNames = FALSE
)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_header,
    rows = 4,
    cols = 1:10,
    gridExpand = TRUE
)

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    integrated_df,
    startRow = 5,
    startCol = 1,
    colNames = FALSE
)

integrated_last_row <- 4 + nrow(integrated_df)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_body_left,
    rows = 5:integrated_last_row,
    cols = 1,
    gridExpand = TRUE
)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_body_center,
    rows = 5:integrated_last_row,
    cols = 2:10,
    gridExpand = TRUE
)

integrated_note_row <- integrated_last_row + 2

openxlsx::writeData(
    wb,
    "Table1_Integrated",
    footnote_text,
    startRow = integrated_note_row,
    startCol = 1
)

openxlsx::mergeCells(
    wb,
    "Table1_Integrated",
    cols = 1:10,
    rows = integrated_note_row
)

openxlsx::addStyle(
    wb,
    "Table1_Integrated",
    style_note,
    rows = integrated_note_row,
    cols = 1:10,
    gridExpand = TRUE
)

openxlsx::setColWidths(
    wb,
    "Table1_Integrated",
    cols = 1,
    widths = 44
)

openxlsx::setColWidths(
    wb,
    "Table1_Integrated",
    cols = c(
        2,
        3,
        5,
        6,
        8,
        9
    ),
    widths = 15
)

openxlsx::setColWidths(
    wb,
    "Table1_Integrated",
    cols = c(
        4,
        7,
        10
    ),
    widths = 11
)

openxlsx::freezePane(
    wb,
    "Table1_Integrated",
    firstActiveRow = 5,
    firstActiveCol = 2
)


# ------------------------------------------------------------------------------
# 23D. Separate Excel comparison sheets
# ------------------------------------------------------------------------------

write_nature_sheet(
    wb,
    "Table1_mRNA",
    "Table 1A. DNMT3A mRNA expression",
    df_mrna,
    high_n = n_mrna["High"],
    low_n = n_mrna["Low"],
    total_n = n_mrna["Total"]
)

write_nature_sheet(
    wb,
    "Table1_Promoter",
    "Table 1B. DNMT3A promoter-region mean methylation",
    df_promoter,
    high_n = n_promoter["High"],
    low_n = n_promoter["Low"],
    total_n = n_promoter["Total"]
)

write_nature_sheet(
    wb,
    "Table1_cg21629895",
    "Table 1C. DNMT3A cg21629895 site-specific methylation",
    df_cg21629895,
    high_n = n_cg21629895["High"],
    low_n = n_cg21629895["Low"],
    total_n = n_cg21629895["Total"]
)


# ------------------------------------------------------------------------------
# 23E. Merged master and mutation QC sheets
# ------------------------------------------------------------------------------

openxlsx::addWorksheet(
    wb,
    "Master_with_mutations",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Master_with_mutations",
    as.data.frame(
        master
    ),
    withFilter = TRUE
)

openxlsx::freezePane(
    wb,
    "Master_with_mutations",
    firstRow = TRUE
)

openxlsx::addWorksheet(
    wb,
    "Mutation_QC",
    gridLines = FALSE
)

openxlsx::writeData(
    wb,
    "Mutation_QC",
    mutation_counts
)

if (length(unmatched_samples) > 0) {

    unmatched_df <- data.frame(
        Unmatched_Sample_Core = sort(
            unmatched_samples
        )
    )

    openxlsx::writeData(
        wb,
        "Mutation_QC",
        unmatched_df,
        startRow = nrow(
            mutation_counts
        ) + 4
    )
}

assert_output_ready(EXCEL_FILE)

openxlsx::saveWorkbook(
    wb,
    EXCEL_FILE,
    overwrite = TRUE
)


# ------------------------------------------------------------------------------
# 24. Save analysis and QC CSVs
# ------------------------------------------------------------------------------

ANALYSIS_FILE <- file.path(
    OUTPUT_DIR,
    "Table_1_Analysis_Data_FINAL.csv"
)

analysis_output <- table_data %>%
    select(
        Patient_ID,
        Sample_Core,
        DNMT3A_mRNA_Group,
        DNMT3A_Promoter_Methylation_Group,
        cg21629895_Group,
        all_of(
            analysis_variables
        )
    )

write.csv(
    analysis_output,
    ANALYSIS_FILE,
    row.names = FALSE,
    na = ""
)


QC_FILE <- file.path(
    OUTPUT_DIR,
    "Mutation_Merge_QC.csv"
)

qc_summary <- data.table::rbindlist(
    list(
        data.table::data.table(
            Metric = c(
                "Master samples",
                "MC3 unique samples",
                "Matched samples",
                "Unmatched master samples"
            ),
            Value = as.character(
                c(
                    length(master_samples),
                    length(all_mc3_samples),
                    length(matched_samples),
                    length(unmatched_samples)
                )
            )
        ),

        data.table::data.table(
            Metric = paste0(
                mutation_counts$Gene,
                " mutant samples"
            ),
            Value = as.character(
                mutation_counts$Mutant
            )
        ),

        data.table::data.table(
            Metric = "Unmatched Sample_Core IDs",
            Value = paste(
                sort(
                    unmatched_samples
                ),
                collapse = "; "
            )
        )
    ),
    fill = TRUE
)

data.table::fwrite(
    qc_summary,
    QC_FILE
)


# ------------------------------------------------------------------------------
# 25. Final console report
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL DNMT3A NATURE-STYLE TABLE OUTPUT COMPLETE\n")
cat("============================================================\n")

cat("\nWord outputs:\n")

cat(
    "  Split portrait: ",
    WORD_SPLIT_FILE,
    "\n",
    sep = ""
)

cat(
    "  Integrated landscape: ",
    WORD_INTEGRATED_FILE,
    "\n",
    sep = ""
)

cat("\nExcel output:\n")

cat(
    "  ",
    EXCEL_FILE,
    "\n",
    sep = ""
)

cat("\nExcel summary sheets:\n")
cat("  Table1_Integrated\n")
cat("  Table1_mRNA\n")
cat("  Table1_Promoter\n")
cat("  Table1_cg21629895\n")
cat("  Master_with_mutations\n")
cat("  Mutation_QC\n")

cat("\nOther outputs:\n")

cat(
    "  Merged master CSV: ",
    MERGED_MASTER_FILE,
    "\n",
    sep = ""
)

cat(
    "  Analysis CSV: ",
    ANALYSIS_FILE,
    "\n",
    sep = ""
)

cat(
    "  Mutation QC CSV: ",
    QC_FILE,
    "\n",
    sep = ""
)

cat("\n============================================================\n")
