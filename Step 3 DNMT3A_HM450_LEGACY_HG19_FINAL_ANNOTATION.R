# =============================================================================
# DNMT3A HumanMethylation450 Probe Annotation
# UCSC Xena TCGA legacy hg19 probeMap + Official Illumina 450K annotation
#
# Purpose
# -------
# Extract all DNMT3A CpG probes present in the TCGA-LGG master file and create
# a manuscript-ready annotation table with:
#
#   Name
#   chr
#   pos
#   UCSC_RefGene_Name
#   UCSC_RefGene_Group
#   Relation_to_Island
#
# Coordinate source:
#   UCSC Xena TCGA legacy HumanMethylation450 hg19 probe map
#   illuminaMethyl450_hg19_GPL16304_TCGAlegacy
#
# Biological annotation source:
#   Illumina HumanMethylation450 v1.2 annotation through
#   IlluminaHumanMethylation450kanno.ilmn12.hg19
#
# IMPORTANT:
#   - chr and pos in the final output come ONLY from the UCSC Xena TCGA legacy
#     hg19 probeMap used with the current legacy HumanMethylation450 dataset.
#   - The analysis remains on hg19; no hg38 liftover is performed.
#   - UCSC_RefGene_Name, UCSC_RefGene_Group, and Relation_to_Island are taken
#     from the official Illumina 450K hg19 annotation.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install/load required packages
# -----------------------------------------------------------------------------

cran_packages <- c(
    "data.table",
    "dplyr",
    "openxlsx"
)

missing_cran <- cran_packages[
    !vapply(
        cran_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing_cran) > 0) {
    install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
    "minfi",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_bioc <- bioc_packages[
    !vapply(
        bioc_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing_bioc) > 0) {
    BiocManager::install(
        missing_bioc,
        ask = FALSE,
        update = FALSE
    )
}

library(data.table)
library(dplyr)
library(openxlsx)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# -----------------------------------------------------------------------------
# 1. Select the TCGA-LGG integrated master file
# -----------------------------------------------------------------------------

cat("\nSelect the TCGA-LGG master CSV file.\n")
master_file <- file.choose()

master <- fread(
    master_file,
    check.names = FALSE
)

cat("\nMaster file loaded.\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")


# -----------------------------------------------------------------------------
# 2. Identify all CpG probe columns in the master file
# -----------------------------------------------------------------------------

probe_columns <- grep(
    "^cg[0-9]+$",
    names(master),
    value = TRUE
)

if (length(probe_columns) == 0) {
    stop(
        "No CpG probe columns matching the pattern cg######## were found."
    )
}

cat("\nCpG probe columns detected:", length(probe_columns), "\n")


# -----------------------------------------------------------------------------
# 3. Obtain the exact UCSC Xena TCGA legacy hg19 HM450 probeMap
# -----------------------------------------------------------------------------

probe_map_url <- paste0(
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/",
    "probeMap%2FilluminaMethyl450_hg19_GPL16304_TCGAlegacy"
)

probe_map_file <- file.path(
    tempdir(),
    "illuminaMethyl450_hg19_GPL16304_TCGAlegacy"
)

cat("\nDownloading the UCSC Xena TCGA legacy hg19 HM450 probe map...\n")
cat(probe_map_url, "\n")

download_ok <- TRUE

tryCatch(
    {
        utils::download.file(
            probe_map_url,
            destfile = probe_map_file,
            mode = "wb",
            quiet = FALSE
        )
    },
    error = function(e) {
        download_ok <<- FALSE
        cat(
            "\nAutomatic download failed:",
            conditionMessage(e),
            "\n"
        )
    }
)

if (
    !download_ok ||
    !file.exists(probe_map_file) ||
    file.info(probe_map_file)$size == 0
) {
    cat(
        "\nPlease select the locally downloaded ",
        "illuminaMethyl450_hg19_GPL16304_TCGAlegacy probeMap file.\n",
        sep = ""
    )

    probe_map_file <- file.choose()
}

cat("\nReading legacy hg19 probe map...\n")

probe_map <- data.table::fread(
    probe_map_file,
    check.names = FALSE
)

cat("\nProbe map dimensions:\n")
cat("Rows:", nrow(probe_map), "\n")
cat("Columns:", ncol(probe_map), "\n")

cat("\nProbe map column names:\n")
print(names(probe_map))


# -----------------------------------------------------------------------------
# 4. Detect key columns in the legacy hg19 probeMap
# -----------------------------------------------------------------------------

find_first_column <- function(
    column_names,
    candidates,
    required = TRUE
) {

    exact_match <- candidates[
        candidates %in% column_names
    ]

    if (length(exact_match) > 0) {
        return(exact_match[1])
    }

    lower_names <- tolower(column_names)
    lower_candidates <- tolower(candidates)

    idx <- match(
        lower_candidates,
        lower_names
    )

    idx <- idx[
        !is.na(idx)
    ]

    if (length(idx) > 0) {
        return(column_names[idx[1]])
    }

    if (required) {
        stop(
            paste0(
                "Could not identify a required probeMap column. ",
                "Candidates checked: ",
                paste(candidates, collapse = ", ")
            )
        )
    }

    return(NA_character_)
}


probe_id_col <- find_first_column(
    names(probe_map),
    c(
        "id",
        "#id",
        "Name",
        "name",
        "probe",
        "Probe_ID",
        "IlmnID"
    )
)

chr_col <- find_first_column(
    names(probe_map),
    c(
        "chrom",
        "chr",
        "Chromosome",
        "chromosome"
    )
)

pos_col <- find_first_column(
    names(probe_map),
    c(
        "chromStart",
        "pos",
        "position",
        "Position",
        "start",
        "Start"
    )
)

cat("\nDetected legacy hg19 probeMap columns:\n")
cat("Probe ID:", probe_id_col, "\n")
cat("Chromosome:", chr_col, "\n")
cat("Position:", pos_col, "\n")


# -----------------------------------------------------------------------------
# 5. Build legacy hg19 coordinate table
# -----------------------------------------------------------------------------

hg19_map <- probe_map %>%
    dplyr::transmute(
        Name = as.character(.data[[probe_id_col]]),
        chr = as.character(.data[[chr_col]]),
        pos_raw = suppressWarnings(
            as.numeric(.data[[pos_col]])
        )
    ) %>%
    dplyr::filter(
        !is.na(Name),
        Name != ""
    ) %>%
    dplyr::distinct(
        Name,
        .keep_all = TRUE
    )


# -----------------------------------------------------------------------------
# 6. Retain legacy hg19 coordinates exactly as provided by the Xena probeMap
# -----------------------------------------------------------------------------

hg19_map <- hg19_map %>%
    dplyr::mutate(
        pos = as.integer(pos_raw)
    ) %>%
    dplyr::select(
        Name,
        chr,
        pos
    )

cat(
    "\nLegacy hg19 coordinates were retained exactly as provided by the ",
    "UCSC Xena TCGA legacy probeMap. No liftover or coordinate offset was applied.\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 7. Retrieve official Illumina HumanMethylation450 annotation
# -----------------------------------------------------------------------------

cat("\nLoading official Illumina HumanMethylation450 annotation...\n")

illumina_annotation <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19
)

illumina_annotation <- as.data.frame(
    illumina_annotation,
    stringsAsFactors = FALSE
)

illumina_annotation$Name <- rownames(
    illumina_annotation
)


# -----------------------------------------------------------------------------
# 8. Confirm required Illumina annotation fields
# -----------------------------------------------------------------------------

required_illumina_columns <- c(
    "Name",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Group",
    "Relation_to_Island"
)

missing_illumina_columns <- setdiff(
    required_illumina_columns,
    names(illumina_annotation)
)

if (length(missing_illumina_columns) > 0) {
    stop(
        paste(
            "Missing required Illumina annotation columns:",
            paste(
                missing_illumina_columns,
                collapse = ", "
            )
        )
    )
}


# -----------------------------------------------------------------------------
# 9. Create the official 450K biological annotation table
#
# hg19 coordinate fields from this package are deliberately excluded.
# -----------------------------------------------------------------------------

illumina_table <- illumina_annotation %>%
    dplyr::transmute(
        Name = as.character(Name),
        UCSC_RefGene_Name = as.character(UCSC_RefGene_Name),
        UCSC_RefGene_Group = as.character(UCSC_RefGene_Group),
        Relation_to_Island = as.character(Relation_to_Island)
    ) %>%
    dplyr::distinct(
        Name,
        .keep_all = TRUE
    )


# -----------------------------------------------------------------------------
# 10. Restrict annotation to probes actually present in the master file
# -----------------------------------------------------------------------------

master_probe_table <- data.frame(
    Name = probe_columns,
    stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 11. Combine legacy hg19 coordinates with Illumina 450K annotation
# -----------------------------------------------------------------------------

final_annotation <- master_probe_table %>%
    dplyr::left_join(
        hg19_map,
        by = "Name"
    ) %>%
    dplyr::left_join(
        illumina_table,
        by = "Name"
    )


# -----------------------------------------------------------------------------
# 12. Add data-completeness QC from the actual master cohort
# -----------------------------------------------------------------------------

probe_qc <- lapply(
    probe_columns,
    function(probe_id) {

        values <- suppressWarnings(
            as.numeric(master[[probe_id]])
        )

        data.frame(
            Name = probe_id,
            N_with_beta = sum(!is.na(values)),
            N_missing = sum(is.na(values)),
            Completeness_percent = round(
                mean(!is.na(values)) * 100,
                2
            ),
            stringsAsFactors = FALSE
        )
    }
)

probe_qc <- dplyr::bind_rows(
    probe_qc
)

final_annotation_qc <- final_annotation %>%
    dplyr::left_join(
        probe_qc,
        by = "Name"
    )


# -----------------------------------------------------------------------------
# 13. Sort probes by genomic location
# -----------------------------------------------------------------------------

final_annotation <- final_annotation %>%
    dplyr::arrange(
        chr,
        pos
    )

final_annotation_qc <- final_annotation_qc %>%
    dplyr::arrange(
        chr,
        pos
    )


# -----------------------------------------------------------------------------
# 14. QC checks
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ANNOTATION QC\n")
cat("============================================================\n")

cat("\nProbes in master:", length(probe_columns), "\n")

cat(
    "Probes matched to legacy hg19 probeMap:",
    sum(!is.na(final_annotation$pos)),
    "/",
    nrow(final_annotation),
    "\n"
)

cat(
    "Probes matched to Illumina RefGene annotation:",
    sum(
        !is.na(final_annotation$UCSC_RefGene_Name) &
        final_annotation$UCSC_RefGene_Name != ""
    ),
    "/",
    nrow(final_annotation),
    "\n"
)

cat(
    "Probes with CpG island annotation:",
    sum(
        !is.na(final_annotation$Relation_to_Island) &
        final_annotation$Relation_to_Island != ""
    ),
    "/",
    nrow(final_annotation),
    "\n"
)


# -----------------------------------------------------------------------------
# 15. Identify unmatched probes
# -----------------------------------------------------------------------------

unmatched_hg19 <- final_annotation %>%
    dplyr::filter(
        is.na(chr) |
        is.na(pos)
    ) %>%
    dplyr::select(
        Name
    )

unmatched_illumina <- final_annotation %>%
    dplyr::filter(
        is.na(UCSC_RefGene_Name) |
        UCSC_RefGene_Name == ""
    ) %>%
    dplyr::select(
        Name
    )

if (nrow(unmatched_hg19) > 0) {
    cat("\nProbes not matched to legacy hg19 probeMap:\n")
    print(unmatched_hg19)
}

if (nrow(unmatched_illumina) > 0) {
    cat("\nProbes without RefGene annotation:\n")
    print(unmatched_illumina)
}


# -----------------------------------------------------------------------------
# 16. Final manuscript table
# -----------------------------------------------------------------------------

paper_table <- final_annotation %>%
    dplyr::select(
        Name,
        chr,
        pos,
        UCSC_RefGene_Name,
        UCSC_RefGene_Group,
        Relation_to_Island
    )


# -----------------------------------------------------------------------------
# 17. Save CSV
# -----------------------------------------------------------------------------

output_dir <- dirname(
    master_file
)

csv_file <- file.path(
    output_dir,
    "DNMT3A_HM450_LEGACY_HG19_PROBE_ANNOTATION_PAPER_TABLE.csv"
)

fwrite(
    paper_table,
    csv_file,
    na = ""
)


# -----------------------------------------------------------------------------
# 18. Save Excel workbook
# -----------------------------------------------------------------------------

xlsx_file <- file.path(
    output_dir,
    "DNMT3A_HM450_LEGACY_HG19_PROBE_ANNOTATION_PAPER_TABLE.xlsx"
)

wb <- createWorkbook()

addWorksheet(
    wb,
    "Paper_Table"
)

addWorksheet(
    wb,
    "QC"
)

addWorksheet(
    wb,
    "Annotation_Source"
)


# -----------------------------------------------------------------------------
# 18A. Paper table
# -----------------------------------------------------------------------------

writeData(
    wb,
    sheet = "Paper_Table",
    x = paper_table,
    startRow = 1
)

header_style <- createStyle(
    fontName = "Arial",
    fontSize = 10,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom"
)

body_style <- createStyle(
    fontName = "Arial",
    fontSize = 9,
    valign = "center"
)

addStyle(
    wb,
    sheet = "Paper_Table",
    style = header_style,
    rows = 1,
    cols = 1:ncol(paper_table),
    gridExpand = TRUE
)

addStyle(
    wb,
    sheet = "Paper_Table",
    style = body_style,
    rows = 2:(nrow(paper_table) + 1),
    cols = 1:ncol(paper_table),
    gridExpand = TRUE
)

setColWidths(
    wb,
    sheet = "Paper_Table",
    cols = 1:6,
    widths = c(
        15,
        10,
        14,
        38,
        48,
        22
    )
)

freezePane(
    wb,
    sheet = "Paper_Table",
    firstRow = TRUE
)

addFilter(
    wb,
    sheet = "Paper_Table",
    row = 1,
    cols = 1:ncol(paper_table)
)


# -----------------------------------------------------------------------------
# 18B. QC table
# -----------------------------------------------------------------------------

writeData(
    wb,
    sheet = "QC",
    x = final_annotation_qc
)

addStyle(
    wb,
    sheet = "QC",
    style = header_style,
    rows = 1,
    cols = 1:ncol(final_annotation_qc),
    gridExpand = TRUE
)

setColWidths(
    wb,
    sheet = "QC",
    cols = 1:ncol(final_annotation_qc),
    widths = "auto"
)

freezePane(
    wb,
    sheet = "QC",
    firstRow = TRUE
)


# -----------------------------------------------------------------------------
# 18C. Annotation source information
# -----------------------------------------------------------------------------

source_table <- data.frame(
    Field = c(
        "Methylation platform",
        "Coordinate assembly",
        "Coordinate source",
        "Probe map URL",
        "Gene annotation source",
        "Coordinate policy",
        "Main output columns"
    ),
    Description = c(
        "Illumina HumanMethylation450 BeadChip",
        "GRCh37 / hg19",
        "UCSC Xena TCGA legacy illuminaMethyl450_hg19_GPL16304_TCGAlegacy",
        probe_map_url,
        paste0(
            "Illumina HumanMethylation450 v1.2 annotation via ",
            "IlluminaHumanMethylation450kanno.ilmn12.hg19"
        ),
        paste0(
            "chr and pos are derived exclusively from the GDC/UCSC Xena ",
            "legacy hg19 probeMap. hg19 coordinates from the Illumina annotation ",
            "package are not included."
        ),
        paste(
            names(paper_table),
            collapse = "; "
        )
    ),
    stringsAsFactors = FALSE
)

writeData(
    wb,
    sheet = "Annotation_Source",
    x = source_table
)

addStyle(
    wb,
    sheet = "Annotation_Source",
    style = header_style,
    rows = 1,
    cols = 1:2,
    gridExpand = TRUE
)

setColWidths(
    wb,
    sheet = "Annotation_Source",
    cols = c(1, 2),
    widths = c(28, 100)
)

wrap_style <- createStyle(
    wrapText = TRUE,
    valign = "top",
    fontName = "Arial",
    fontSize = 9
)

addStyle(
    wb,
    sheet = "Annotation_Source",
    style = wrap_style,
    rows = 2:(nrow(source_table) + 1),
    cols = 1:2,
    gridExpand = TRUE
)

saveWorkbook(
    wb,
    xlsx_file,
    overwrite = TRUE
)


# -----------------------------------------------------------------------------
# 19. Save a Word-compatible HTML table
#
# This can be opened directly in Microsoft Word and then saved as DOCX.
# -----------------------------------------------------------------------------

html_file <- file.path(
    output_dir,
    "DNMT3A_HM450_LEGACY_HG19_PROBE_ANNOTATION_PAPER_TABLE.html"
)

html_header <- paste0(
    "<html><head><meta charset='UTF-8'>",
    "<style>",
    "body{font-family:Arial;font-size:10pt;}",
    "table{border-collapse:collapse;width:100%;}",
    "th,td{border:1px solid #000;padding:4px;vertical-align:top;}",
    "th{font-weight:bold;text-align:center;}",
    "</style></head><body>",
    "<p><b>Supplementary Table. DNMT3A CpG probe annotation on the ",
    "Illumina HumanMethylation450 platform.</b></p>",
    "<p>Genomic coordinates are GRCh37/hg19 and were obtained from the ",
    "UCSC Xena TCGA legacy illuminaMethyl450_hg19_GPL16304_TCGAlegacy probeMap. ",
    "UCSC RefGene and CpG-island annotations were obtained from the ",
    "Illumina HumanMethylation450 annotation.</p>"
)

html_table <- paste0(
    "<table><thead><tr>",
    paste0(
        "<th>",
        names(paper_table),
        "</th>",
        collapse = ""
    ),
    "</tr></thead><tbody>"
)

for (i in seq_len(nrow(paper_table))) {

    row_values <- paper_table[i, ]

    row_values <- lapply(
        row_values,
        function(x) {
            x <- as.character(x)
            x[is.na(x)] <- ""
            x
        }
    )

    html_table <- paste0(
        html_table,
        "<tr>",
        paste0(
            "<td>",
            unlist(row_values),
            "</td>",
            collapse = ""
        ),
        "</tr>"
    )
}

html_table <- paste0(
    html_table,
    "</tbody></table>"
)

html_footer <- "</body></html>"

writeLines(
    c(
        html_header,
        html_table,
        html_footer
    ),
    con = html_file,
    useBytes = TRUE
)


# -----------------------------------------------------------------------------
# 20. Final console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("DNMT3A HM450 ANNOTATION COMPLETED\n")
cat("============================================================\n")

cat("\nPaper table CSV:\n")
cat(csv_file, "\n")

cat("\nPaper table Excel:\n")
cat(xlsx_file, "\n")

cat("\nWord-compatible HTML table:\n")
cat(html_file, "\n")

cat("\nFinal paper-table columns:\n")
print(
    names(paper_table)
)

cat(
    "\nIMPORTANT:\n",
    "The final chr and pos columns are hg38 coordinates from the exact ",
    "GDC/UCSC Xena HM450 probeMap used for the methylation mapping.\n",
    sep = ""
)

cat(
    "UCSC_RefGene_Name, UCSC_RefGene_Group, and Relation_to_Island ",
    "come from the official Illumina HumanMethylation450 annotation.\n"
)
