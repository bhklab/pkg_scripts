#' load, convert and process the GDS data into a gDR object
#'
#' @param qcs QCS object exported from GeneDataScreenR
#' @param Gnumbers vector of Gnumber that should be exported from QCS
#' @param clids vector of CLids that should be exported from QCS
#' @param exclude_combo boolean indicating whether or not to exclude combination data.
#' Defaults to \code{FALSE}.
#' @param base_url base_url of REST API
#'
#' @return SummarizedExperiment
#'
#' @export
#'
convert_GDS_qcs_to_gDR_SE <- function(qcs,
                                      Gnumbers = NULL,
                                      clids = NULL,
                                      exclude_combo = FALSE,
                                      base_url = gDRwrapper::getConfig("base_url")) {

  # extract data and convert columns names
  df_ <- process_normdata_from_qcs(qcs,
                                   Gnumbers = Gnumbers,
                                   clids = clids,
                                   exclude_combo = exclude_combo)
  futile.logger::flog.warn("Reading '%i' wells from '%s'\n", nrow(df_), id(qcs))

  # correct Gnumbers in df_
  df_ <- correct_gnumbers(df = df_, base_url = base_url)

  # clean up to fit gDR format
  df_normed <- gDRcore::cleanup_metadata(df_)

  # clean up NA in CellLineName
  df_normed$CellLineName <- ifelse(is.na(df_normed$CellLineName), df_normed$clid, df_normed$CellLineName)
  
  # remove unnecessary cols in single-agent data
  df_normed <- cleanup_single_agent_data(df_normed)
  
  # set nested identifiers
  data_type <- gDRcore::data_model(df_normed)
  
  nested_identifiers <- if ((data_type == "combo" && !is(qcs, "MatrixQCSession")) || data_type == "single-agent") {
    gDRutils::get_env_identifiers("concentration")
  } else if (data_type == "combo" && is(qcs, "MatrixQCSession")) {
    unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"),
                                         simplify = FALSE))
  } else {
    stop("Unrecognized data type")
  }
  
  # create SE and process the data as gDR data
  se <- convert_df_GDSnormed_to_SE(df_normed, qcs, nested_identifiers = nested_identifiers)
  futile.logger::flog.warn("Data converted and normalized
                           into a '%i' by '%i' SummarizedExperiment object", nrow(se), ncol(se))


  se <- gDRcore::average_SE(se, series_identifiers = nested_identifiers)
  futile.logger::flog.warn("Averaging successful")
  
  se <- if (is(qcs, "MatrixQCSession")) {
    gDRcore::fit_SE.combinations(se = se, series_identifiers = nested_identifiers)
  } else {
    gDRcore::fit_SE(se = se, nested_identifiers = nested_identifiers)
  }
  futile.logger::flog.warn("Fitting successful")

  # replace fit results with original metrics
  se <- add_GDS_metrics(se, qcs)

  se <- gDRcore::add_codrug_group_SE(se)
  futile.logger::flog.warn("Co-treatment identification successful")
  gDRutils::validate_SE(se, expect_single_agent = exclude_combo)
  se # ready to be used or pushed to the database
}


#' process_normdata_from_qcs: get the well-level data from a qcs and
#'   format it as a data.frame comparible with gDR
#'
#' @param qcs QCS object exported from GeneDataScreenR
#' @param Gnumbers vector of Gnumber that should be exported from QCS
#' @param clids vector of CLids that should be exported from QCS
#' @param exclude_combo boolean
#'
#' @return data.frame with well-level data
#' @export
#'
process_normdata_from_qcs <- function(qcs,
                                      Gnumbers = NULL,
                                      clids = NULL,
                                      exclude_combo = FALSE) {

  # columns headers
  colnames_selected_substitutions <- list(
    c("series.CLID", "clid"),     # project <= 38 #nolint
    c("CLID", "clid"),            # project 39 #nolint
    c("series.Drug1", "Gnumber"), # project <= 38 #nolint
    c("d1_id", "Gnumber"),        # project 39 #nolint
    c("d1_conc", "Concentration"),
    c("d2_conc", "Concentration_2"),
    c("series.Drug2", "Gnumber_2"), # project <= 38 #nolint
    c("d2_id", "Gnumber_2"),        # project >=39 #nolint
    c("activities", "RelativeViability", function(x) 1 + (as.numeric(x) / 100)),
    # converts activities (range -100 to 0) to RV (range 0 to 1)
    c("masked", "masked"),
    c("compoundId", "compoundId"),   #  --> create an issue because assigned as row
    # identifier whereas it is a 'unit' identifier
    c("", "Duration")
  )

  # getting the project number based on experiment name
  exp_name <- name(qcs)
  project_number <- as.numeric(substr(exp_name,
                                     attr(regexpr("^P[A-Za-z]*", exp_name), "match.length") + 1,
                                     attr(regexpr("^P[A-Za-z]*\\d\\d?_", exp_name), "match.length") - 1))

  # get the data
  df_ <- as.data.frame(wells(qcs))

  # convert into hours for gDR
  df_$Duration <- as.numeric(duration(qcs)) * 24

  # checks on the data
  if (all(df_$masked) || all(is.na(df_$activities)) ||
      all(df_$activities %in% c("NA", "na", "NaN"))) {
    return(df_[FALSE, ]) # empty dataframe
  }

  # rename columns and select the ones of interest
    for (sel_sub in colnames_selected_substitutions) {
    if (sel_sub[[1]] %in% colnames(df_) && !(sel_sub[[2]] %in% colnames(df_))) {
      colnames(df_)[colnames(df_) == sel_sub[[1]]] <- sel_sub[[2]]
      if (length(sel_sub) == 3) {
        df_[, sel_sub[[2]]] <- sel_sub[[3]](df_[, sel_sub[[2]]]) # for activities conversion
      }
    }
  }

  if (exclude_combo) {
    df_ <- df_[df_$Concentration_2 %in% "NA" |
      is.na(df_$Concentration_2), ]
  }

  # remove data with NA in Gnumbers or in clid
  df_[df_ == "NA"] <- NA
  df_$clid <- ifelse(is.na(df_$clid), df_$CellLine, df_$clid)
  df_ <- df_[!is.na(df_$clid) & !is.na(df_$Gnumber), ]
  
  # convert the Gnumbers
  df_$Gnumber <- as.character(df_$Gnumber) %>%
    ifelse(is.na(stringr::str_extract(., "G[0-9].*")), ., stringr::str_extract(., "G[0-9].*"))

  df_$Gnumber_2 <- as.character(df_$Gnumber_2) %>%
    ifelse(is.na(stringr::str_extract(., "G[0-9].*")), ., stringr::str_extract(., "G[0-9].*"))

  df_$clid <- as.character(df_$clid)
  df_$clid[grepl("^\\d+$", df_$clid)] <- paste0("CL", df_$clid[grepl("^\\d+$", df_$clid)])

  # Check the correction of Concentration/Concentration_2 format
  if (any(!grepl("\\d+\\.*\\d*$|ug", na.omit(c(df_$Concentration, df_$Concentration_2))))) {
    stop("Uncorrected format of Concentration/Concentration_2")
  }
  
  # convert the concentration in numeric values (from factors) and extract only numeric values
  df_$Concentration <- as.numeric(as.character(stringr::str_extract(df_$Concentration, "\\d+\\.*\\d*")))
  df_$Concentration_2 <- as.numeric(as.character(stringr::str_extract(df_$Concentration_2, "\\d+\\.*\\d*")))

  # select the data based on Gnumber and clids
  df_ <- df_[, intersect(sapply(colnames_selected_substitutions, "[[", 2), colnames(df_))]

  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")[1]

  # replace NA with untreated for second agent
  df_$Concentration_2[is.na(df_$Concentration_2)] <- 0

  # swap Gnumber_2 and Gnumber if Concentration == 0
  if (any(df_$Concentration %in% 0)) {
    df_$Gnumber[df_$Concentration %in% 0] <- untreated_tag
  }

  # deal with the issue of annotation in P22 if required:
  if (expFolder(qcs) == "/FunctionalGenomics/gCSI_Published/P22") {
    df_ <- add_cotreatment_P22(df_, qcs)
  }

  # filter data if selected
  if (!is.null(Gnumbers)) {
    df_ <- df_[
      (df_$Gnumber %in% c(Gnumbers, untreated_tag) &
         df_$Gnumber_2 %in% c(Gnumbers, untreated_tag)) |
        (substr(df_$Gnumber, 1, 9) %in% c(Gnumbers, untreated_tag) &
           substr(df_$Gnumber_2, 1, 9) %in% c(Gnumbers, untreated_tag)), ]
  }

  if (!is.null(clids)) {
    df_ <- df_[df_$clid %in% clids, ]
  }

  # swap for getting avoiding untreated as Gnumber (primary)
  if (any(df_$Gnumber == untreated_tag)) {
    temp_df <- df_[df_$Gnumber == untreated_tag, ]
    temp_df$Gnumber <- temp_df$Gnumber_2
    temp_df$Concentration <- temp_df$Concentration_2
    temp_df$Gnumber_2 <- untreated_tag
    temp_df$Concentration_2 <- 0
    df_ <- rbind(df_[df_$Gnumber != untreated_tag, ], temp_df)
  }

  # deal with replicated data (same clid/Gnumber but different compoundId)
  if (!("compoundId" %in% colnames(df_)))  {
    # nothing to do here, returning the data
    return(df_)
  }

  if ("Concentration_2" %in% colnames(df_)) {
    df_conditions <- unique(df_[df_$Concentration_2 == 0, c("compoundId", "clid", "Gnumber")])
  } else {
    df_conditions <- unique(df_[, c("compoundId", "clid", "Gnumber")])
  }

  dupl_id <- reshape2::melt(table(df_conditions[, -1]))
  dupl_id <- as.character(unique(dupl_id[dupl_id[, 3] > 1, "Gnumber"]))
  dupl_conditions <- df_[df_$Concentration_2 == 0 & df_$Gnumber %in% dupl_id, ]

  if (nrow(dupl_conditions) == 0) {
    # nothing to do here, returning the data
    return(df_)
  }

  dupl_conc_top_conc <- aggregate(dupl_conditions[, "Concentration", drop = FALSE],
    by = as.list(dupl_conditions[, c("compoundId", "Gnumber", "clid")]),
    max)
  sel_duplicate <- aggregate(dupl_conc_top_conc[, "Concentration", drop = FALSE],
    by = as.list(dupl_conc_top_conc[, c("Gnumber", "clid")]),
    function(x) length(unique(x)))

  # only remove one of the duplicates that have divergent top concentrations
  sel_duplicate <- sel_duplicate[sel_duplicate$Concentration > 1, c("Gnumber", "clid")]

  # select a single compoundId
  if (nrow(sel_duplicate) > 1) {
    for (i in seq(nrow(sel_duplicate))) {
      temp <- dupl_conc_top_conc[dupl_conc_top_conc$Gnumber == sel_duplicate$Gnumber[i] &
              dupl_conc_top_conc$clid == sel_duplicate$clid[i], ]
      # select the compoundId code (P##S##_*****) with the the highest values
      temp$compoundId_value <- as.numeric(gsub("^\\w*_", "", temp$compoundId))
      sel_duplicate$discarded_compoundId[i] <- setdiff(temp$compoundId,
              temp$compoundId[which.max(temp$compoundId_value)])
    }

    # remove the duplicated compoundId
    df_ <- df_[!(df_$compoundId %in% sel_duplicate$discarded_compoundId), ]

    # warnings
    futile.logger::flog.warn("%i drugs (%s) have multiple replicates with different top concentrations",
        length(unique(sel_duplicate$Gnumber)), paste(unique(sel_duplicate$Gnumber), collapse = ", "))
    futile.logger::flog.warn("%i compoundId for %s are discarded: \n\t%s",
        length(unique(sel_duplicate$discarded_compoundId)),
        gsub("_\\d*$", "_*", sel_duplicate$discarded_compoundId[1]),
        paste(gsub("^\\w*_", "", unique(sel_duplicate$discarded_compoundId)), collapse = ", "))
  }
  df_
}


add_cotreatment_P22 <- function(df_, qcs) {

  # add some checks: the df_ should be matching the qcs
  if (nrow(wells(qcs)) != nrow(df_) || max(abs(wells(qcs)$activities - 100 * (df_$RelativeViability - 1)))
      > 1e-10) {
    futile.logger::flog.error("Data in df_ seem to not match the data in qcs ; please check")
  }

  # P22_mapping_2020_0806.tsv is a copy of:
  # https://urldefense.com/v3/__https://docs.google.com/spreadsheets/d/1Di8GYbuKLrMIvrfcXK8nRzwZydK9jiLqi3pMBtti0bk/edit*gid=411979340__;Iw!!CjcC7IQ!bo9Ei7d3tS_DMJeIypFgnkrK5JbajV2yRdMCNAhs65W5Y1mpPyGkDC1ig7VITLwu8_apPvM5OdI$ [docs[.]google[.]com]
  df_mapping <- read.csv(system.file("P22_mapping_2020_0806.tsv", package = "gDRinternal"),
                  sep = "\t")
  # check if any information is available
  if (!(id(qcs) %in% df_mapping$Ref.ID.1)) {
      futile.logger::flog.warn("Identified a QCS from P22 but no mapping information found in reference file")
      return(df_)
  } else {
      futile.logger::flog.warn("Identified a QCS from P22: adding co-treatment information from mapping file")
  }

  df_mapping <- df_mapping[df_mapping$Ref.ID.1 == id(qcs), ]

  df_$clid <- df_mapping[match(wells(qcs)$CellLine,
                df_mapping$CL.name.used.in.CMT), "CLID"]

  df_$Gnumber_2 <- df_mapping[match(wells(qcs)$CellLine,
                df_mapping$CL.name.used.in.CMT), "drug.2.G."]
  df_$Gnumber_2[df_$Gnumber_2 %in% ""] <- "vehicle"

  df_$Concentration_2 <- df_mapping[match(wells(qcs)$CellLine,
                                          # `μM` is not recognized on Rosalind so instead
                                          # of subsetting by colnames let's use grep
                df_mapping$CL.name.used.in.CMT), grep("supplement", colnames(df_mapping))]
  df_$Concentration_2[is.na(df_$Concentration_2)] <- 0

  return(df_)
}


#' Convert a data.frame into a SE object at the normalized data stage
#'
#' @param df_normed data.frame will well-level data
#' @param qcs QCS object
#' @param nested_identifiers character vector of column names to include in the data.frames
#' in the assays of the resulting \code{SummarizedExperiment} object.
#'
#' @return initial SummarizedExperiment
#' @export
#'
convert_df_GDSnormed_to_SE <- function(df_normed, qcs = NULL, nested_identifiers = NULL) {

  df_raw_data <- df_normed

  # convert to characters (necessary for gDR operations)
  df_raw_data$DrugName <- as.character(df_raw_data$DrugName)
  df_raw_data$Gnumber <- as.character(df_raw_data$Gnumber)

  if ("DrugName_2" %in% colnames(df_raw_data)) {
    df_raw_data$DrugName_2 <- as.character(df_raw_data$DrugName_2)
    df_raw_data$Gnumber_2 <- as.character(df_raw_data$Gnumber_2)
  }

  # create dummy untreated for fitting gDR scheme
  df_ctrl <- unique(df_raw_data[, !(colnames(df_raw_data) %in%
        c("ReadoutValue", "BackgroundValue", "CorrectedReadout", "Concentration",
          "RelativeViability", "DrugName", "Gnumber", "masked"))])
  df_ctrl$RelativeViability <- 1
  df_ctrl$Concentration <- 0
  df_ctrl$Gnumber <- df_ctrl$DrugName <- df_ctrl$drug_moa <- gDRutils::get_env_identifiers("untreated_tag")[2]
  df_ctrl$masked <- FALSE

  if ("Gnumber_2" %in% colnames(df_raw_data)) {
    df_ctrl <- df_ctrl[df_ctrl$Gnumber_2 %in% gDRutils::get_env_identifiers("untreated_tag"), ]
  }

  if ("compoundId" %in% colnames(df_raw_data)) {
    df_ctrl$compoundId <- "ctrl"
  }

  # stack the data and contols
  df_raw_data <- data.table::rbindlist(list(df_raw_data, df_ctrl), use.names = TRUE, fill = TRUE)
  data.table::setDF(df_raw_data)
  colnames(df_raw_data)[colnames(df_raw_data) == "RelativeViability"] <- "ReadoutValue"
  df_raw_data$BackgroundValue <- 0

  # remove the compoundId column because it confused the assignment of row and col in the SE
  temp_df <- df_raw_data
  if ("compoundId" %in% colnames(df_raw_data)) {
    temp_df$compoundId <- NULL
  }
  
  # convert into a SE with the right formating (to test for combo experiments)
  # skip nested_confounders due to the lack of barcodes
  se <- gDRcore::create_and_normalize_SE(temp_df,
                                         nested_identifiers = nested_identifiers,
                                         nested_confounders = NULL)

  # add the qcs compoundId back
  S4Vectors::metadata(se)$df_raw_data <- df_raw_data

  if (!is.null(qcs)) {
    se <- add_qcs_info_SE(se, qcs)
  }

  se
}

add_qcs_info_SE <- function(se, qcs) {
  .Deprecated(new = "add_qcs_metadata_to_SE")
  add_qcs_metadata_to_SE(se = se, qcs = qcs)
}


#' add_qcs_info_SE
#'
#' @param se SummarizedExperiment with the response data
#' @param qcs QCS object from GeneDataScreenR
#'
#' @return SummarizedExperiment with added QCS metadata
#' @export
#'
add_qcs_metadata_to_SE <- function(se, qcs) {
  qcs_name <- name(qcs)
  qcs_id <- id(qcs)
  se <- gDRutils::set_SE_experiment_metadata(se, list(
    description = sprintf(paste("gCSI experiment %s (%s - %s); imported from folder %s.",
                                "Done in %s-well plate on %s (update %s)"),
                          qcs_name, qcs_id, expId(qcs), expFolder(qcs),
                          plateFormat(qcs), creationDate(qcs), updatedDate(qcs)),
    name = paste(qcs_id, qcs_name),
    experimentalist = creator(qcs)
  ))
  se
}


#' add_GDS_metrics
#'
#' @param se SummarizedExperiment with response data
#' @param qcs QCS object
#' @param capping_fold level of capping of new data
#'
#' @return SummarizedExperiment with GDS metrics
#'
#' @details Adding metrics are only supported for QCSessions of class \code{"CGSQCSession"} and \code{"gCSIQCSession"}.
#' @export
#'
add_GDS_metrics <-
  function(se,
           qcs,
           capping_fold = 5,
           ic50_cols = c("IC50", "meanviability.ic50"),
           mv_cols = c("MV", "meanviability.meanviability")) {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_multi_class(qcs, c("MatrixQCSession", "gCSIQCSession", "CGSQCSession"))
  checkmate::assert_number(capping_fold)

  # get the original values
  df_metrics <- treatments(qcs)

  # skip it if no "compoundId" which is used for Matrix synergy experiments
  if (!has_valid_compound_id(df_metrics)) {
    futile.logger::flog.warn("column `compoundId` not valid in treatments(qcs) --> metric replacement skipped")
    return(se)
  }

  df_raw_data <- S4Vectors::metadata(se)$df_raw_data
  metrics_mx <- assay(se, "Metrics")
  cdata <- colData(se)
  rdata <- rowData(se)
  for (i in rownames(metrics_mx)) {
    for (j in colnames(metrics_mx)) {
      nested <- metrics_mx[i, j][[1]]
      if (is.null(nested) || nrow(nested) == 0) {
        next
      }
      # get the compoundId for each condition in the SE
      rowIdx <- array(TRUE, nrow(df_raw_data))
      for (mdata in intersect(colnames(rdata),
                              colnames(df_raw_data))) {
        if (!is.na(rdata[i, mdata])) {
          rowIdx <- rowIdx &
            (df_raw_data[, mdata] %in% rdata[i, mdata])
        }
      }

      for (mdata in intersect(colnames(cdata),
                              colnames(df_raw_data))) {
        if (!is.na(cdata[j, mdata])) {
          rowIdx <- rowIdx &
            (df_raw_data[, mdata] %in% cdata[j, mdata])
        }
      }

      compoundId_match <- unique(df_raw_data$compoundId[rowIdx])

      original_metrics <- df_metrics[df_metrics$compoundId %in% as.character(compoundId_match) &
                                      !(df_metrics$fitModelName %in% "DRCInvalidFitResult"), , drop = FALSE]

      # skip replacement if this is a drug matrix experiment
      if (any(original_metrics$fitStrategy %in% "Synergy Model")) {
        next
      }

      # clean up if there is not one and only one match
      if (nrow(original_metrics) == 0) {
        futile.logger::flog.warn("no matching condition found in the qcs for replacing metrics for %s, %s", i, j)
        next

      } else if (nrow(original_metrics) > 1) {
        # get x_mean to display
        mv_metric <- intersect(names(original_metrics), mv_cols)
        MV_GDS <- if (length(mv_metric) > 0) {
          standardize_metric(original_metrics, mv_metric[1])[["v"]]
        } else {
          NA
        }
        # calc ic50 to display
        ic50_metric <- intersect(names(original_metrics), ic50_cols)
        ic50_GDS <- if (length(ic50_metric) > 0) {
          # add missing cols (just for display)
          if (!"max_conc" %in% colnames(original_metrics)) {
            original_metrics[["max_conc"]] <- max(10 ^ nested$maxlog10Concentration)
          }
          standardize_metric(original_metrics, ic50_metric[1])
        } else {
          NA
        }

        # post warning message
        futile.logger::flog.warn(paste0("More than one conditions found in the qcs for replacing metrics for %s, %s:\n",
            "\t MV_GDS = %s, MV_gDR = %.4f\n\t ic50_GDS = %s, ic50_gDR = %.4f\n --> ignored"),
            i, j, paste(MV_GDS, collapse = ", "), nested[nested$normalization_type == "RV", "x_mean"],
                paste(ic50_GDS, collapse = ", "), nested[nested$normalization_type == "RV", "xc50"])
        next
      }

      # create new row
      row_idx <- nrow(nested) + 1
      nested[row_idx, ] <- NA
      nested[row_idx, "fit_source"] <- "GDS"
      nested[row_idx, "normalization_type"] <- "RV"
      nested[row_idx, "maxlog10Concentration"] <- nested[1, "maxlog10Concentration"]
      nested[row_idx, "N_conc"] <- nested[1, "N_conc"]

      # replace x_mean
      mv_metric <- names(original_metrics)[names(original_metrics) %in% mv_cols]
      if (length(mv_metric) > 0) {
        nested[row_idx, "x_mean"] <- standardize_metric(original_metrics, mv_metric[1])[["v"]]
      }

      # replace x_AOC
      nested[row_idx, "x_AOC"] <- 1 - nested[row_idx, "x_mean"]
      original_metrics[["max_conc"]] <- 10 ^ nested$maxlog10Concentration[1]

      if (original_metrics$fitModelName %in% c("DRCConstantFitResult", "DRCInvalidFitResult")) {
        # need to handle the cases of a flat fit
        #   --> copied from the gDR function (TODO: make a separate function to be called for default values)
        nested[row_idx, "ec50"] <- 0
        nested[row_idx, "h"] <- 0.0001
        nested[row_idx, "xc50"] <- ifelse(nested[row_idx, "x_mean"] > .5,
                                                  Inf, -Inf)
        nested[row_idx, "x_inf"] <- nested[row_idx, "x_mean"]
        nested[row_idx, "x_0"] <- 1
        nested[row_idx, "fit_type"] <- "DRCConstantFitResult"

      } else {
        # there is a proper fit
        for (metric in get_gds_metrics()[["main"]]) {
          if (metric %in% names(original_metrics)) {
            ml <- standardize_metric(original_metrics, metric)
            nested[row_idx, ml[["n"]]] <- ml[["v"]]
          }
        }

        # Add xc50 = +/-Inf for any curves that do not reach RelativeViability = 0.5
        if (is.na(nested[row_idx, "xc50"])) {
          nested[row_idx, "xc50"] <- ifelse(nested[row_idx, "x_inf"] > 0.5, Inf, -Inf)
        }
      }

    metrics_mx[i, j][[1]] <- nested
    }
  }

  assay(se, "Metrics") <- metrics_mx
  se
}

#' Correct wrong Gnumbers in QCS
#'
#' @param df_ data frame produced by `process_normdata_from_qcs()`
#' @param base_url base_url of REST API
#'
#' @return data.frame with corrected gnumbers
#' @export
#'
#' @examples
correct_gnumbers <- function(df,
                             base_url = gDRwrapper::getConfig("base_url")) {
  drugsCols <- grep("Gnumber", colnames(df), value = TRUE)
  drugsDB <- gDRwrapper::get_drugs(base_url)
  for (col in drugsCols) { # remove non-alphanumeric signs from the end of gnumbers
    df[[col]] <- gsub("[^[:alnum:] ]$", "", df[[col]])
  }
  for (col in drugsCols) {
    df[[col]] <- unlist(lapply(df[[col]], function(x) {
      if (grepl("^G[0-9]+", x) || x %in% gDRutils::get_env_identifiers("untreated_tag")) {
        return(x)
      } else {
        corrected <- drugsDB$gnumber[agrep(gsub("\\..*", "", x), drugsDB$drug_name, ignore.case = TRUE)[1]]
        ifelse(is.na(corrected), x, corrected)
      }
    }))
  }
  df
}


#' validate compound_id data
#'
#' check if given `S4Vectors::DataFrame` has valid compound_id data
#' @param compound_id_col string with the name of the column to be checked
#' @param compound_regex string with the regex to be used for data validation
#' @return logical flag, TRUE if column with compound-id data is present and all
#' data is passing the expected regex, FALSE otherwise
#' @export
has_valid_compound_id <-
  function(df_metrics,
           compound_id_col = "compoundId",
           compound_id_regex = "^P\\d+S\\d+_\\d*$") {
    checkmate::assert_true(inherits(df_metrics, c("DataFrame", "data.frame")))
    checkmate::assert_string(compound_id_col)
    checkmate::assert_string(compound_id_regex)

    compound_id_col %in% colnames(df_metrics) &&
      all(grepl(compound_id_regex, df_metrics[[compound_id_col]]))
  }

#' @export
calculate_ic50 <-
  function(x,
           var_name,
           log_factor = 6,
           capping_fold = 5,
           ic50_cols = c("IC50", "meanviability.ic50")) {

  checkmate::assert_true(inherits(x, c("DataFrame", "data.frame")))
  checkmate::assert_string(var_name)
  checkmate::assert_number(log_factor)
  checkmate::assert_number(capping_fold)
  checkmate::assert_character(ic50_cols)
  checkmate::assert_choice(var_name, ic50_cols)

  # gladkia: based on the code I assume that following columns are required in 'x'
  req_cols <- c(var_name, "fitModelName", "qAC50Mode", "SInf", "S0", "linearAC50", "hillCoefficient", "max_conc")
    if (length(setdiff(req_cols, colnames(x))) > 0) {
      errMsg <-
        sprintf("following columns have not been found in 'x' but are required: '%s'",
                toString(setdiff(req_cols, colnames(x))))
      stop(errMsg)
    }

    if (x[["fitModelName"]] == "DRCConstantFitResult") {
      return(gsub("<", -Inf, gsub(">", Inf, x[["qAC50Mode"]])))
    }

    qic50 <- as.numeric(x[[var_name]]) * (10 ^ log_factor)
    if (is.na(qic50)) {
      if (as.numeric(x[["SInf"]]) > -50) {
        qic50 <- Inf
      } else if (as.numeric(x[["S0"]]) < -50) {
        qic50 <- -Inf
      } else {
        # multiply by 1e6 because GDS uses M as unit whereas we use µM in gDR
        qic50 <- 1e6 * as.numeric(x[["linearAC50"]]) * (((as.numeric(x[["S0"]]) / 100) -
                                                           (as.numeric(x[["SInf"]]) / 100)) /
                                                          (-0.5 - (as.numeric(x[["SInf"]]) / 100)) - 1) ^
          (1 / as.numeric(x[["hillCoefficient"]]))
      }
    }

    # the capping is a bit arbitrary. Ideally, it should be based on the lowest tested concentration
    # but we don't necessarily record that, so I put 5 orders below the highest dose
    qic50 <- ifelse(qic50 > as.numeric(x[["max_conc"]]) * capping_fold, Inf, qic50)  # capping
    qic50 <- ifelse(qic50 < as.numeric(x[["max_conc"]]) / (capping_fold * 1e5), -Inf, qic50) # capping
    qic50
  }

#' @export
standardize_metric <- function(x, metric) {
  checkmate::assert_true(inherits(x, c("DataFrame", "data.frame")))
  checkmate::assert_choice(metric, get_gds_metrics()[["main"]])

  switch(metric,
    linearAC50 =  list(v = as.numeric(x[["linearAC50"]]) * 1e6, n = "ec50"),
    S0 = list(v = 1 + (as.numeric(x[["S0"]]) / 100), n = "x_0"),
    SInf = list(v = 1 + (as.numeric(x[["SInf"]]) / 100), n = "x_inf"),
    hillCoefficient = list(v = as.numeric(x[["hillCoefficient"]]), n = "h"),
    meanviability.meanviability = list(v = x[["meanviability.meanviability"]], n = "x_mean"),
    MV = list(v = x[["MV"]], n = "x_mean"),
    meanviability.ic50 = list(v = calculate_ic50(x, "meanviability.ic50"), n = "xc50"),
    IC50 = list(v = calculate_ic50(x, "IC50"), n = "xc50"),
    fitModelName = list(v = as.character(x[["fitModelName"]]), n = "fit_type"),
    rSquare = list(v = x[["rSquare"]], n = "r2")
  )
}

#' @export
get_gds_metrics <- function() {
  list(
    main = c(
      "linearAC50",
      "S0",
      "SInf",
      "hillCoefficient",
      "meanviability.meanviability",
      "MV",
      "meanviability.ic50",
      "IC50",
      "fitModelName",
      "rSquare"
    ),
    auxiliary = c("qAC50Mode",
                  "max_conc")
  )
}

#' @export
cleanup_single_agent_data <- function(df) {
  combo_cols <- grep("_2", colnames(df), value = TRUE)
  if (all(df[gDRutils::get_env_identifiers("concentration2")] == 0)) {
    df[, combo_cols] <- NULL
  }
  df
}
