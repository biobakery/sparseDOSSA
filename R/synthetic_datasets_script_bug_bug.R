#' A function to generate a bug-bug correlation matrix with a
#'     fixed number of spikes.
#'
#' @param num_features is the number of features in the dataset
#'     (so the dimension of the correlation matrix)
#' @param num_spikes is the number of non-zero correlations to
#'     generate
#' @param log_corr_values is a scalar or vector of correlation values. If
#'     one value, all spikes will have that correlation value. If a vector
#'     of length > \code{num_features}, it will be sampled without
#'     replacement. If a vector of length < \code{num_features}, it will be
#'     recycled.
#' @return A list contains the microbe-microbe spiked abundance matrix and the
#'     information about the pairs of microbes that are correlated when only
#'     number of correlated microbes are specified.
#' @keywords internal
func_get_log_corr_mat_from_num <- function(num_features,
                                           num_spikes,
                                           log_corr_values) {
  message("start func_get_corr_mat_from_num")
  if (length(log_corr_values) == 1) {
    log_corr_values <- rep(log_corr_values, num_spikes)
  } else if (length(log_corr_values) > num_spikes) {
    log_corr_values <- sample(log_corr_values, num_spikes)
  } else if (length(log_corr_values) < num_spikes) {
    log_corr_values_rep <- rep(log_corr_values,
                               ceiling(num_spikes / length(log_corr_values)))
    log_corr_values <- log_corr_values_rep[1:num_spikes]
  }
  
  if (num_spikes >= num_features) {
    stop(
      paste0(
        "num_spikes (",
        num_spikes,
        ") >= num_features (",
        num_features,
        "); cannot spike any association"
      )
    )
  }
  
  if (num_spikes > 0) {
    domain_features <- sample(num_features, num_spikes)
    possible_range_features <-
      setdiff(seq_len(num_features), domain_features)
    
    range_features <- sample(possible_range_features, num_spikes)
    
    log_corr_mat <- func_generate_spike_structure(
      range_features = range_features,
      domain_features = domain_features,
      num_features = num_features,
      log_corr_values = log_corr_values
    )
    
    range_vec <- unique(range_features)
    ss <- lapply(range_vec, function(r) {
      list(
        domain = domain_features[range_features == r],
        range = r,
        log_corr_values = log_corr_values[range_features == r]
      )
    })
    domain_features_list <- lapply(ss, '[[', 'domain')
    log_corr_values_list <- lapply(ss, '[[', "log_corr_values")
    
    num_associations_str <- as.character(num_spikes)
    num_domain_features_str <-
      func_list_to_str(lapply(domain_features_list, length))
    range_features_str <- paste(range_vec, collapse = "; ")
    domain_features_str <- func_list_to_str(domain_features_list)
    assoc_dir_str <-
      func_list_to_str(lapply(log_corr_values_list, sign))
    log_corr_values_str <- func_list_to_str(log_corr_values_list)
    
    to_return <-
      list(mdLogCorr = log_corr_mat,
           RangeBugs = range_features,
           DomainBugs = domain_features)
  } else {
    to_return <- list(mdLogCorr = diag(num_features))
    num_associations_str <- "0"
    num_domain_features_str <- "NA"
    range_features_str <- "NA"
    domain_features_str <- "NA"
    assoc_dir_str <- "NA"
    log_corr_values_str <- "NA"
  }
  param_list <- list(
    strNumberOfAssociations = num_associations_str,
    strNumberCorrDomainBugs = num_domain_features_str,
    strIdxCorrRangeBugs = range_features_str,
    strIdxCorrDomainBugs = domain_features_str,
    strDirAssociations = assoc_dir_str,
    strLogCorrValues = log_corr_values_str
  )
  
  message("end func_get_corr_mat_from_num")
  return(c(to_return, param_list))
}

#' Get the log-basis correlation matrix from a file
#'
#' @inheritParams func_get_log_corr_mat_from_num
#' @param file_name The name of the file where the correlation values are
#'   stored. Should have fields `Domain`, `Range`, and `Correlation`.
#' @return A list contains the microbe-microbe spiked abundance matrix and the
#'     information about the pairs of microbes that are correlated when the pairs
#'     are specified through a user provided file.
#' @keywords internal
func_get_log_corr_mat_from_file <-
  function(num_features, file_name) {
    message("start func_get_corr_mat_from_file")
    corr_data <-
      read.table(file_name, stringsAsFactors = FALSE, header = TRUE)
    
    domain_features <- as.numeric(corr_data$Domain)
    range_features  <- as.numeric(corr_data$Range)
    log_corr_values <- as.numeric(corr_data$Correlation)
    log_corr_mat <- func_generate_spike_structure(
      range_features = range_features,
      domain_features = domain_features,
      num_features = num_features,
      log_corr_values = log_corr_values
    )
    
    range_vec <- unique(range_features)
    ss <- lapply(range_vec, function(r) {
      list(domain = domain_features[range_features == r],
           range = r,
           log_corr_values = log_corr_values[range_features == r])
    })
    domain_features_list <- lapply(ss, '[[', 'domain')
    log_corr_values_list <- lapply(ss, '[[', "log_corr_values")
    
    num_domain_features_str <-
      func_list_to_str(lapply(domain_features_list, length))
    range_features_str <- paste(range_vec, collapse = "; ")
    domain_features_str <- func_list_to_str(domain_features_list)
    assoc_dir_str <-
      func_list_to_str(lapply(log_corr_values_list, sign))
    log_corr_values_str <- func_list_to_str(log_corr_values_list)
    
    to_return <-
      list(mdLogCorr = log_corr_mat,
           RangeBugs = range_features,
           DomainBugs = domain_features)
    
    param_list <-
      list(
        strNumberOfAssociations = as.character(nrow(corr_data)),
        strNumberCorrDomainBugs = num_domain_features_str,
        strIdxCorrRangeBugs = range_features_str,
        strIdxCorrDomainBugs = domain_features_str,
        strDirAssociations = assoc_dir_str,
        strLogCorrValues = log_corr_values_str
      )
    
    message("end func_get_corr_mat_from_file")
    return(c(to_return, param_list))
  }

#' Convert a list of vectors into an appropriately-delimited single string
#'
#' @param mylist A list of vectors
#' @return A string concatenates all elements in a list by semi-colons.
#' @keywords internal
func_list_to_str <- function(mylist) {
  list_elts <- unlist(lapply(mylist, paste, collapse = ","))
  paste(list_elts, collapse = "; ")
}

#' Generate the spike structure as a list, and the correlation matrix
#'
#' @param range_features A vector of features which will be replaced with
#'   spiked data
#' @param domain_features A vector which contains the features that will be
#'   used to generate the spiked data
#' @inheritParams func_get_log_corr_mat_from_num
#' @return A correlation matrix corresponding to the spike structure.
#' @keywords internal
func_generate_spike_structure <-
  function(range_features,
           domain_features,
           num_features,
           log_corr_values) {
    message("start func_generate_spike_structure")
    log_corr_mat <- diag(num_features)
    log_corr_mat[cbind(range_features, domain_features)] = log_corr_values
    log_corr_mat[cbind(domain_features, range_features)] = log_corr_values
    # for (s in seq_along(range_features)){
    #   log_corr_mat[range_features[s],  domain_features[s]] <- log_corr_values[s]
    #   log_corr_mat[domain_features[s], range_features[s]]  <- log_corr_values[s]
    # }
    
    message("end func_generate_spike_structure")
    return(log_corr_mat)
  }
