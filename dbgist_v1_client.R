#!/usr/bin/env Rscript

# Minimal R client for dbGIST unified API v1.
# Intended to be sourced directly:
#
#   source("https://raw.githubusercontent.com/youngfly93/test_repo/main/dbgist_v1_client.R")
#   client <- dbgist_client("https://www.dbgist.com")
#   client$health()
#   client$transcriptomics_summary("KIT")

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
})

dbgist_api_error <- function(status, code, message, payload = NULL) {
  structure(
    list(status = status, code = code, message = message, payload = payload),
    class = c("dbgist_api_error", "error", "condition")
  )
}

conditionMessage.dbgist_api_error <- function(c) {
  sprintf("dbGIST API error [%s %s]: %s", c$status, c$code, c$message)
}

dbgist_request <- function(base_url, path, method = "GET", payload = NULL, timeout = 120) {
  url <- paste0(sub("/+$", "", base_url), if (startsWith(path, "/")) path else paste0("/", path))

  req <- request(url) |>
    req_method(method) |>
    req_timeout(timeout) |>
    req_options(noproxy = "*") |>
    req_headers(Accept = "application/json") |>
    req_error(is_error = function(resp) FALSE)

  if (!is.null(payload)) {
    req <- req |>
      req_headers(`Content-Type` = "application/json") |>
      req_body_json(payload, auto_unbox = TRUE)
  }

  resp <- tryCatch(
    req_perform(req),
    error = function(e) {
      stop(
        dbgist_api_error(
          status = 503,
          code = "NETWORK_ERROR",
          message = conditionMessage(e),
          payload = NULL
        )
      )
    }
  )

  body_text <- resp_body_string(resp)
  parsed <- tryCatch(
    fromJSON(body_text, simplifyVector = FALSE),
    error = function(e) NULL
  )

  if (is.null(parsed) || !is.list(parsed)) {
    stop(
      dbgist_api_error(
        status = resp_status(resp),
        code = "INVALID_RESPONSE",
        message = "Response is not a JSON object",
        payload = NULL
      )
    )
  }

  if (resp_status(resp) >= 400 || identical(parsed$ok, FALSE)) {
    err <- parsed$error %||% list()
    stop(
      dbgist_api_error(
        status = resp_status(resp),
        code = err$code %||% "API_ERROR",
        message = err$message %||% "Unknown API error",
        payload = parsed
      )
    )
  }

  parsed
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_features <- function(features) {
  vals <- trimws(as.character(features))
  vals[nzchar(vals)]
}

dbgist_client <- function(base_url = "http://127.0.0.1:8000", timeout = 120) {
  module_get <- function(module, endpoint) {
    dbgist_request(base_url, sprintf("/api/v1/%s/%s", module, endpoint), timeout = timeout)
  }

  module_post <- function(module, endpoint, payload) {
    dbgist_request(
      base_url,
      sprintf("/api/v1/%s/%s", module, endpoint),
      method = "POST",
      payload = payload,
      timeout = timeout
    )
  }

  list(
    base_url = base_url,
    health = function() dbgist_request(base_url, "/api/v1/health", timeout = timeout),
    modules = function() dbgist_request(base_url, "/api/v1/modules", timeout = timeout),

    module_health = function(module) module_get(module, "health"),
    module_ready = function(module) module_get(module, "ready"),
    module_capabilities = function(module) module_get(module, "capabilities"),

    transcriptomics_health = function() module_get("transcriptomics", "health"),
    transcriptomics_ready = function() module_get("transcriptomics", "ready"),
    transcriptomics_capabilities = function() module_get("transcriptomics", "capabilities"),
    transcriptomics_summary = function(feature) {
      module_post("transcriptomics", "summary", list(feature = feature))
    },
    transcriptomics_features_check = function(features) {
      module_post("transcriptomics", "features/check", list(features = normalize_features(features)))
    },
    transcriptomics_clinical = function(feature) {
      module_post("transcriptomics", "analysis/clinical", list(feature = feature))
    },
    transcriptomics_survival = function(feature, survival_type = "OS", cutoff = "Auto") {
      module_post(
        "transcriptomics",
        "analysis/survival",
        list(feature = feature, type = survival_type, cutoff = cutoff)
      )
    },
    transcriptomics_drug_response = function(feature) {
      module_post("transcriptomics", "analysis/drug-response", list(feature = feature))
    },

    noncoding_health = function() module_get("noncoding", "health"),
    noncoding_ready = function() module_get("noncoding", "ready"),
    noncoding_capabilities = function() module_get("noncoding", "capabilities"),
    noncoding_summary = function(feature) {
      module_post("noncoding", "summary", list(feature = feature))
    },
    noncoding_features_check = function(features) {
      module_post("noncoding", "features/check", list(features = normalize_features(features)))
    },
    noncoding_clinical = function(feature) {
      module_post("noncoding", "analysis/clinical", list(feature = feature))
    },
    noncoding_drug_response = function(feature) {
      module_post("noncoding", "analysis/drug-response", list(feature = feature))
    },

    proteomics_health = function() module_get("proteomics", "health"),
    proteomics_ready = function() module_get("proteomics", "ready"),
    proteomics_capabilities = function() module_get("proteomics", "capabilities"),
    proteomics_summary = function(feature) {
      module_post("proteomics", "summary", list(feature = feature))
    },
    proteomics_features_check = function(features) {
      module_post("proteomics", "features/check", list(features = normalize_features(features)))
    },
    proteomics_clinical = function(feature) {
      module_post("proteomics", "analysis/clinical", list(feature = feature))
    },
    proteomics_survival = function(feature, survival_type = "OS", cutoff = "Auto") {
      module_post(
        "proteomics",
        "analysis/survival",
        list(feature = feature, type = survival_type, cutoff = cutoff)
      )
    },
    proteomics_drug_response = function(feature) {
      module_post("proteomics", "analysis/drug-response", list(feature = feature))
    },

    phosphoproteomics_health = function() module_get("phosphoproteomics", "health"),
    phosphoproteomics_ready = function() module_get("phosphoproteomics", "ready"),
    phosphoproteomics_capabilities = function() module_get("phosphoproteomics", "capabilities"),
    phosphoproteomics_summary = function(feature) {
      module_post("phosphoproteomics", "summary", list(feature = feature))
    },
    phosphoproteomics_features_check = function(features) {
      module_post("phosphoproteomics", "features/check", list(features = normalize_features(features)))
    },
    phosphoproteomics_clinical = function(feature) {
      module_post("phosphoproteomics", "analysis/clinical", list(feature = feature))
    },
    phosphoproteomics_survival = function(feature, survival_type = "OS", cutoff = "Auto") {
      module_post(
        "phosphoproteomics",
        "analysis/survival",
        list(feature = feature, type = survival_type, cutoff = cutoff)
      )
    },

    singlecell_health = function() module_get("singlecell", "health"),
    singlecell_ready = function() module_get("singlecell", "ready"),
    singlecell_capabilities = function() module_get("singlecell", "capabilities"),
    singlecell_summary = function(feature) {
      module_post("singlecell", "summary", list(feature = feature))
    },
    singlecell_features_check = function(features) {
      module_post("singlecell", "features/check", list(features = normalize_features(features)))
    },
    singlecell_clinical = function(feature) {
      module_post("singlecell", "analysis/clinical", list(feature = feature))
    },
    singlecell_submit_job = function(feature, mode = "full", analysis = NULL) {
      payload <- list(feature = feature, mode = mode)
      if (!is.null(analysis)) payload$analysis <- analysis
      module_post("singlecell", "jobs", payload)
    },
    singlecell_job_status = function(job_id) {
      module_get("singlecell", sprintf("jobs/%s", job_id))
    },
    singlecell_job_result = function(job_id) {
      module_get("singlecell", sprintf("jobs/%s/result", job_id))
    }
  )
}
