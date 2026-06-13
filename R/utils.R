# Internal utilities shared across the package.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

deepseekcell_version <- function() {
  "0.1.0"
}

is_blank <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    !nzchar(trimws(as.character(x[1])))
}

split_marker_text <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(character())
  }

  genes <- unlist(strsplit(as.character(x), "[,;[:space:]]+"))
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  unique(genes)
}

normalize_marker_list <- function(markers, max_genes = 30) {
  if (!is.list(markers) || length(markers) == 0) {
    stop("markers must be a non-empty list.", call. = FALSE)
  }

  if (!is.numeric(max_genes) || length(max_genes) != 1 || is.na(max_genes) || max_genes < 1) {
    stop("max_genes must be a positive numeric scalar.", call. = FALSE)
  }

  if (is.null(names(markers)) || any(!nzchar(trimws(names(markers))))) {
    names(markers) <- paste0("Cluster", seq_along(markers))
  }

  names(markers) <- make.unique(trimws(names(markers)))

  markers <- lapply(markers, function(x) {
    genes <- if (is.character(x) && length(x) == 1) {
      split_marker_text(x)
    } else {
      unique(trimws(as.character(x)))
    }

    genes <- genes[!is.na(genes) & nzchar(genes)]
    head(unique(genes), max_genes)
  })

  markers[lengths(markers) > 0]
}

as_confidence <- function(x, default = 0.5) {
  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- default
  out[out > 1 & out <= 100] <- out[out > 1 & out <= 100] / 100
  pmin(pmax(out, 0), 1)
}

as_flag <- function(x, default = FALSE) {
  if (is.logical(x)) {
    out <- x
  } else if (is.numeric(x)) {
    out <- x != 0
  } else {
    normalized <- tolower(trimws(as.character(x)))
    out <- normalized %in% c("true", "t", "yes", "y", "1", "mixed", "doublet")
    out[normalized %in% c("false", "f", "no", "n", "0", "single", "")] <- FALSE
  }

  out[is.na(out)] <- default
  out
}

extract_json_payload <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(NULL)
  }

  cleaned <- trimws(as.character(text))
  cleaned <- gsub("^```(?:json)?\\s*", "", cleaned, perl = TRUE)
  cleaned <- gsub("\\s*```$", "", cleaned, perl = TRUE)

  object_payload <- .extract_balanced_json(cleaned, "{", "}")
  if (!is.null(object_payload)) {
    return(object_payload)
  }

  .extract_balanced_json(cleaned, "[", "]")
}

.extract_balanced_json <- function(text, open_char, close_char) {
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  start <- which(chars == open_char)[1]

  if (is.na(start)) {
    return(NULL)
  }

  depth <- 0
  in_string <- FALSE
  escaped <- FALSE

  for (i in seq.int(start, length(chars))) {
    ch <- chars[[i]]

    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, "\"")) {
        in_string <- FALSE
      }
      next
    }

    if (identical(ch, "\"")) {
      in_string <- TRUE
      next
    }

    if (identical(ch, open_char)) {
      depth <- depth + 1
    } else if (identical(ch, close_char)) {
      depth <- depth - 1
      if (depth == 0) {
        return(paste(chars[start:i], collapse = ""))
      }
    }
  }

  NULL
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

first_env <- function(candidates) {
  for (candidate in candidates) {
    value <- Sys.getenv(candidate, unset = "")
    if (nzchar(value)) {
      return(value)
    }
  }

  ""
}
