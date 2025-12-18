edit_utils <- function() {
  file.edit("utils.R")
}

head1 <- function(..., n = c(6, 6)) {
  head(..., n = n)
}

get_manifest <- function(platform = c("EPICv1", "EPICv2", "MSA"), root = here::here()) {
  platform <- unique(match.arg(platform, several.ok = TRUE))
  path <- glue::glue("{root}/data/manifest/{platform}/{platform}_clean.pq")
  stopifnot(all(file.exists(path)))
  return(as.character(path))
}

get_DNAm <- function(platform = NULL, chr = NULL, root = here::here()) {
  if (is.null(platform)) {
    message("Available platforms:")
    print(GSE_meta)
    return(invisible(NULL))
  }
  if (!platform %in% GSE_meta$platform) {
    message("Available platforms:")
    print(GSE_meta)
    return(invisible(NULL))
  }

  raw <- if (platform == "EPICv1") {
    file.path(root, "GSE286313", "EPICv1")
  } else if (platform == "EPICv2") {
    file.path(root, "GSE286313", "EPICv2")
  } else if (platform == "MSA") {
    file.path(root, "GSE264438")
  } else {
    stop("platform must be either 'EPICv1', 'EPICv2', or 'MSA'")
  }
  if (!dir.exists(raw)) {
    stop("Directory not found: ", raw)
  }
  all_files <- list.files(raw, full.names = FALSE)

  if (is.null(chr)) {
    message("Available chromosomes:")
    chr_pattern <- "chr[0-9XYM]+"
    chr_matches <- stringr::str_extract(all_files, chr_pattern)
    available_chrs <- unique(na.omit(chr_matches))
    if (length(available_chrs) > 0) {
      print(sort(available_chrs))
      message("Use chr = 'all' to get all chromosomes")
    } else {
      message("No chromosome-based files found. All files:")
      print(all_files)
    }
    return(invisible(NULL))
  }

  if (length(chr) == 1 && tolower(chr) == "all") {
    chr_pattern <- "_chr[0-9XYM]+\\."
    chr_files <- all_files[stringr::str_detect(all_files, chr_pattern)]
    if (length(chr_files) == 0) {
      message("No chromosome-split files found in ", raw)
      return(tibble::tibble(file = character(0), path = character(0)))
    }
    return(tibble::tibble(
      file = chr_files,
      path = file.path(raw, chr_files)
    ))
  }

  # Handle specific chr(s)
  chr <- as.character(chr)
  chr <- ifelse(stringr::str_detect(chr, "^chr"), chr, paste0("chr", chr))
  patterns <- paste0("_", chr, "\\.")
  combined_pattern <- paste(patterns, collapse = "|")
  matching_files <- all_files[stringr::str_detect(all_files, combined_pattern)]
  if (length(matching_files) == 0) {
    message("No files matching chr(s) ", paste(chr, collapse = ", "), " in ", raw)
    return(tibble::tibble(file = character(0), path = character(0)))
  }
  return(tibble::tibble(
    file = matching_files,
    path = file.path(raw, matching_files)
  ))
}


if (!file.exists(here::here("data", "metadata.csv"))) {
  stop("The metadata file can be generated from `00.99.metadata.qmd`")
}

GSE_meta <- suppressMessages(readr::read_csv(here::here("data", "metadata.csv")))
