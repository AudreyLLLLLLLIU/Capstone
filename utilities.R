
#' Resolve file or URL
#'
#' @param path File path or URL.
#'
#' @details If \code{path} is a local file path, assert that it
#'   exists, then just return it. If \code{path} is a URL, see if it
#'   has been downloaded before by checking if the file
#'   \code{basename(path)} exists. If it exists, just return
#'   \code{basename(path)}. Otherwise, download the file at
#'   \code{path} to \code{basename(path)} and return
#'   \code{basename(path)}.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom urltools scheme
#' @importFrom utils download.file
download_if_url <- function(path) {
  assert_that(is.string(path))
  
  is_url <- function(path) {
    grepl("^https?://", path)
  }
  
  if (!is_url(path)) {
    assert_that(file.exists(path))
    return(path)
  }
  
  if (file.exists(basename(path))) {
    message(sprintf("Skipping download: file %s already exists", basename(path)))
    return(basename(path))
  }
  
  message(sprintf("Downloading to %s...", basename(path)))
  download.file(path, basename(path), "auto")
  return(basename(path))
}


#' Decompress a gzipped tarball
#'
#' @param path File path to gzipped tarball.
#'
#' @details The contents of the tarball are extracted into a directory which has
#'   directory name equals to the basename of the \code{path}. For example,
#'   \code{a.tar.gz} containing \code{b.txt} will generate directory \code{./a/}
#'   and file \code{./a/b.txt}.
#'
#' @importFrom assertthat assert_that is.string is.dir
#' @importFrom stringr str_extract
#' @importFrom utils untar
extract_if_targz <- function(path) {
  assert_that(is.string(path))
  
  is_targz <- grepl("\\.tar\\.gz$", path)
  dirn <- ifelse(is_targz, stringr::str_extract(path, ".*(?=\\.tar\\.gz$)"), path)
  
  if (is_targz & !dir.exists(dirn)) {
    message(sprintf("Extracting %s to %s...", path, dirn))
    utils::untar(path, exdir=dirn)
  } else if (is_targz & dir.exists(dirn)) {
    message(sprintf("Skip extracting %s, because %s exists", path, dirn))
  } else { }  # !is_targz, so nothing to do.
  
  assert_that(is.dir(dirn))
  dirn
}
