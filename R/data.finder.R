getTarBallPath <- function() {
  system.file("extdata", "filtered_gene_bc_matrices", "GRCh38", "data.tar.gz", package = "tmetutorial", mustWork = TRUE)
}

extract_data <- function() {
  tar_name = tmetutorial::getTarBallPath()
  tmp_dir = tempdir()
  untar(tar_name, exdir = tmp_dir)
  return(tmp_dir)
}

#' @export
getSingleCellDataPath <- function() {
  return(extract_data())
  #system.file("extdata", "filtered_gene_bc_matrices", "GRCh38", package = "tmetutorial", mustWork = TRUE)
}

#' @export
open_tutorial <-function() {
  fname = system.file("extdata", "tme_tutorial.R", package = "tmetutorial", mustWork = TRUE)
  file.edit(fname)
}
