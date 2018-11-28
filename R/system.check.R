check_packages <- function() {
  version = paste(R.version$major, R.version$minor, sep=".")
  version = substr(version, 1, 3)
  pkgs = installed.packages()
  SEURAT = "Seurat"
  if(!(SEURAT %in% rownames(pkgs))) {
    stop("Seurat is not installed which is needed for the tutorial.")
  }
  # Check for packages not built under this version of R
  details = pkgs[SEURAT,]
  versions = sapply(pkgs[,"Built"], function(s) substr(s, 1, 3))
  v_match = versions != version
  if(any(v_match)) {
    # Ensure mismatched package is Seurat dependency
    requirements = paste(details["Imports"], details["Depends"], details["Suggests"], sep="")
    bad_pkgs = names(v_match[v_match])
    is_dep = sapply(bad_pkgs, function(s) grepl(s, requirements))
    if(any(is_dep)) {
      bad = paste(bad_pkgs[is_dep], collapse = "\n")
      msg=paste("You have packages built under different versions of R than ",
                    "you are currently running.  This can lead to unexpected ",
                    "behavior and possible crashes.  Please reinstall/update the following packages ",
                    "before continuing:\n", bad, sep="", collapse = " ")
      stop(msg)
    }
  }
}
