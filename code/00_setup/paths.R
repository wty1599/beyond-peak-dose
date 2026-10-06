release_config <- function() {
  p <- Sys.getenv('BEYOND_PATHS_CONFIG', unset='config/paths.yml')
  if (!file.exists(p)) stop('Copy paths.example.yml to paths.yml and supply private input locations.')
  jsonlite::fromJSON(p, simplifyVector=TRUE)
}
release_path <- function(key) {
  value <- release_config()[[key]]
  if (is.null(value) || !is.character(value) || length(value)!=1L ||
      !nzchar(value) || startsWith(value,'REPLACE_')) stop(paste('Configure',key))
  value
}
private_output <- function(stage) {
  root <- release_path('private_output')
  p <- file.path(root,stage)
  dir.create(p,recursive=TRUE,showWarnings=FALSE)
  p
}
