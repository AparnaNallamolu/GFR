Error <- function(message, ...) {
  stop('[x] ', message, ..., '\n', call. = FALSE)
}

Message <- function(message, show = TRUE, ...) {
  if(show){ message('[i] ', message, ...)}
}

Warning <- function(message, ...) {
  warning('[!] ',message, ...)
}
