Error <- function(message, ...) {
  stop('[x] ', message, ..., '\n', call. = FALSE)
}

Message <- function(message, ...) {
  message('[i] ', message, ...)
}

Warning <- function(message, ...) {
  warning('[!] ',message, ...)
}
