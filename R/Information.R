Error <- function(message, ...) {
  stop('[x] ', message, ..., '\n', call. = FALSE)
}

Message <- function(message, ...) {
  message('[i] ', message, ..., '\n')
}

Warning <- function(message, ...) {
  warning('[!] ',message, ..., '\n')
}
