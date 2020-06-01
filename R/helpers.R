
inv_order <- function(ord) {
  stopifnot(is.integer(ord),
            is.vector(ord))

  io <- integer(length(ord))
  for (i in seq_along(ord)) {
    io[ord[i]] <- i
  }
  return(io)
}

log_add <- function (x, y) {
  .max <- pmax(x, y)
  .min <- pmin(x, y)
  .max + log1p(exp(.min - .max))
}

#' @importFrom purrr reduce
log_sum <- function(x) {
  if (length(x) <= 1) {
    x
  } else {
    reduce(sort(x, TRUE), function (a, b) a + log1p(exp(b-a)))
  }
}

cosine <- function(x, y) {
  stopifnot(is.vector(x),
            is.vector(y),
            length(x) == length(y))
  return(c(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))))
}

inv_arr_ind <- function(row_idx, col_idx, nrow) {

  stopifnot(length(row_idx) == length(col_idx),
            max(row_idx) <= nrow,
            rlang::is_integerish(row_idx),
            rlang::is_integerish(col_idx),
            rlang::is_scalar_integerish(nrow))

  (col_idx - 1L) * nrow + row_idx
}

