
#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
#' @importFrom abind abind
hap_mix_em2 <- function(read_allele_mat,
                        var_pos,
                        haps,
                        hap_prop,
                        ts_rate,
                        ts_max = 2,
                        var_error_rate = 1e-2,
                        read_error_rate = 1e-2,
                        max_var_error_rate = 0.5,
                        max_iter = 20,
                        fixed_hap_prop = FALSE,
                        epsilon = 1e-3,
                        read_weight = NULL) {

  # bring back read_count to speed up search

  stopifnot(is.matrix(read_allele_mat),
            is_integer(read_allele_mat),
            is_integer(var_pos),
            length(var_pos) == nrow(read_allele_mat),
            is.matrix(haps) && ncol(haps) >= 1L,
            !any(is.na(haps)),
            nrow(haps) == nrow(read_allele_mat),
            is_double(hap_prop),
            length(hap_prop) == ncol(haps),
            all(hap_prop > 0),
            is_scalar_double(ts_rate),
            ts_rate > 0 && ts_rate < 1,
            is_double(var_error_rate),
            length(var_error_rate) == 1 | length(var_error_rate) == nrow(read_allele_mat),
            all(var_error_rate >= 0 && var_error_rate <= 1),
            is_double(var_miss_rate),
            length(var_miss_rate) == 1 | length(var_miss_rate) == nrow(read_allele_mat),
            all(var_miss_rate >= 0 && var_miss_rate <= 1),
            is_double(var_miss_rate),
            length(read_error_rate) == 1 | length(read_error_rate) == ncol(read_allele_mat),
            all(read_error_rate >= 0 && read_error_rate <= 1),
            is_double(read_miss_rate),
            length(read_miss_rate) == 1 | length(read_miss_rate) == ncol(read_allele_mat),
            all(read_miss_rate >= 0 && read_miss_rate <= 1),
            is_bool(fixed_hap_prop),
            is_scalar_double(epsilon),
            epsilon > 0,
            is.null(read_weight) || length(read_weight) == ncol(read_allele_mat))

  read_allele_mat <- unname(read_allele_mat)

  hap_prop <- hap_prop / sum(hap_prop)
  n_hap <- ncol(haps)
  if (n_hap == 1) {
    # chimera's are irrelevant for single haplotype
    ts_rate <- 0
    ts_max <- 0
  }
  n_var <- length(var_pos)
  n_read <- ncol(read_allele_mat)
  var_width <- last(var_pos) - first(var_pos)

  if (length(var_error_rate) != n_var) {
    var_error_rate <- rep(var_error_rate[1], n_var)
  }

  if (length(read_error_rate) != n_read) {
    read_error_rate <- rep(read_error_rate[1], n_read)
  }

  read_error_rate <- rep(0.01, n_read)
  read_miss_rate <- rep(0.01, n_read)
  var_miss_rate <- rep(0.01, n_var)
  var_error_rate <- rep(0.01, n_var)

  if (is.null(read_weight)) {
    read_weight <- rep(1, n_read)
  }
  read_weight_sum <- sum(read_weight)

  chims <- enum_chimeras(haps, var_pos = var_pos, ts_max = ts_max)
  n_chim <- nrow(chims$chime_space)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)

  lvls <- c('match', 'mismatch', 'missing')
  chim_read_var_state <- array(integer(), dim = c(n_chim, n_read, n_var))
  for (v in seq_len(n_var)) {
    for (r in seq_len(n_read)) {
      chim_read_var_state[, r, v] <-
        if_else(chim_alleles[v, ] == read_allele_mat[v, r],
                1L, 2L, missing = 3L)
    }
  }
  var_read_missing <- map(seq_len(n_var), function(vi) {
    list(yes = which(is.na(read_allele_mat[vi, ])),
         no = which(!is.na(read_allele_mat[vi, ])))
  })
  read_var_missing <- map(seq_len(n_read), function(ri) {
    list(yes = which(is.na(read_allele_mat[, ri])),
         no = which(!is.na(read_allele_mat[, ri ])))
  })
  # read_chim_dist <- read_chimera_dist(read_allele_mat, chim_alleles)
  # read_chim_ident <- n_var - read_chim_dist
  # read_chim_match <- read_chimera_match2(read_allele_mat, chim_alleles)

  n_iter <- 0L
  LH <- -Inf
  search <- tibble(iter = integer(),
                   LH = double(),
                   ts_rate = double(),
                   hap_prop = list(),
                   var_error_rate = list(),
                   var_miss_rate = list(),
                   read_error_rate = list(),
                   read_miss_rate = list())

  while(n_iter < max_iter) {

    n_iter <- n_iter + 1L

    ##### Expecation #####
    chim_lh <- chimera_lh(chims, hap_prop, ts_rate, var_width = var_width)

    # LH of missing/mismatch/match at all read vars
    missing_0 <- # p(A or B) = p(A) + p(B) - p(A) * p(B)
      matrix(var_miss_rate, nrow = n_var, ncol = n_read) %>%
      (function(x) sweep(x, 2, read_miss_rate, `+`, check.margin = FALSE) -
         sweep(x, 2, read_miss_rate, `*`, check.margin = FALSE))

    mismatch_0 <-
      matrix(var_error_rate, nrow = n_var, ncol = n_read) %>%
      (function(x) sweep(x, 2, read_error_rate, `+`, check.margin = FALSE) -
         sweep(x, 2, read_error_rate, `*`, check.margin = FALSE))

    missing <- log(missing_0)
    mismatch <- log(1 - missing_0) + log(mismatch_0)
    match <- log(1 - missing_0) + log(1 - mismatch_0)

    read_chim_lh <-
      map(seq_len(n_var), function(vi) {
        map(seq_len(n_read), function(ri) {
          c(match[vi, ri], mismatch[vi, ri], missing[vi, ri])[chim_read_var_state[, ri, vi]]
        }) %>% do.call('cbind', .)
      }) %>%
      reduce(`+`) + chim_lh$lh

    read_chim_posterior <-
      read_chim_lh %>%
      exp() %>%
      apply(2, function(x) x / sum(x)) %>%
      matrix(ncol = ncol(read_chim_lh))

    chim_posterior <-
      sweep(read_chim_posterior, 2, read_weight, '*') %>%
      rowSums() / read_weight_sum

    var_err_posterior <- exp(sweep(-mismatch, 1, log(var_error_rate),  `+`))
    var_miss_posterior <- exp(sweep(-missing, 1, log(var_miss_rate),  `+`))
    read_err_posterior <- exp(sweep(-mismatch, 2, log(read_error_rate),  `+`))
    read_miss_posterior <- exp(sweep(-missing, 2, log(read_miss_rate),  `+`))

    LH_last <- LH

    LH <-
      read_chim_lh %>%
      exp() %>%
      colSums() %>%
      log() %>%
      (function(x) sum(x * read_weight))

    # record results
    search <- add_row(search,
                      iter = n_iter,
                      LH = LH,
                      ts_rate = ts_rate,
                      hap_prop = list(hap_prop),
                      var_error_rate = list(var_error_rate),
                      var_miss_rate = list(var_miss_rate),
                      read_error_rate = list(read_error_rate),
                      read_miss_rate = list(read_miss_rate))

    ## need to calculate probailities under the model of
    # each haplotype or chimera
    # this can be pulled from chim_lh straigthforwardly
    # should include column here for is chimera

    if (epsilon > LH - LH_last) {
      break
    }

    ##### Maximisation #####

    if (!fixed_hap_prop) {
      hap_prop <-
        select(chims$chime_space, starts_with('p')) %>%
        map_dbl(~ weighted.mean(., chim_posterior))
    }

    if (n_hap > 1) {
      ts_rate <- weighted.mean(chim_lh$en, chim_posterior) / var_width
    }

    var_miss_rate <-
      map_dbl(seq_len(n_var), function(vi) {
        weighted.mean(
          replace(numeric(n_read), var_read_missing[[vi]]$yes, var_miss_posterior[vi, var_read_missing[[vi]]$yes]),
          read_weight)
      })

    read_miss_rate <-
      map_dbl(seq_len(n_read), function(ri) {
       replace(numeric(n_var), read_var_missing[[ri]]$yes,  read_miss_posterior[read_var_missing[[ri]]$yes, ri]) %>%
          mean()
      })

    var_error_rate <-
      map_dbl(seq_len(n_var), function(vi) {
        map_dbl(var_read_missing[[vi]]$no, function(ri) {
          weighted.mean(
            c(0, var_err_posterior[vi, ri])[chim_read_var_state[, ri, vi]],
            read_chim_posterior[, ri])
        }) %>% weighted.mean(w = read_weight[var_read_missing[[vi]]$no])
      })

    read_error_rate <-
      map_dbl(seq_len(n_read), function(ri) {
        map_dbl(read_var_missing[[ri]]$no, function(vi) {
          weighted.mean(
            c(0, read_err_posterior[vi, ri])[chim_read_var_state[, ri, vi]],
            read_chim_posterior[, ri])
        }) %>% mean()
      })
  }

  read_hap_lh <-
    tibble(read_id = seq_len(n_read)) %>%
    bind_cols(
      select(chims$chime_space, starts_with('p')) %>%
        mutate(pchimera = if_else(rowSums(. > 0 & . < 1) > 0, 1, 0)) %>%
        map2_dfc(str_remove(names(.), 'p'), ., function(n, p) {
          tibble(!!str_c('raw_llh_', n) := log(colSums(exp(read_chim_lh[which(p == 1), , drop = FALSE]))),
                 !!str_c('post_lh_', n) := colSums(read_chim_posterior[which(p == 1), , drop = FALSE]))
        })
    ) %>%
    pivot_longer(-read_id,
                 names_to = c('.value', "haplotype"),
                 names_pattern = '(.+)_([^_]+)')

  read_hap_lh_smry <-
    read_hap_lh %>%
    arrange(read_id, desc(post_lh)) %>%
    group_by(read_id) %>%
    slice(1) %>%
    ungroup()

  return(list(search = search,
              read_hap_lh = read_hap_lh,
              read_hap_lh_smry = read_hap_lh_smry))
}

miss_em <- function(x, max_iter = 100L, epsilon = 1e-6) {
  # note that this missing rate is independant of other rates
  # therefore it can be estimated separately
  # can also be used for independant row/col filtering

  stopifnot(is.matrix(x), length(x) > 1)

  if (!any(is.na(x))) {
    # easy case
    return(list(row_rate = rep(0, nrow(x)),
                col_rate = rep(0, ncol(x)),
                n_iter = 0L,
                converged = NA,
                LH = 0))
  }

  n_row <- nrow(x)
  n_col <- ncol(x)
  is_na_x <- is.na(c(x))
  row_col_missing <- map(seq_len(n_row), function(i) which(is.na(x[i, ])))
  col_row_missing <- map(seq_len(n_col), function(j) which(is.na(x[, j])))
  row_rate <- lengths(row_col_missing) / n_col
  col_rate <- lengths(col_row_missing) / n_row
  iter <- 0L

  search <- tibble(iter = integer(),
                   LH = numeric(),
                   row_rate = list(),
                   col_rate = list())

  while(iter < max_iter) {
    iter <- iter + 1

    ## expectation
    miss_lh <-
      matrix(row_rate, nrow = n_row, ncol = n_col) %>%
      { sweep(., 2, col_rate, `+`, check.margin = FALSE) - sweep(., 2, col_rate, `*`, check.margin = FALSE) }

    non_miss_lh <- 1 - miss_lh

    LH <- sum(log(c(miss_lh[is_na_x], (1 - miss_lh)[!is_na_x])))

    search <- add_row(search,
                      iter = iter,
                      LH = LH,
                      row_rate = list(row_rate),
                      col_rate = list(col_rate))

    row_post <- sweep(1/miss_lh, 1, row_rate, `*`)
    col_post <- sweep(1/miss_lh, 2, col_rate, `*`)

    ## maximisation
    row_rate <-
      map_dbl(seq_len(n_row), function(i) {
        mean(replace(numeric(n_col), row_col_missing[[i]], row_post[i, row_col_missing[[i]]]))
      })

    col_rate <-
      map_dbl(seq_len(n_col), function(j) {
        mean(replace(numeric(n_row), col_row_missing[[j]], col_post[col_row_missing[[j]], j]))
      })
  }

}

estimate_marginal_rates <- function(read_allele_mat) {


  dissim <-
    t(read_allele_mat) %>%
    as.data.frame() %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')

  clustering <- unname(fpc::pamk(dissim)$pamobject$clustering)

}
