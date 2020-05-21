
#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
#' @importFrom abind abind
hap_mix_em2 <- function(allele_mat,
                        var_pos,
                        haps,
                        hap_prop,
                        ts_rate = 1e-4,
                        ts_max = 2,
                        max_var_error_rate = 0.5,
                        max_iter = 20,
                        fixed_hap_prop = FALSE,
                        epsilon = 1e-3,
                        read_weight = NULL) {

  # bring back read_count to speed up search

  stopifnot(is.matrix(allele_mat),
            is_integer(allele_mat),
            is_integer(var_pos),
            length(var_pos) == nrow(allele_mat),
            is.matrix(haps) && ncol(haps) >= 1L,
            !any(is.na(haps)),
            nrow(haps) == nrow(allele_mat),
            is_double(hap_prop),
            length(hap_prop) == ncol(haps),
            all(hap_prop > 0),
            is_scalar_double(ts_rate),
            ts_rate > 0 && ts_rate < 1,
            is_bool(fixed_hap_prop),
            is_scalar_double(epsilon),
            epsilon > 0,
            is.null(read_weight) || length(read_weight) == ncol(allele_mat))

  allele_mat <- unname(allele_mat)

  hap_prop <- hap_prop / sum(hap_prop)
  n_hap <- ncol(haps)
  if (n_hap == 1) {
    # chimera's are irrelevant for single haplotype
    ts_rate <- 0
    ts_max <- 0
  }
  n_var <- length(var_pos)
  n_read <- ncol(allele_mat)
  is_fixed_var <- map_lgl(seq_len(n_var), ~ all(haps[., ] == haps[., 1]))
  which_fixed <- which(is_fixed_var)
  which_unfixed <- which(!is_fixed_var)
  n_unfixed_var <- length(which_unfixed)
  n_fixed_var <- length(which_fixed)
  unfixed_pos <- var_pos[which_unfixed]
  unfixed_width <- last(unfixed_pos) - first(unfixed_pos)
  consensus_alleles <- haps[which_fixed, 1]
  fixed_alleles <- allele_mat[which_fixed, ]
  unfixed_alleles <- allele_mat[which_unfixed, ]
  error_est <- est_err_rates(allele_mat, haps)
  read_error_rate <- error_est$read_error_rate
  var_error_rate <- error_est$var_error_rate
  fixed_var_error_rate <- var_error_rate[which_fixed]
  unfixed_var_error_rate <- var_error_rate[which_unfixed]

  if (is.null(read_weight)) {
    read_weight <- rep(1, n_read)
  }
  read_weight_sum <- sum(read_weight)

  chims <- enum_chimeras(haps[which_unfixed, , drop = FALSE], var_pos = unfixed_pos, ts_max = 2)
  n_chim <- nrow(chims$chime_space)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)

  # state coding - 1:'match', 2:'mismatch', 3:'missing'
  chim_read_var_state <-
    map(seq_len(n_unfixed_var), function(v) {
      map(seq_len(n_read), function(r) {
        if_else(chim_alleles[v, ] == unfixed_alleles[v, r], 1L, 2L, missing = 3L)
      }) %>% do.call('cbind', .)
    }) %>% abind::abind(along = 3)

  fixed_read_var_state <-
    map(seq_len(n_fixed_var), function(v) {
      if_else(consensus_alleles[v] == fixed_alleles[v, ], 1L, 2L, missing = 3L)
    }) %>%
    do.call('rbind', .)

  n_iter <- 0L
  LH <- -Inf
  search <- tibble(iter = integer(),
                   LH = double(),
                   ts_rate = double(),
                   hap_prop = list(),
                   read_error_rate = list(),
                   var_error_rate = list())

  while(n_iter < max_iter) {
    n_iter <- n_iter + 1L
    message(n_iter)
    ##### Expecation #####
    chim_lh <- chimera_lh(chims, hap_prop, ts_rate, var_width = unfixed_width)

    fixed_mismatch_e <-
      matrix(fixed_var_error_rate, nrow = n_fixed_var, ncol = n_read) %>%
      (function(x) sweep(x, 2, read_error_rate, `+`, check.margin = FALSE) -
         sweep(x, 2, read_error_rate, `*`, check.margin = FALSE))

    fixed_match_lh <- log(1-fixed_mismatch_e)
    fixed_mismatch_lh <- log(fixed_mismatch_e)

    read_fixed_lh <-
      map_dbl(seq_len(n_read), function(ri) {
        map_dbl(seq_len(n_fixed_var), function(vi) {
          c(fixed_match_lh[vi, ri], fixed_mismatch_lh[vi, ri], 0)[fixed_read_var_state[vi, ri]]
        }) %>% sum()
      })

    unfixed_mismatch_e <-
      matrix(unfixed_var_error_rate, nrow = n_unfixed_var, ncol = n_read) %>%
      (function(x) sweep(x, 2, read_error_rate, `+`, check.margin = FALSE) -
         sweep(x, 2, read_error_rate, `*`, check.margin = FALSE))

    unfixed_match_lh <- log(1-unfixed_mismatch_e)
    unfixed_mismatch_lh <- log(unfixed_mismatch_e)

    read_chim_lh <-
      map(seq_len(n_unfixed_var), function(vi) {
        map(seq_len(n_read), function(ri) {
          c(unfixed_match_lh[vi, ri], unfixed_mismatch_lh[vi, ri], 0)[chim_read_var_state[, ri, vi]]
        }) %>% do.call('cbind', .)
      }) %>%
      reduce(`+`) + chim_lh$lh

    # posterior likelihood of each chimeric state for each read
    read_chim_posterior <-
      read_chim_lh %>%
      exp() %>%
      apply(2, function(x) x / sum(x)) %>%
      matrix(ncol = ncol(read_chim_lh))

    # posterior likelihood of each chimeric state across all reads
    chim_posterior <-
      sweep(read_chim_posterior, 2, read_weight, '*') %>%
      rowSums() / read_weight_sum

    # posterior error likelihoods given error
    fixed_read_err_post <- exp(sweep(-fixed_mismatch_lh, 2, log(read_error_rate), `+`, check.margin = FALSE))
    fixed_var_err_post <- exp(sweep(-fixed_mismatch_lh, 1, log(fixed_var_error_rate), `+`, check.margin = FALSE))
    unfixed_read_err_post <- exp(sweep(-unfixed_mismatch_lh, 2, log(read_error_rate), `+`, check.margin = FALSE))
    unfixed_var_err_post <- exp(sweep(-unfixed_mismatch_lh, 1, log(unfixed_var_error_rate), `+`, check.margin = FALSE))

    LH_last <- LH

    LH <-
      read_chim_lh %>%
      exp() %>% colSums() %>% log() %>%
      (function(x) x + read_fixed_lh ) %>%
      (function(x) sum(x * read_weight))

    # record results
    search <- add_row(search,
                      iter = n_iter,
                      LH = LH,
                      ts_rate = ts_rate,
                      hap_prop = list(hap_prop),
                      read_error_rate = list(read_error_rate),
                      var_error_rate = list(var_error_rate))

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
      ts_rate <- weighted.mean(chim_lh$en, chim_posterior) / unfixed_width
    }

    fixed_var_error_rate <-
      map_dbl(seq_len(n_fixed_var), function(vi) {
        map_dbl(seq_len(n_read), function(ri) {
          c(0, fixed_var_err_post[vi, ri], NA_real_)[fixed_read_var_state[vi, ri]]
        }) %>% weighted.mean(read_weight, na.rm = TRUE)
      }) %>%
      pmin(0.5)

    unfixed_var_error_rate <-
      map_dbl(seq_len(n_unfixed_var), function(vi) {
        map_dbl(seq_len(n_read), function(ri) {
          weighted.mean(
            c(0, unfixed_var_err_post[vi, ri], NA_real_)[chim_read_var_state[, ri, vi]],
            read_chim_posterior[, ri])
        }) %>% weighted.mean(read_weight, na.rm = TRUE)
      }) %>%
      pmin(0.5)

    var_error_rate[which_fixed] <- fixed_var_error_rate
    var_error_rate[which_unfixed] <- unfixed_var_error_rate

    read_error_rate <-
      map_dbl(seq_len(n_read), function(ri) {
        c(
          map_dbl(seq_len(n_fixed_var), function(vi) {
            c(0, fixed_read_err_post[vi, ri], NA_real_)[fixed_read_var_state[vi, ri]]
          }),
          map_dbl(seq_len(n_unfixed_var), function(vi) {
            weighted.mean(
              c(0, unfixed_read_err_post[vi, ri], NA_real_)[chim_read_var_state[, ri, vi]],
              read_chim_posterior[, ri])
          })) %>%
          mean(na.rm = TRUE)
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

est_err_rates <- function(allele_mat, haps) {
  # todo: allow read weights here for working with dereplicated read sets

    stopifnot(is.matrix(allele_mat),
            is.integer(allele_mat),
            is.matrix(haps),
            is.integer(haps),
            nrow(haps) == nrow(allele_mat))

  n_hap <- ncol(haps)
  n_read <- ncol(allele_mat)

  clust <-
    t(cbind(haps, allele_mat)) %>%
    as_tibble() %>%
    mutate_all(~ replace_na(., -1L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower') %>%
    as.matrix() %>%
    (function(x) {
      x[-seq_len(n_hap), seq_len(n_hap)] %>%
        unname() %>%
        apply(1, which.max)
    })

  is_err <-
    map(seq_len(n_read), function(i) {
      allele_mat[, i] != haps[, clust[i]]
    }) %>%
    do.call('cbind', .)

  n_row <- nrow(is_err)
  n_col <- ncol(is_err)
  global_rate <- sum(is_err, na.rm = T) / sum(!is.na(is_err))
  # correction for over estimation of missingness when marginalising
  adj <- (1-sqrt(1-global_rate)) / global_rate
  row_col_which <- map(seq_len(n_row), function(i) which(is_err[i, ]))
  row_col_size <- map_int(seq_len(n_row), function(i) sum(!is.na(is_err[i, ])))
  col_row_which <- map(seq_len(n_col), function(j) which(is_err[, j]))
  col_row_size <- map_int(seq_len(n_col), function(j) sum(!is.na(is_err[, j])))
  row_rate <- adj*(lengths(row_col_which) / row_col_size)
  col_rate <- adj*(lengths(col_row_which) / col_row_size)

  return(list(var_error_rate = row_rate,
              read_error_rate = col_rate))
}

marginal_rate_em <- function(x, max_iter = 100L, epsilon = 1e-3) {

  # EM-framework
  # observed variables: sites true
  # latent variables: reason true (v/s/vs)
  # params: row-wise and col-wise rates

  stopifnot(is.matrix(x),
            is.logical(x),
            all(!is.na(x)),
            length(x) > 1)

  if (!any(x)) {
    # easy case
    return(list(row_rate = rep(0, nrow(x)),
                col_rate = rep(0, ncol(x)),
                n_iter = 0L,
                converged = NA,
                LH = 0))
  }
  x <- unname(x)
  n_row <- nrow(x)
  n_col <- ncol(x)
  global_rate <- sum(x) / length(x)
  # correction for over estimation of missingness when marginalising
  adj <- (1-sqrt(1-global_rate)) / global_rate
  row_col_which <- map(seq_len(n_row), function(i) which(x[i, ]))
  col_row_which <- map(seq_len(n_col), function(j) which(x[, j]))
  row_rate <- adj*(lengths(row_col_which) / n_col)
  col_rate <- adj*(lengths(col_row_which) / n_row)
  iter <- 0L
  LH_last <- -Inf

  search <- tibble(iter = integer(),
                   LH = numeric(),
                   row_rate = list(),
                   col_rate = list())

  while(iter < max_iter) {
    iter <- iter + 1

    ## expectation
    x_lh <-
      matrix(row_rate, nrow = n_row, ncol = n_col) %>%
      { sweep(., 2, col_rate, `+`, check.margin = FALSE) - sweep(., 2, col_rate, `*`, check.margin = FALSE) }

    LH <- sum(log(c(x_lh[x], (1 - x_lh)[!x])))

    search <- add_row(search,
                      iter = iter,
                      LH = LH,
                      row_rate = list(row_rate),
                      col_rate = list(col_rate))

    if (epsilon > (LH - LH_last)) {
      break
    }

    LH_last <- LH

    ## maximisation
    row_post <- sweep(1/x_lh, 1, row_rate, `*`)
    col_post <- sweep(1/x_lh, 2, col_rate, `*`)

    row_rate <-
      map_dbl(seq_len(n_row), function(i) {
        mean(replace(numeric(n_col), row_col_which[[i]], row_post[i, row_col_which[[i]]]))
      })

    col_rate <-
      map_dbl(seq_len(n_col), function(j) {
        mean(replace(numeric(n_row), col_row_which[[j]], col_post[col_row_which[[j]], j]))
      })
  }

  return(list(row_rate = row_rate, col_rate = col_rate, search = search))
}
