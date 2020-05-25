
#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
#' @importFrom abind abind
hap_mix_em <- function(allele_mat,
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

  if (is.null(read_weight)) {
    read_weight <- rep(1, n_read)
  }
  read_weight_sum <- sum(read_weight)

  error_est <- est_err_rates(allele_mat, haps, read_weight)
  read_error_rate <- error_est$read_error_rate
  var_error_rate <- error_est$var_error_rate
  fixed_var_error_rate <- var_error_rate[which_fixed]
  unfixed_var_error_rate <- var_error_rate[which_unfixed]

  chims <- enum_chimeras(haps[which_unfixed, , drop = FALSE], var_pos = unfixed_pos, ts_max = ts_max)
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

    if (n_hap > 1) {
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

    } else {
      read_chim_lh <- matrix(0, nrow = n_chim, ncol = n_read)
      read_chim_posterior <- matrix(1, nrow = n_chim, ncol = n_read)
      chim_posterior <- 1
    }

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

#' @importFrom digest digest
enum_chimeras <- function(haps, var_pos, ts_max) {

  n_hap <- ncol(haps)
  n_var <- length(var_pos)

  if (n_hap == 1L) {
    hap_space <- tibble(ns = 0L,
                        haps = list(1L),
                        hsid = 1L)
    ts_space <- tibble(tsid = 1L,
                       ns = 0L,
                       vars = list(integer()),
                       width_non_zero = list(integer(0)),
                       width_zero = list(1L),
                       state_span = list(n_var - 1L))
    chime_space <- tibble(chid = 1L,
                         ids = list(tibble(hsid = 1L, tsid = 1L)),
                         ns = 0L,
                         alleles = list(c(haps)),
                         p1 = 1,
                         is_chimera = FALSE)
  } else {

    var_set <- seq_along(var_pos)
    var_dist <- as.integer(as.matrix(dist(var_pos))) %>%
      { matrix(., ncol = n_var)}

    hap_space <-
      seq_len(ts_max + 1) %>%
      map_df(function(nh) {
        list(seq_len(n_hap)) %>%
          rep(nh) %>%
          setNames(., seq_along(.))%>%
          do.call('expand_grid', .) %>%
          mutate(id = seq_len(n())) %>%
          gather(-id, key = 'index', value = 'haps') %>%
          arrange(id, index) %>%
          select(-index) %>%
          mutate(ns = nh - 1) %>%
          chop(haps) %>%
          select(-id)
      }) %>%
      mutate(hsid = seq_len(n()))

    ts_space <-
      seq_len(ts_max) %>%
      map_df(function(ns) {
        vars <- t(combn(n_var - 1L, ns))
        var_index <- seq_len(nrow(vars))

        nrc <-  2 * (ns + 1)
        ranges <- matrix(0L, nrow = nrow(vars), ncol = nrc)
        ranges[, 1] <- 1L
        ranges[,  nrc] <- n_var
        ranges[, 2 * seq_len(ns)] <- vars
        ranges[, 2 * seq_len(ns) + 1] <- vars + 1L
        non_zero_at <- seq.int(2, nrc-1, 2)
        zero_at <- seq.int(1, nrc-1, 2)

        spans <- ranges[, 1 + seq_len(nrc - 1)] - ranges[, seq_len(nrc - 1)]
        state_spans <- matrix(0L, nrow = nrow(vars), ncol = ns + 1)
        state_spans[, seq_len(ns)] <- spans[, non_zero_at] + spans[, zero_at[-(ns+1)]]
        state_spans[, ns + 1] <- spans[, last(zero_at)]

        widths <-
          var_dist[inv_arr_ind(ranges[, seq_len(nrc - 1)],
                               ranges[, 1 + seq_len(nrc - 1)],
                               nrow(var_dist))] %>%
          matrix(ncol = nrc - 1)

        tibble(ns = ns,
               vars = map(var_index, ~ vars[., ]),
               width_non_zero = map(var_index, ~ widths[., non_zero_at]),
               width_zero = map(var_index, ~ { widths[., zero_at] %>% keep( ~ . > 0) }),
               state_span = map(var_index, ~ state_spans[., ]))

      }) %>%
      mutate(tsid = seq_len(n()) + 1L) %>%
      { bind_rows(tibble(tsid = 1L,
                         ns = 0L,
                         vars = list(integer()),
                         width_non_zero = list(integer()),
                         width_zero = list(var_dist[1, n_var]),
                         state_span = list(n_var - 1L)),
                  .) } %>%
      select(tsid, ns, everything())

    ts_state <-
      map(ts_space$state_span, function(ss) rep(seq_along(ss), ss)) %>%
      do.call('cbind', .) %>%
      { rbind(1L, .) }

    chime_space <-
      hap_space %>%
      full_join(select(ts_space, tsid, ns), by = 'ns') %>%
      mutate(hap_hash = map2_chr(haps, tsid, function(h, i) digest(h[ts_state[, i]]))) %>%
      group_by(hap_hash) %>%
      (function(x) {
        full_join(
          select(x, -haps) %>% nest(ids = c(hsid, tsid)) %>% ungroup(),
          slice(x, 1) %>% mutate(hap_state = map2(haps, tsid, function(h, i) h[ts_state[, i]])) %>%
            ungroup() %>% select(hap_hash, hap_state),
          by = 'hap_hash')
      }) %>%
      mutate(chid = seq_len(n())) %>%
      mutate(alleles = map(hap_state, function(hs) {
        haps[inv_arr_ind(seq_along(hs), hs, length(hs))]
      })) %>%
      (function(x) {
        seq_len(n_hap) %>%
          setNames(., str_c('p', seq_along(.))) %>%
          map_dfc(function(h) {
            map_int(x$hap_state, ~ sum(. == h)) / n_var
          }) %>%
          mutate(is_chimera = rowSums(. > 0 & . < 1) > 1) %>%
          bind_cols(x, .)
      }) %>%
      select(chid, ids, ns, alleles, starts_with('p'), is_chimera)
  }

  return(list(hap_space = hap_space, ts_space = ts_space, chime_space = chime_space))
}

chimera_lh <- function(chimeras, hap_prop, ts_rate, var_width, rescale_lh = TRUE) {

  hap_lh <-
    chimeras$hap_space %>%
    mutate(hslh = map_dbl(haps, ~ sum(log(hap_prop[.])))) %>%
    select(hsid, hslh)

  ts_lh <-
    bind_rows(chimeras$ts_space %>%
                select(tsid, width = width_zero) %>%
                unnest(width) %>%
                mutate(lh = dbinom(0, width, ts_rate, log = TRUE),
                       en = 0),
              chimeras$ts_space %>%
                select(tsid, width = width_non_zero) %>%
                unnest(width) %>%
                mutate(lh = pbinom(0, width, ts_rate, lower.tail = FALSE, log.p = TRUE),
                       en = exp(log(ts_rate) + log(width) - lh))) %>%
    group_by(tsid) %>%
    summarise(tslh = sum(lh, na.rm = T), en = sum(en))


  if (rescale_lh) {
    # ts_rate and ts_lh will be under-estimated if we exclude events with rate > ts_max
    # we can rescale en to account for this simplification
    ts_lh$tslh %<>% exp() %>% { . / sum(.) } %>% log()
    en_sub <- with(ts_lh, sum(exp(tslh) * en))
    ts_lh$en %<>% { . * ( ts_rate * var_width / en_sub) }
  }

  chime_lh <-
    chimeras$chime_space %>%
    select(chid, ids) %>%
    unnest(ids) %>%
    left_join(hap_lh, 'hsid') %>%
    left_join(ts_lh, 'tsid') %>%
    group_by(chid) %>%
    summarise(lh = log_sum(hslh + tslh),
              en = weighted.mean(en, exp(tslh)))

  return(chime_lh)
}

est_err_rates <- function(allele_mat, haps, read_weight) {

  stopifnot(is.matrix(allele_mat),
            is.integer(allele_mat),
            is.matrix(haps),
            is.integer(haps),
            nrow(haps) == nrow(allele_mat))

  n_hap <- ncol(haps)
  n_read <- ncol(allele_mat)
  n_var <- nrow(allele_mat)

  clust <-
    t(cbind(haps, allele_mat)) %>%
    as_tibble() %>%
    mutate_all(~ replace_na(., -1L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower') %>%
    as.matrix() %>%
    (function(x) {
      x[-seq_len(n_hap), seq_len(n_hap), drop = FALSE] %>%
        unname() %>%
        apply(1, which.max)
    })

  is_err <-
    map(seq_len(n_read), function(i) allele_mat[, i] != haps[, clust[i]]) %>%
    do.call('cbind', .)

  is_err_count <- sweep(is_err, 2, read_weight, `*`)
  isnt_na_count <- sweep(!is.na(is_err), 2, read_weight, `*`)

  global_rate <- sum(is_err_count, na.rm = T) / sum(isnt_na_count)
  # simple correction for over estimation of missingness when marginalising
  adj <- (1-sqrt(1-global_rate)) / global_rate
  row_rate <- adj * rowSums(is_err_count, na.rm = T) / rowSums(isnt_na_count)
  col_rate <- adj * colSums(is_err_count, na.rm = T) / colSums(isnt_na_count)

  return(list(var_error_rate = row_rate, read_error_rate = col_rate))
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
