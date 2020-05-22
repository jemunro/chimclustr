
#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
hap_mix_em <- function(read_allele_mat,
                       var_pos,
                       haps,
                       hap_prop,
                       ts_rate,
                       ts_max = 3,
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
            nrow(haps) == nrow(read_allele_mat),
            is_double(hap_prop),
            length(hap_prop) == ncol(haps),
            all(hap_prop > 0),
            is_scalar_double(ts_rate),
            ts_rate > 0 && ts_rate < 1,
            is_scalar_double(var_error_rate),
            is_scalar_double(read_error_rate),
            var_error_rate > 0 && var_error_rate < 1,
            read_error_rate > 0 && read_error_rate < 1,
            is_bool(fixed_hap_prop),
            is_scalar_double(epsilon),
            epsilon > 0,
            is.null(read_weight) || length(read_weight) == ncol(read_allele_mat))

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

  if (length(var_error_rate) < n_var) {
    var_error_rate <- rep(var_error_rate[1], n_var)
  }

  if (length(read_error_rate) < n_read) {
    read_error_rate <- rep(read_error_rate[1], n_read)
  }

  if (is.null(read_weight)) {
    read_weight <- rep(1, n_read)
  }
  read_weight_sum <- sum(read_weight)

  chims <- enum_chimeras(haps, var_pos = var_pos, ts_max = ts_max)
  n_chim <- nrow(chims$chime_space)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)

  read_chim_dist <- read_chimera_dist(read_allele_mat, chim_alleles)
  read_chim_ident <- n_var - read_chim_dist
  read_chim_match <- read_chimera_match2(read_allele_mat, chim_alleles)

  n_iter <- 0L
  LH <- -Inf
  search <- tibble(iter = integer(),
                   LH = double(),
                   ts_rate = double(),
                   var_error_rate = list(),
                   hap_prop = list(),
                   read_error_rate = list())

  while(n_iter < max_iter) {

    n_iter <- n_iter + 1L

    ## Expecation
    chim_lh <- chimera_lh(chims, hap_prop, ts_rate, var_width = var_width)

    re_mat <- matrix(read_error_rate, nrow = n_var,  ncol = n_read, byrow = TRUE)
    ve_mat <- matrix(var_error_rate, nrow = n_var, ncol = n_read)
    p_err <- log(re_mat + ve_mat - re_mat * ve_mat)
    inv_p_err <- log(1 - exp(p_err))

    read_chim_lh <-
      future_map(seq_len(n_var), function(vi) {
        map(seq_len(n_read), function(ri) {
          rep(inv_p_err[vi, ri], n_chim) %>%
            replace(which(!read_chim_match[[vi]][, ri]), p_err[vi, ri])
        }) %>% do.call('cbind', .)
      }) %>%
      purrr::reduce(`+`) +
      chim_lh$lh

    read_chim_posterior <-
      read_chim_lh %>%
      exp() %>%
      apply(2, function(x) x / sum(x)) %>%
      matrix(ncol = ncol(read_chim_lh))

    chim_posterior <-
      sweep(read_chim_posterior, 2, read_weight, '*') %>%
      rowSums() / read_weight_sum

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
                      var_error_rate = list(var_error_rate),
                      hap_prop = list(hap_prop),
                      read_error_rate = list(read_error_rate))

    ## need to calculate probailities under the model of
    # each haplotype or chimera
    # this can be pulled from chim_lh straigthforwardly
    # should include column here for is chimera

    if (epsilon > LH - LH_last) {
      break
    }

    ## Maximisation
    if (!fixed_hap_prop) {
      hap_prop <-
        select(chims$chime_space, starts_with('p')) %>%
        map_dbl(~ weighted.mean(., chim_posterior))
    }

    if (n_hap > 1) {
      ts_rate <- weighted.mean(chim_lh$en, chim_posterior) / var_width
    }

    var_err_posterior <- exp(log(ve_mat) - p_err)

    var_error_rate <-
      future_map_dbl(seq_len(n_var), function(i) {
        map_dbl(seq_len(n_read), function(j) {
          if_else(read_chim_match[[i]][, j], 0, var_err_posterior[i, j]) %>%
            weighted.mean(w = read_chim_posterior[, j])
        }) %>% weighted.mean(w = read_weight)
      }) %>%
      pmin(0.5)

    read_err_posterior <- exp(log(re_mat) - p_err)

    read_error_rate <-
      future_map_dbl(seq_len(n_read), function(j) {
        map_dbl(seq_len(n_var), function(i) {
          if_else(read_chim_match[[i]][, j], 0, read_err_posterior[i, j]) %>%
            weighted.mean(w = read_chim_posterior[, j])
        }) %>% mean()
      }) %>%
      pmin(0.5)
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

read_chimera_dist <- function(read_alleles,
                              chimera_alleles,
                              read_derep = NULL,
                              chimera_derep = NULL) {

  n_read <- ncol(read_alleles)
  n_chim <- ncol(chimera_alleles)
  n_var <- nrow(read_alleles)

  dist_mat <- matrix(0L, ncol = n_read, nrow = n_chim)

  for (i in seq_len(n_var)) {
    var_set <- union(unique(read_alleles[i, ]),
                     unique(chimera_alleles[i, ]))
    for (v in var_set) {
      read_set <- which(read_alleles[i,] == v)
      chim_set <- which(chimera_alleles[i,] == v)
      dist_mat[chim_set, read_set] <- dist_mat[chim_set, read_set] + 1L
    }
  }

  dist_mat <- n_var - dist_mat

  if (!is.null(read_derep)) {
    dist_mat_2 <- matrix(0L, ncol = sum(lengths(read_derep)), nrow = nrow(dist_mat))
    for (i in seq_along(read_derep)) {
      dist_mat_2[, read_derep[[i]]] <- dist_mat[, i]
    }
    dist_mat <- dist_mat_2
  }

  if (!is.null(chimera_derep)) {
    dist_mat_2 <- matrix(0L, ncol = sum(lengths(chimera_derep)), nrow = ncol(dist_mat))
    for (i in seq_along(chimera_derep)) {
      dist_mat_2[, chimera_derep[[i]]] <- dist_mat[i,]
    }
    dist_mat <- t(dist_mat_2)
  }

  return(dist_mat)
}

read_chimera_match <- function(read_alleles,
                               chimera_alleles) {

  n_read <- ncol(read_alleles)
  n_chim <- ncol(chimera_alleles)
  n_var <- nrow(read_alleles)

  map(seq_len(n_var), function(i) {
    var_set <- union(unique(read_alleles[i, ]),
                     unique(chimera_alleles[i, ]))
    map(var_set, function(v) {
      ri <- which(read_alleles[i,] == v) %>% unname()
      ci <- which(chimera_alleles[i,] == v) %>% unname()
      if (length(ri) > 0 && length(ci) > 0) {
        list(ri = ri, ci = ci)
      }
    }) %>% keep(~ !is.null(.))
  })
}

read_chimera_match_sub <- function(match_list,
                                   n_read,
                                   n_chim,
                                   x0, x1) {
  res_mat <- matrix(x0, ncol = n_read, nrow = n_chim)
  for (it in match_list) {
    res_mat[it$ci, it$ri] <- x1
  }
  return(res_mat)
}

read_chimera_match2 <- function(read_alleles,
                                chimera_alleles) {

  n_read <- ncol(read_alleles)
  n_chim <- ncol(chimera_alleles)
  n_var <- nrow(read_alleles)

  map(seq_len(n_var), function(i) {
    mat <- matrix(FALSE, nrow = n_chim, ncol = n_read)
    var_set <- union(unique(read_alleles[i, ]),
                     unique(chimera_alleles[i, ]))
    for (v in var_set) {
      rs <- which(read_alleles[i,] == v) %>% unname()
      cs <- which(chimera_alleles[i,] == v) %>% unname()
      if (length(rs) > 0 && length(cs) > 0) {
        mat[cs, rs] <- TRUE
      }
    }
    mat
  })
}
