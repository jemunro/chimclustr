#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
hap_mix_em2 <- function(read_allele_mat,
                       var_pos,
                       haps,
                       hap_prop,
                       ts_rate,
                       ts_max = 3,
                       var_error_rate,
                       flat_var_error_rate = FALSE,
                       max_var_error_rate = 0.5,
                       max_iter = 50,
                       fixed_hap_prop = FALSE,
                       epsilon = 1e-3) {
  require(tidyverse)

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
            var_error_rate > 0 && var_error_rate < 1,
            is_bool(fixed_hap_prop),
            is_bool(flat_var_error_rate),
            is_scalar_double(epsilon),
            epsilon > 0)

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

  if (! flat_var_error_rate) {
    var_error_rate <- rep(var_error_rate, n_var)
  }
  var_error_rate <- rep(0.05, n_var)
  read_error_rate <- rep(0.05, n_read)
  # note that this is the most time consuming step
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

    read_chim_lh <- matrix(chim_lh$lh, nrow = n_chim, ncol = n_read)
    for (i in seq_len(n_var)) {
      for (j in seq_len(n_read)) {
        read_chim_lh[which(read_chim_match[[i]][, j]), j] %<>% { . + inv_p_err[i, j] }
        read_chim_lh[which(!read_chim_match[[i]][, j]), j] %<>% { . + p_err[i, j] }
      }
    }

    read_chim_posterior <-
      read_chim_lh %>%
      exp() %>%
      apply(2, function(x) x / sum(x)) %>%
      matrix(ncol = ncol(read_chim_lh))

    chim_posterior <- rowSums(read_chim_posterior) / n_read

    LH_last <- LH

    LH <-
      read_chim_lh %>%
      exp() %>%
      colSums() %>%
      log() %>% sum()

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
        }) %>% mean()
      })

    read_err_posterior <- exp(log(re_mat) - p_err)

    read_error_rate <-
      future_map_dbl(seq_len(n_read), function(j) {
        map_dbl(seq_len(n_var), function(i) {
          if_else(read_chim_match[[i]][, j], 0, read_err_posterior[i, j]) %>%
            weighted.mean(w = read_chim_posterior[, j])
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
                 !!str_c('rel_lh_', n) := colSums(read_chim_posterior[which(p == 1), , drop = FALSE]))
        })
    ) %>%
    pivot_longer(-read_id,
                 names_to = c('.value', "haplotype"),
                 names_pattern = '(.+)_([^_]+)')

  read_hap_lh_smry <-
    read_hap_lh %>%
    arrange(read_id, desc(rel_lh)) %>%
    group_by(read_id) %>%
    slice(1) %>%
    ungroup()

  return(list(search = search,
              read_hap_lh = read_hap_lh,
              read_hap_lh_smry = read_hap_lh_smry))
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
