# obs. vars. (y) = read variant alleles
# lat. vars. (z) = haplotype at each variant
# parameters     = p_err, p_chim, p_hap_i

require(tidyverse)
require(magrittr)

log_add <- function (x, y) {
  .max <- pmax(x, y)
  .min <- pmin(x, y)
  .max + log1p(exp(.min - .max))
}

#' @importFrom purrr reduce
log_sum <- function(x) {
  reduce(sort(x, TRUE), function (a, b) a + log1p(exp(b-a)))
}


inv_arr_ind <- function(row_idx, col_idx, nrow) {

  stopifnot(length(row_idx) == length(col_idx),
            max(row_idx) <= nrow,
            rlang::is_integerish(row_idx),
            rlang::is_integerish(col_idx),
            rlang::is_scalar_integerish(nrow))

  (col_idx - 1L) * nrow + row_idx
}


test_data <- function() {
  require(SeqVarTools)
  require(SeqVarTools)
  require(tidyverse)
  gds <- seqOpen('test_data/SM-NA10005.AM-CYP2D6_6k_v2.htc_1.gds', allow.duplicate = T)

  read_tsv('test_data/SM-NA10005.AM-CYP2D6_6k_v2.haplotag.tsv.gz',
           col_types = cols(PS = col_integer(), HP = col_integer())) %>%
    mutate(vid = as.integer(PS + 1)) %>%
    select(rid = name, vid, allele = HP) %>%
    nest(read_data = c(rid, allele)) %>%
    left_join(variantInfo(gds) %>%
                as_tibble() %>%
                select(vid = variant.id, pos),
              by='vid') %>%
    select(vid, pos, read_data) %>%
    (function(x) {
      test_matrix(x) %>%
        as_tibble() %>%
        mutate(vid = seq_len(n())) %>%
        gather(-vid, key = 'rid', value = 'allele') %>%
        group_by(vid) %>%
        summarise(n_allele = n_distinct(allele)) %>%
        { left_join(x, ., 'vid') }
    }) %>%
    filter(n_allele > 1) %>%
    mutate(vid = as.integer(as.factor(vid))) %>%
    arrange(vid, pos)
}

test_matrix <- function(data = test_data()) {
  data %>%
    select(vid, read_data) %>%
    unnest(read_data) %>%
    spread(rid, allele) %>%
    arrange(vid) %>%
    { as.matrix(select(., -vid)) }
}

## also consider genotyping with gatk ploidy set to max CN,
# & only keep hets

##
# given a set of haplotypes, we only care about variants that differ between these haplotypes
# i.e. no need to enumerate over others

##
# another simplification
# when haplotype remains the same, the position of the switch is not important
# i.e. can be take over whole range
# should speed up calc

#' @importFrom digest digest
enum_chimeras <- function(haps, var_pos, ts_max) {

  n_hap <- ncol(haps)
  n_var <- length(var_pos)
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

  return(list(hap_space = hap_space,
              ts_space = ts_space,
              chime_space = chime_space))
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
    summarise(tslh = sum(lh), en = sum(en))


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

derep_allele_set <- function(allele_mat, id = 'id') {
  tibble(index = seq_len(ncol(allele_mat))) %>%
    mutate(hash = map_chr(index, ~ digest(allele_mat[, .]))) %>%
    chop(index) %>%
    mutate(rep = map_int(index, first),
           !! id := seq_len(n()),
           n = lengths(index)) %>%
    select(!! id, rep, n, indices = index)
}


#' @importFrom rlang is_bool is_integer is_scalar_double
#' @importFrom dplyr first last mutate
hap_mix_em <- function(read_allele_mat,
                       var_pos,
                       haps,
                       hap_prop,
                       ts_rate,
                       ts_max = 3,
                       error_rate,
                       flat_error_rate = FALSE,
                       max_error_rate = 0.5,
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
            is_scalar_double(error_rate),
            error_rate > 0 && error_rate < 1,
            is_bool(fixed_hap_prop),
            is_bool(flat_error_rate),
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

  if (! flat_error_rate) {
    error_rate <- rep(error_rate, n_var)
  }
  # note that this is the most time consuming step
  chims <- enum_chimeras(haps, var_pos = var_pos, ts_max = ts_max)
  n_chim <- nrow(chims$chime_space)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)

  read_chim_dist <- read_chimera_dist(read_allele_mat, chim_alleles)
  read_chim_ident <- n_var - read_chim_dist
  read_chim_match <- read_chimera_match(read_allele_mat, chim_alleles)

  n_iter <- 0L
  LH <- -Inf
  search <- tibble(iter = integer(),
                LH = double(),
                ts_rate = double(),
                error_rate = list(),
                hap_prop = list())

  while(n_iter < max_iter) {
    n_iter <- n_iter + 1L

    ## Expecation
    chim_lh <- chimera_lh(chims, hap_prop, ts_rate, var_width = var_width)

    # log lh
    if (flat_error_rate || all(error_rate[1] == error_rate)) {
      read_chim_lh <-
        (read_chim_dist * log(error_rate[1])) +
        (read_chim_ident * log(1-error_rate[1])) +
        chim_lh$lh
    } else {
      # slow - good candidate for optimisation
      read_chim_lh <-
        map(seq_len(n_var), function(i) {
          read_chimera_match_sub(match_list = read_chim_match[[i]],
                                 n_read = n_read,
                                 n_chim = n_chim,
                                 x0 = log(error_rate[i]),
                                 x1 = log(1 - error_rate[i]))
        }) %>%
        reduce(`+`) +
        chim_lh$lh
    }

    read_chim_lh_norm <-
      read_chim_lh %>%
      exp() %>%
      apply(2, function(x) x / sum(x)) %>%
      matrix(ncol = ncol(read_chim_lh))

    # lh
    read_chim_lh_norm_rs <- rowSums(read_chim_lh_norm)

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
                      error_rate = list(error_rate),
                      hap_prop = list(hap_prop))

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
        map_dbl(~ weighted.mean(., read_chim_lh_norm_rs))
    }

    if (n_hap > 1) {
      ts_rate <- weighted.mean(chim_lh$en, read_chim_lh_norm_rs) / var_width
    }

    if (flat_error_rate) {
      error_rate <- weighted.mean(read_chim_dist, read_chim_lh_norm) / n_var
    } else {
      # slow - good candidate for optimisation
      error_rate <- map_dbl(seq_len(n_var), function(i) {
        read_chimera_match_sub(match_list = read_chim_match[[i]],
                               n_read = n_read,
                               n_chim = n_chim,
                               x0 = 1,
                               x1 = 0) %>%
          weighted.mean(read_chim_lh_norm)
      })
    }
    # error_rate > 0.5 is meaningless
    error_rate <- pmin(max_error_rate, error_rate)
  }

  read_hap_lh <-
    tibble(read_id = seq_len(n_read)) %>%
    bind_cols(
      select(chims$chime_space, starts_with('p')) %>%
        mutate(pchimera = if_else(rowSums(. > 0 & . < 1) > 0, 1, 0)) %>%
        map2_dfc(str_remove(names(.), 'p'), ., function(n, p) {
          tibble(!!str_c('raw_llh_', n) := log(colSums(exp(read_chim_lh[which(p == 1), , drop = FALSE]))),
                 !!str_c('rel_lh_', n) := colSums(read_chim_lh_norm[which(p == 1), , drop = FALSE]))
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


enum_hap_ratios <- function(n_haps,
                            max_ploidy,
                            min_ploidy = n_haps) {


  if (n_haps == 1) {
    tibble(ratio = list(1L))
  } else {
    min_ploidy <- max(min_ploidy, n_haps)
    map_df(seq.int(min_ploidy, max_ploidy), function(ploidy) {
      rep(list(seq_len(n_haps)), ploidy) %>%
        setNames(seq_len(ploidy)) %>%
        do.call(expand_grid, .) %>%
        mutate(id = seq_len(n())) %>%
        pivot_longer(-id, names_to = 'copy', values_to = 'hap') %>%
        group_by(id, hap) %>%
        summarise(hap_count = n()) %>%
        group_by(id) %>%
        filter(n_distinct(hap) == n_haps,
               sum(hap_count) == ploidy) %>%
        mutate(ratio = `if`(n() == 1, 1L, hap_count / reduce(hap_count, pracma::gcd))) %>%
        ungroup() %>%
        select(id,hap, ratio) %>%
        pivot_wider(names_from = hap, values_from = ratio) %>%
        select(-id) %>%
        distinct() }) %>%
      distinct() %>%
      mutate(id = seq_len(n())) %>%
      pivot_longer(-id, names_to = 'hap', values_to = 'ratio') %>%
      mutate(hap = as.integer(hap), ratio = as.integer(ratio)) %>%
      chop(c(hap, ratio)) %>%
      transmute(ratio = map2(ratio, hap, ~.x[.y]))
  }
}


chimclustr <- function(read_allele_mat,
                       var_pos,
                       max_ploidy = 4L) {
  #
  read_allele_df <-
    t(unname(read_allele_mat)) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_all(as.factor) %>%
    as.data.frame() %>%
    column_to_rownames()

  read_dissim <- cluster::daisy(read_allele_df, metric = 'gower')
  read_dissim_mat <- as.matrix(read_dissim)
  ro <- hclust(read_dissim)$order

  t_read_allele_df <-
    unname(read_allele_mat) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_all(as.factor) %>%
    as.data.frame() %>%
    column_to_rownames()

  var_dissim <- cluster::daisy(t_read_allele_df, metric = 'gower')
  vo <- hclust(var_dissim)$order

  pam_meds <-
    seq_len(max_ploidy) %>%
    map_df(function(k) {
      tibble(k = k,
             medoid = as.integer(cluster::pam(read_dissim, k = k)$medoids))
    }) %>%
    mutate(medoid = map_int(medoid, ~ min(which(read_dissim_mat[., ] == 0))))

  var_set <-
    read_allele_mat[, unique(pam_meds$medoid)] %>%
    apply(1, function(x) n_distinct(x) > 1) %>%
    which()

  model_search <-
    pam_meds %>%
    chop(medoid) %>%
    mutate(ratios = map(k, ~ enum_hap_ratios(n_haps = .,
                                             max_ploidy = max_ploidy))) %>%
    unnest(ratios) %>%
    mutate(em_res = pmap(., function(medoid, ratio, ...) {
      message('ratio: ', str_c(ratio, collapse = ':'))
      chimclustr:::hap_mix_em(read_allele_mat = read_allele_mat[var_set, , drop = FALSE],
                             var_pos = var_pos[var_set],
                             haps = read_allele_mat[var_set, medoid, drop = FALSE],
                             hap_prop = ratio / sum(ratio),
                             ts_rate = 1e-4,
                             ts_max = 2,
                             error_rate = 1e-1,
                             max_iter = 20,
                             epsilon = 0.1,
                             fixed_hap_prop = TRUE,
                             flat_error_rate = FALSE)
    })) %>%
    mutate(LH = map_dbl(em_res, ~ last(.$search$LH))) %>%
    saveRDS('model_search.rds')

  model_search <-
    read_rds('model_search.rds') %>%
    mutate(multinom_lh = map2_dbl(ratio, em_res, function(ratio, em_res) {
      hap_counts <-
        em_res$read_hap_lh_smry %>%
        filter(haplotype != 'chimera') %>%
        count(haplotype) %>%
        mutate(haplotype = as.integer(haplotype)) %>%
        arrange(haplotype) %>%
        pull(n)

      chim_count <-
        em_res$read_hap_lh_smry %>%
        filter(haplotype == 'chimera') %>%
        nrow()

      size <- sum(hap_counts) + chim_count

      chim_prop <- chim_count / (chim_count + sum(hap_counts))
      hap_prop <- ratio * ((1 - chim_prop) / sum(ratio))
      exp_hap_count <- round(hap_prop * size)
      delta_exp <- round(abs(hap_counts - exp_hap_count))

      suppressWarnings(
        pmultinom::pmultinom(lower = c(exp_hap_count - delta_exp, -Inf),
                             upper = c(exp_hap_count + delta_exp - 1, Inf),
                             size = size,
                             probs = c(hap_prop, chim_prop),
                             method = 'exact') %>%
          { log (1-.) }
      )
    })) %>%
    mutate(multinom_lh_2 = map2_dbl(ratio, em_res, function(ratio, em_res) {
      hap_counts <-
        em_res$read_hap_lh %>%
        filter(haplotype != 'chimera') %>%
        group_by(read_id) %>%
        summarise(haplotype = haplotype[which.min(raw_llh)]) %>%
        count(haplotype) %>%
        mutate(haplotype = as.integer(haplotype)) %>%
        arrange(haplotype) %>%
        pull(n)

      size <- sum(hap_counts)
      hap_prop <- ratio / sum(ratio)
      exp_hap_count <- round(hap_prop * size)
      delta_exp <- round(abs(hap_counts - exp_hap_count))

      suppressWarnings(
        pmultinom::pmultinom(lower = c(exp_hap_count - delta_exp),
                             upper = c(exp_hap_count + delta_exp - 1),
                             size = size,
                             probs = c(hap_prop),
                             method = 'exact') %>%
          { log (1-.) }
      )
    }))


  read_sample <- sample(ro, 500) %>% as.character()

  read_heatmap <-
    read_dissim_mat %>%
    as.data.frame() %>%
    rownames_to_column('row_id') %>%
    as_tibble() %>%
    gather(-row_id, key = 'col_id', value = 'sim') %>%
    mutate(sim = 1 - sim) %>%
    filter(row_id %in% read_sample, col_id %in% read_sample) %>%
    mutate(col_id = factor(col_id, levels = as.character(ro)),
           row_id = factor(row_id, levels = as.character(rev(ro)))) %>%
    ggplot() +
    geom_tile(aes(row_id, col_id, fill = sim)) +
    scale_fill_viridis_c(option = 'E') +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.x = element_blank(), legend.position = 'top') +
    ylab('read') +
    xlab('read')
  # read_heatmap

  ra_plot <-
    read_allele_df %>%
    { set_colnames(., 1:ncol(.)) } %>%
    as_tibble() %>%
    mutate_all(as.integer) %>%
    mutate(read_id = as.character(seq_len(n()))) %>%
    filter(read_id %in% read_sample) %>%
    gather(-read_id, key = 'var_id', value = 'allele') %>%
    mutate(read_id = factor(read_id, levels = as.character(ro)),
           var_id = factor(var_id, levels = as.character(rev(vo)))) %>%
    mutate(allele = as.character(allele)) %>%
    ggplot() +
    geom_tile(aes(var_id, read_id, fill = allele)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
    guides(fill = F) +
    xlab('variant') +
    scale_fill_manual(values =  c(`1` = "#FFEA46FF", `2`= "#00204DFF"))

  for( i in (seq_along(model_search$em_res))) {
    hap_plot <-
      model_search$em_res[[i]]$read_hap_lh_smry %>%
      mutate(read_id = as.character(read_id)) %>%
      filter(read_id %in% read_sample) %>%
      mutate(read_id = factor(read_id, ro)) %>%
      ggplot() +
      geom_tile(aes(x=1, y=read_id, fill = haplotype)) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title.y = element_blank(), legend.position = 'top',
            axis.ticks.x = element_blank()) +
      scale_fill_brewer(palette = 'Dark2') +
      xlab('haplotype')

    hap_lh_plot <-
      model_search$em_res[[i]]$read_hap_lh %>%
      mutate(read_id = as.character(read_id)) %>%
      filter(read_id %in% read_sample) %>%
      mutate(read_id = factor(read_id, ro)) %>%
      ggplot() +
      geom_tile(aes(x=haplotype, y=read_id, fill = raw_llh)) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      scale_fill_viridis_c() +
      guides(fill = guide_colorbar(title = 'log likelihood')) +
      ylab('read') +
      xlab('haplotype')


    plot_grid(read_heatmap, ra_plot, hap_plot, hap_lh_plot,
              rel_widths = c(15, 5, 3, 6),
              nrow = 1, align = 'h', axis = 'tblr')

    ggsave(str_c(i, '.png'), width = 11, height = 8)
  }




}

run_test <- function() {

  # would be good to filter outliers here...
  # based on dist from medoids
  # but we need selected reads to be present across different clusterings for LH to be comparable (similarly all vars)
  # soo... whithin 1/0.5 * max_center_dist

  # add modes for fixed hap prop - this will be better for assigning discrete ploidies
  # output should be read_id + ML state ? (e.g. cluster_1:0.95, cluster_2:0.05, error: 0.25)
  #     or fuzzy state - e.g. p(cluster_1), p(cluster_2), p(chimera)
  #     or simply cluster_id + cluster_dist

  # might neeed to include a "noise" haplotype
  # where allele probabilities is estimated from obs alleles???

  # could use ‘poisbinom’ to calculate likelhood of x errors - useful for filtering ???
  # is in necessary to filter at all, or just assign reads?

  # can we also include a read error rate in the model??
  # this would allow some reads to be --silenced--

  data <- test_data()
  read_allele_mat <- test_matrix(data)
  read_allele_mat[is.na(read_allele_mat)] <- -1L

  read_allele_tbl <-
    t(read_allele_mat) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_all(as.factor) %>%
    as.data.frame() %>%
    column_to_rownames()

  disim <- cluster::daisy(read_allele_tbl, metric = 'gower')
  pam <- cluster::pam(disim, k = 2)

  meds <- pam$medoids

  vars_diff <-
    read_allele_mat[, meds] %>%
    apply(1, function(x) n_distinct(x) > 1) %>%
    which()

  read_allele_mat <- read_allele_mat[vars_diff, ]

  haps <-
    as.integer(read_allele_mat[, meds]) %>%
    matrix(nrow = nrow(read_allele_mat))

  var_pos <- data$pos[vars_diff]


  em_res_1 <- hap_mix_em(read_allele_mat = read_allele_mat,
                         var_pos = var_pos,
                         haps = haps,
                         hap_prop = c(0.5, 0.5),
                         ts_rate = 1e-4,
                         error_rate = 1e-1,
                         max_iter = 10,
                         epsilon = 1e-1,
                         fixed_hap_prop = TRUE,
                         flat_error_rate = FALSE)

  # disim_1 <- cluster::daisy(read_allele_tbl[em_res_1$read_hap_lh_smry$haplotype == '1', ],
                            # metric = 'gower')

  em_res_mono <- hap_mix_em(read_allele_mat = read_allele_mat,
                            var_pos = var_pos,
                            haps = haps[, 1, drop = FALSE],
                            hap_prop = 1,
                            ts_rate = 1e-4,
                            error_rate = 1e-1,
                            max_iter = 10,
                            epsilon = 1e-1,
                            fixed_hap_prop = TRUE,
                            flat_error_rate = FALSE)


}
