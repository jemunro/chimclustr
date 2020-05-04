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
    map_df(function(ne) {
      list(seq_len(n_var-1)) %>%
        rep(ne) %>%
        setNames(., seq_along(.)) %>%
        do.call('expand_grid', .) %>%
        mutate(id = seq_len(n())) %>%
        gather(-id, key = 'index', value = 'vars') %>%
        arrange(id, index) %>%
        select(-index) %>%
        chop(vars) %>%
        select(-id) %>%
        filter(!map_lgl(vars, is.unsorted))
    }) %>%
    mutate(tsid = seq_len(n()) + 1,
           ts_count = map(vars, ~as.integer(table(.))),
           vars = map(vars, unique),
           ns = lengths(vars)) %>%
    mutate(vd = map2(vars, ts_count, function(v, tsc) {
      tibble(start = v, n = tsc) %>%
        { bind_rows(.,
                    tibble(start = c(1L, .$start + 1), n = 0L))
        } %>%
        arrange(start, n) %>%
        mutate(width = { lead(start) - start } %>% replace(n(), n_var - start[n()]),
               state = as.integer(cumsum(n > 0) + 1)) %>%
        filter(width > 0)
    }))

  ts_space_short <- bind_rows(
    tibble(tsid = 1L, ns = 0),
    select(ts_space, tsid, ns))


  ts_space_long <-
    select(ts_space, tsid, vd) %>%
    unnest(vd) %>%
    mutate(end = start + width,
           size = var_dist[inv_arr_ind(start, end, n_var)])

  ts_state <-
    ts_space_long %>%
    group_by(tsid) %>%
    summarise(state = list(c(1L, rep(state, width)))) %>%
    with(do.call(rbind, state)) %>%
    { rbind(1L, .) }

  ts_space_long <-
    bind_rows(tibble(tsid = 1L, n = 0, size = var_dist[1, n_var]),
              select(ts_space_long, tsid, n, size))

  chime_space <-
    full_join(hap_space, ts_space_short, 'ns') %>%
    mutate(hap_hash = map2_chr(haps, tsid, function(h, i) digest(h[ts_state[i,]]))) %>%
    group_by(hap_hash) %>%
    (function(x) {
      full_join(
        select(x, -haps) %>% nest(ids = c(hsid, tsid)) %>% ungroup(),
        slice(x, 1) %>% mutate(hap_state = map2(haps, tsid, function(h, i) h[ts_state[i,]])) %>%
          ungroup() %>% select(hap_hash, hap_state),
        by = 'hap_hash')
    }) %>%
    mutate(chid = seq_len(n())) %>%
    mutate(alleles = map(hap_state, function(hs) {
      haps[inv_arr_ind(seq_along(hs), hs, length(hs))]
    }))

  for (h in seq_len(n_hap)) {
    v <- str_c('p', h)
    chime_space <- mutate(chime_space, !!v := map_int(hap_state, ~ sum(. == h)) / n_var)
  }

  chime_space <- select(chime_space, chid, ids, ns, alleles, starts_with('p'))

  return(list(hap_space = hap_space,
              ts_space = ts_space_long,
              chime_space = chime_space))
}

chimera_lh <- function(chimeras, hap_prop, ts_rate) {

  hap_lh <-
    chimeras$hap_space %>%
    mutate(hslh = map_dbl(haps, ~ sum(log(hap_prop[.])))) %>%
    select(hsid, hslh)

  ts_lh <-
    chimeras$ts_space %>%
    mutate(tslh = dbinom(n, size, ts_rate, log = TRUE)) %>%
    group_by(tsid) %>%
    summarise(tslh = sum(tslh))

  chime_lh <-
    chimeras$chime_space %>%
    select(chid, ids) %>%
    unnest(ids) %>%
    left_join(hap_lh, 'hsid') %>%
    left_join(ts_lh, 'tsid') %>%
    group_by(chid) %>%
    summarise(lh = log_sum(hslh + tslh))

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
hap_mix_em <- function(read_allele_mat,
                       var_pos,
                       haps,
                       hap_prop,
                       ts_rate,
                       error_rate,
                       max_iter = 50,
                       fixed_hap_prop = FALSE,
                       flat_error_rate = FALSE,
                       epsilon = 1e-3) {

  stopifnot(is.matrix(read_allele_mat),
            is_integer(read_allele_mat),
            is_integer(var_pos),
            length(var_pos) == nrow(read_allele_mat),
            is.matrix(haps) && ncol(haps) > 1L,
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

  n_var <- length(var_pos)
  n_read <- ncol(read_allele_mat)
  var_width <- last(var_pos) - first(var_pos)
  if (! flat_error_rate) {
    error_rate <- rep(error_rate, n_var)
  }
  # note that this is the most time consuming step
  chims <- enum_chimeras(haps, var_pos = var_pos, ts_max = 3)
  n_chim <- nrow(chims$chime_space)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)

  read_chim_dist <- read_chimera_dist(read_allele_mat,
                                      chim_alleles)
  read_chim_ident <- n_var - read_chim_dist
  read_chim_match <- read_chimera_match(read_allele_mat, chim_alleles)

  n_iter <- 0L
  LH <- -Inf
  res <- tibble(iter = integer(),
                LH = double(),
                ts_rate = double(),
                error_rate = list(),
                hap_prop = list())

  while(n_iter < max_iter) {
    n_iter <- n_iter + 1L

    ## Expecation
    chim_lh <- chimera_lh(chims, hap_prop, ts_rate)

    # log lh
    if (flat_error_rate || all(error_rate[1] == error_rate)) {
      read_chim_lh <-
        (read_chim_dist * log(error_rate[1])) +
        (read_chim_ident * log(1-error_rate[1])) +
        chim_lh$lh
    } else {
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
      apply(2, function(x) x / sum(x))

    # lh
    read_chim_lh_norm_rs <- rowSums(read_chim_lh_norm)

    LH_last <- LH

    LH <-
      read_chim_lh %>%
      exp() %>%
      colSums() %>%
      log() %>% sum()

    # record results
    res <- add_row(res,
                   iter = n_iter,
                   LH = LH,
                   ts_rate = ts_rate,
                   error_rate = list(error_rate),
                   hap_prop = list(hap_prop))

    if (epsilon > LH - LH_last) {
      break
    }

    ## Maximisation
    if (!fixed_hap_prop) {
      hap_prop <-
        select(chims$chime_space, starts_with('p')) %>%
        map_dbl(~ weighted.mean(., read_chim_lh_norm_rs))
    }

    ts_rate <- weighted.mean(chims$chime_space$ns, read_chim_lh_norm_rs) / var_width

    if (flat_error_rate) {
       error_rate <- weighted.mean(read_chim_dist, read_chim_lh_norm) / n_var
    } else {
      error_rate <- map_dbl(seq_len(n_var), function(i) {
        read_chimera_match_sub(match_list = read_chim_match[[i]],
                               n_read = n_read,
                               n_chim = n_chim,
                               x0 = 1,
                               x1 = 0) %>%
          weighted.mean(read_chim_lh_norm)
      })
    }
  }

  return(res)
}

#' @importFrom rlang is_integerish
#' @importFrom magrittr '%<>%'
a_mult_d3 <- function(a, x) {
  stopifnot(is.array(a),
            is.vector(x),
            length(dim(a)) >= 3,
            dim(a)[3] == length(x))
  for (i in seq_along(x)) {
    a[ , , i] %<>% { . * x[i] }
  }
  return(a)
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

  # em_res_1 <- hap_mix_em(read_allele_mat = read_allele_mat,
  #                        var_pos = var_pos,
  #                        haps = haps,
  #                        hap_prop = c(0.5, 0.5),
  #                        ts_rate = 1e-4,
  #                        error_rate = 1e-1,
  #                        max_iter = 50,
  #                        fixed_hap_prop = FALSE,
  #                        flat_error_rate = FALSE)
  #
  # em_res_1f <- hap_mix_em(read_allele_mat = read_allele_mat,
  #                         var_pos = var_pos,
  #                         haps = haps,
  #                         hap_prop = c(0.5, 0.5),
  #                         ts_rate = 1e-4,
  #                         error_rate = 1e-1,
  #                         max_iter = 50,
  #                         fixed_hap_prop = TRUE,
  #                         flat_error_rate = TRUE)

  em_res_1v <- hap_mix_em(read_allele_mat = read_allele_mat,
                         var_pos = var_pos,
                         haps = haps,
                         hap_prop = c(0.5, 0.5),
                         ts_rate = 1e-4,
                         error_rate = 1e-1,
                         max_iter = 50,
                         fixed_hap_prop = FALSE,
                         flat_error_rate = FALSE)

  em_res_1fv <- hap_mix_em(read_allele_mat = read_allele_mat,
                           var_pos = var_pos,
                           haps = haps,
                           hap_prop = c(0.5, 0.5),
                           ts_rate = 1e-4,
                           error_rate = 1e-1,
                           max_iter = 50,
                           fixed_hap_prop = TRUE,
                           flat_error_rate = FALSE)

  em_res_1fv_21 <- hap_mix_em(read_allele_mat = read_allele_mat,
                              var_pos = var_pos,
                              haps = haps,
                              hap_prop = c(0.67, 0.33),
                              ts_rate = 1e-4,
                              error_rate = 1e-1,
                              max_iter = 50,
                              fixed_hap_prop = TRUE,
                              flat_error_rate = FALSE)

  em_res_1fv_12 <- hap_mix_em(read_allele_mat = read_allele_mat,
                              var_pos = var_pos,
                              haps = haps,
                              hap_prop = c(0.33, 0.67),
                              ts_rate = 1e-4,
                              error_rate = 1e-1,
                              max_iter = 50,
                              fixed_hap_prop = TRUE,
                              flat_error_rate = FALSE)

}
