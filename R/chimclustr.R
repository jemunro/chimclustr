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
    select(-pos) %>%
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
    select(-ns) %>%
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

  chime_space <- select(chime_space, chid, ids, alleles, starts_with('p'))

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
                              chimera_alleles) {

  expand_grid(r_index = seq_len(ncol(read_alleles)),
              c_index = seq_len(ncol(chimera_alleles))) %>%
    mutate(dist = map2_int(r_index, c_index,
                           ~ sum(read_alleles[, .x] != chimera_alleles[, .y])))

}

derep_allele_set <- function(allele_mat) {
  tibble(index = seq_len(ncol(allele_mat))) %>%
    mutate(hash = map_chr(index, ~ digest(allele_mat[, .]))) %>%
    group_by(hash) %>%
    mutate(group = index[1]) %>%
    add_count() %>%
    ungroup() %>%
    select(-hash)
}


run_test <- function() {

  ## init
  data <- test_data()
  read_allele_mat <- test_matrix(data)
  read_allele_mat[is.na(read_allele_mat)] <- -1


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

  vars_diff <- which(read_allele_mat[, meds[1]] != read_allele_mat[, meds[2]])
  read_allele_mat <- read_allele_mat[vars_diff, ]
  read_allele_derep_key <- derep_allele_set(read_allele_mat)
  read_allele_derep <- read_allele_mat[, unique(read_allele_derep_key$group)]

  haps <-
    as.integer(read_allele_mat[, meds]) %>%
    matrix(nrow = nrow(read_allele_mat))

  var_pos <- data$pos[vars_diff]
  hap_prop <- c(0.49, 0.51)
  ts_rate <- 1/10000
  ts_max <- 3
  error_rate <- 1/100

  chims <- enum_chimeras(haps, var_pos = var_pos, ts_max = 2)
  chim_alleles <- chims$chime_space$alleles %>% do.call('cbind', .)
  chim_allele_derep_key <- derep_allele_set(chim_alleles)
  chim_allele_derep <- chim_alleles[, unique(chim_allele_derep_key$group)]

  read_chim_dist <- read_chimera_dist(read_allele_derep, chim_allele_derep)

  ## Expecation
  read_chim_lh <- mutate(read_chim_dist, lh = dist * log(error_rate))
  chims_lh <- chimera_lh(chims, hap_prop, ts_rate)
  # expand to matrix n_read by n_chimera of complete lh

  ## Maximisation


  ## compare M step speed with Rcgmin, optim(GC) and LBFGS
}
