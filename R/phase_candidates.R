
get_phase_candidates <- function(allele_matrix,
                                 max_phases = 4,
                                 min_ploidy = 1,
                                 min_sim = 0.97,
                                 bimera_max_ratio = 0.5) {

  allele_matrix <- unname(allele_matrix)
  n_read <- ncol(allele_matrix)

  derep_1 <-
    tibble(read_id = seq_len(n_read),
           hash = map_chr(read_id, ~ digest(allele_matrix[, .]))) %>%
    group_by(hash) %>%
    summarise(rep_id = min(read_id),
              n = n(),
              any_na = any(is.na(allele_matrix[, rep_id])),
              read_id = list(read_id)) %>%
    arrange(desc(n)) %>%
    select(-hash)

  # imputate counts at NAs
  rep_id_no_na <- filter(derep_1, !any_na) %>% pull(rep_id)
  n_derep <- length(rep_id_no_na)
  derep_2 <-
    left_join(
      filter(derep_1, !any_na),
      filter(derep_1, any_na) %>%
        pmap_df(function(rep_id, n, ...) {
          na_at <- which(is.na(allele_matrix[, rep_id]))
          hash <- digest(allele_matrix[-na_at, rep_id])
          hashes <- map_chr(rep_id_no_na, ~ digest(allele_matrix[-na_at, .]))
          ident_at <- which(hash == hashes)
          tibble(rep_id = rep_id_no_na[ident_at],
                 n_imp = derep_1$n[match(rep_id, derep_1$rep_id)] %>%  { n * . / sum(.) }) }) %>%
        na.omit() %>%
        group_by(rep_id) %>%
        summarise(n_imp = sum(n_imp)),
      by = 'rep_id') %>%
    select(-any_na) %>%
    mutate(n = n + replace_na(n_imp, 0),
           p = n / n_read) %>%
    arrange(desc(p)) %>%
    mutate(state = map(rep_id, ~ allele_matrix[, .]))

  # choose set of top states to be candidate phase set
  phase_set <-
    derep_2 %>%
    filter(p > (max(p) / (max_phases * 1.5))) %>%
    (function(x) mutate(x, is_bimera = is_bimera(x, bimera_max_ratio = bimera_max_ratio))) %>%
    filter(!is_bimera) %>%
    pull(rep_id)

  # enumerate and filter combinations of phases base on cosine similarity to possible ratios
  # must include most abundant phases
  phases <-
    map_df(seq_len(min(max_phases, length(phase_set))), function(nc) {
      tibble(nc = nc, rep_id = phase_set[seq_len(nc)])
    }) %>%
    left_join(select(derep_2, rep_id, p), 'rep_id') %>%
    group_by(nc) %>%
    mutate(p_tot = sum(p)) %>%
    group_by(nc) %>%
    filter(if_else(nc == 1, p_tot == max(p_tot), p_tot > max(p_tot) / 2)) %>%
    select(-p_tot) %>%
    ungroup() %>%
    chop(c(rep_id, p)) %>%
    mutate(ratio_sim = map2(nc, p, function(nc, p) {
      enum_hap_ratios(nc, min_ploidy = max(nc, min_ploidy), max_ploidy = max_phases) %>%
        rowwise() %>%
        mutate(sim = cosine(ratio, p)) %>%
        slice(which.max(sim))
    })) %>%
    unnest(ratio_sim) %>%
    filter(sim >= min_sim) %>%
    select(nc, rep_id, ratio) %>%
    unnest(c(rep_id, ratio)) %>%
    left_join(select(derep_2, rep_id, state), 'rep_id')

  return(list(phases = phases, derep_state = derep_2))

}

enum_hap_ratios <- function(n_haps, min_ploidy = n_haps, max_ploidy = n_haps) {

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
      mutate(ratio = `if`(n() == 1, 1L, as.integer(hap_count / reduce(hap_count, pracma::gcd)))) %>%
      ungroup() %>%
      select(id, hap, ratio) %>%
      pivot_wider(names_from = hap, values_from = ratio) %>%
      select(-id) %>%
      distinct() }) %>%
    distinct() %>%
    mutate(id = seq_len(n())) %>%
    pivot_longer(-id, names_to = 'hap', names_ptypes = list(hap = integer()), values_to = 'ratio') %>%
    arrange(id, hap) %>%
    chop(c(hap, ratio)) %>%
    transmute(ratio = map2(ratio, hap, ~.x[.y]))
}

# return lgl_vector indicating if state is a potential bimera
is_bimera <- function(state_p, bimera_max_ratio = 0.5) {

  stopifnot(is.data.frame(state_p),
            all(c('p', 'state') %in% colnames(state_p)))

  ord <- order(state_p$p, decreasing = T)
  inv_ord <- inv_order(ord)
  state_p <- state_p[ord, ]

  state_mat <- do.call('cbind', state_p$state)
  state_hash <- map_chr(state_p$state, digest)
  n_state <- nrow(state_p)
  n_var <- nrow(state_mat)

  if (n_state <= 2) {
    return(logical(n_state))
  }

  bimera_at <-
    map_df(seq_len(n_state), function(i) {
      map_df(i + seq_len(n_state - i), function(j) {

        diff_at <- which(state_mat[, i] != state_mat[, j])
        if (length(diff_at) <= 1) {
          return(tibble(i = i, j = j))
        }

        bimera_hash <-
          map(seq_len(length(diff_at) - 1), function(at) {
            c(state_mat[seq_len(diff_at[at]), i], state_mat[diff_at[at] + seq_len(n_var - diff_at[at]), j])
          }) %>%
          map_chr(digest)

        map_df(j + seq_len(n_state - j), function(k) {
          tibble(i = i, j = j, k = k, is_bimera = state_hash[k] %in% bimera_hash)
        })
      })
    }) %>%
    filter(is_bimera) %>%
    rowwise() %>%
    mutate(mean_p_parent = mean(state_p$p[c(i, j)]),
           p_bimera = state_p$p[k]) %>%
    ungroup() %>%
    filter(p_bimera < (1 - bimera_max_ratio) * mean_p_parent) %>%
    pull(k) %>%
    sort() %>%
    unique()

  logical(n_state) %>%
    replace(bimera_at, TRUE) %>%
    { .[inv_ord] }
}
