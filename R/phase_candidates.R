
get_phase_candidates <- function(allele_matrix,
                                 max_phases = 4,
                                 min_ploidy = 1,
                                 min_sim = 0.97) {

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
  phase_set <- derep_2 %>% filter(p > (max(p) / (max_phases * 1.5))) %>% pull(rep_id)

  # enumerate and filter combinations of phases base on cosine similarity to possible ratios
  phases <-
    `if`(length(phase_set) == 1,
         tibble(nc = 1L, id = 1L, rep_id = phase_set),
         seq.int(1, min(max_phases, length(phase_set))) %>%
           map_df(function(nc) {
             t(combn(phase_set, nc)) %>%
               set_colnames(str_c('I_', seq_len(ncol(.)))) %>%
               as_tibble() %>%
               mutate(id = seq_len(n())) %>%
               pivot_longer(-id,
                            names_prefix = 'I_',
                            names_ptypes = list(index = integer()),
                            names_to = 'index',
                            values_to = 'rep_id') %>%
               select(-index) %>%
               mutate(nc = nc) %>%
               select(nc, everything())
           })
    ) %>%
    left_join(select(derep_2, rep_id, p), 'rep_id') %>%
    group_by(id, nc) %>%
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
    select(nc, id, rep_id, ratio) %>%
    mutate(id = seq_along(id)) %>%
    unnest(c(rep_id, ratio)) %>%
    left_join(select(derep_2, rep_id, state), 'rep_id')

  return(list(phases = phases, derep_state = derep_2))

 }
