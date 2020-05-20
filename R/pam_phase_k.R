
pam_hap_k <- function(allele_matrix,
                      max_clust = 4,
                      max_k = max_clust + 2,
                      noise_rel_asw = 0.5,
                      max_rel_diff = 0.25) {

  dissim <-
    t(unname(allele_matrix)) %>%
    as_tibble() %>%
    mutate_all(~ replace_na(., -1L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')

  pam_sil <-
    seq.int(2, max_k) %>%
    map_df(function(k) tibble(k = k, pamobj = list(cluster::pam(dissim, k = k)))) %>%
    mutate(sildata = map(pamobj, function(po) {
      as.data.frame(po$silinfo$widths) %>%
        rownames_to_column('rid') %>%
        mutate(rid = as.integer(rid)) %>%
        as_tibble() %>%
        mutate(cluster = as.character(as.integer(cluster))) %>%
        select(-neighbor) %>%
        mutate(medoid = po$medoids[as.integer(cluster)])
    })) %>%
    select(-pamobj) %>%
    unnest(sildata) %>%
    group_by(k, cluster) %>%
    mutate(cluster_asw = mean(sil_width)) %>%
    group_by(k) %>%
    mutate(k_asw = mean(sil_width)) %>%
    ungroup() %>%
    group_by(k, cluster) %>%
    mutate(is_noise = cluster_asw < noise_rel_asw * abs(k_asw)) %>%
    group_by(k) %>%
    mutate(k_asw_filt = mean(sil_width[!is_noise])) %>%
    chop(c(rid, sil_width)) %>%
    arrange(desc(cluster_asw)) %>%
    mutate(cluster = as.character(seq_along(cluster))) %>%
    mutate(prop = if_else(is_noise, 0, lengths(rid) / sum(lengths(rid[!is_noise])))) %>%
    ungroup() %>%
    (function(x) {
      pam_1 <- cluster::pam(dissim, k = 1)
      bind_rows(tibble(k = 1, cluster = '1', medoid = pam_1$medoids,
                       is_noise = FALSE, prop = 1,
                       rid = list(seq_len(ncol(allele_matrix))),
                       sil_width = list(rep(NA_real_, ncol(allele_matrix)))),
                x)
    }) %>%
    select(k, k_asw, k_asw_filt, cluster, medoid, cluster_asw, is_noise, prop, rid, sil_width) %>%
    arrange(k, as.integer(cluster))

  pam_medoid <-
    pam_sil %>%
    filter(k <= k[which.max(k_asw_filt)],
           !is_noise) %>%
    select(k, cluster, medoid, rid) %>%
    mutate(medoid_state = map2(medoid, rid, function(medoid, rid) {
      medoid_impute_missing(allele_matrix[, medoid], allele_matrix[, rid, drop = FALSE])
    })) %>%
    select(-rid) %>%
    mutate(medoid_hash = map_chr(medoid, digest::digest)) %>%
    group_by(k) %>%
    mutate(nc = n(),
           k_hash = sort(medoid_hash) %>% digest::digest()) %>%
    group_by(k_hash) %>%
    filter(k == min(k),
           nc <= max_clust) %>%
    ungroup() %>%
    select(k, cluster, nc,  medoid, medoid_state)

  # test if any ratios are plausible
  pam_ratio <-
    pam_medoid %>%
    select(k, cluster, nc) %>%
    left_join(select(pam_sil, k, cluster, prop), by = c('k', 'cluster')) %>%
    select(k, nc, prop, cluster) %>%
    chop(c(prop, cluster)) %>%
    mutate(ratio = map2(nc, prop, function(nc, prop) {
      enum_hap_ratios(nc) %>%
        rowwise() %>%
        mutate(sim = lsa::cosine(ratio, prop)) %>%
        arrange(desc(sim)) %>%
        pull(ratio) %>%
        first()
    })) %>%
    mutate(norm_ratio = map(ratio, ~ . / sum(.))) %>%
    unnest(c(prop, ratio, cluster, norm_ratio)) %>%
    mutate(rel_diff = abs(prop - norm_ratio) / norm_ratio) %>%
    group_by(k) %>%
    filter(all(rel_diff < max_rel_diff)) %>%
    ungroup()

  # highest asw w/o nosie ???
  pam_best <-
    pam_medoid %>%
    inner_join(select(pam_ratio, k, cluster, ratio),
               by = c('k', 'cluster'))

  return(list(pam_best = pam_best,
              pam_sil = pam_sil))

}

medoid_impute_missing <- function(med_state, state_matrix) {
  # impute with similarity weighted mode

  stopifnot(is.integer(med_state),
            is.matrix(state_matrix),
            is.integer(state_matrix),
            length(med_state) == nrow(state_matrix))

  if (!any(is.na(med_state))) {
    return(med_state)
  }

  sim_med <-
    t(unname(state_matrix)) %>%
    { rbind(med_state, .) } %>%
    as_tibble() %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower') %>%
    as.matrix() %>%
    { 1 - .[-1, 1] }

  repl <-
    map_int(which(is.na(med_state)), function(i) {
      at <- which(!is.na(state_matrix[i, ]))
      rowsum(sim_med[at], state_matrix[i, at]) %>%
        { as.integer(rownames(.))[which.max(.)] }
    })

  return(replace(med_state, which(is.na(med_state)), repl))
}


enum_hap_ratios <- function(n_haps, min_ploidy = n_haps, max_ploidy = 4L) {

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
