
pam_hap_k <- function(allele_matrix,
                      max_k = 4,
                      noise_rel_asw = 0.5,
                      max_rel_diff = 0.25) {

  dissim <-
    t(unname(allele_matrix)) %>%
    as_tibble() %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')

  pam_1_medoid <- cluster::pam(dissim, k = 1)$medoid

  # summarise pam
  pam_k <-
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
    mutate(prop = if_else(is_noise, 0, lengths(rid) / sum(lengths(rid[!is_noise])))) %>%
    select(k, k_asw, cluster, medoid, cluster_asw, is_noise, prop, rid, sil_width) %>%
    ungroup() %>%
    (function(x) {
      pam_1 <- cluster::pam(dissim, k = 1)
      bind_rows(tibble(k = 1, cluster = '1', medoid = pam_1$medoids,
                       is_noise = FALSE, prop = 1,
                       rid = list(seq_len(ncol(allele_matrix))),
                       sil_width = list(rep(NA_real_, ncol(allele_matrix)))),
                x)
    }) %>%
    select(k, k_asw, cluster, medoid, cluster_asw, is_noise, prop, rid, sil_width)

  # pam_k_medoid_imp <-
  pam_k %>%
    filter(!is_noise) %>%
    select(k, cluster, medoid)

  # test if any ratios are plausible
  pam_ratios <-
    pam_k %>%
    filter(!is_noise) %>%
    select(k, prop, cluster) %>%
    group_by(k) %>%
    mutate(ek = n()) %>%
    ungroup() %>%
    chop(c(prop, cluster)) %>%
    mutate(ratio = map2(ek, prop, function(ek, prop) {
      enum_hap_ratios(ek) %>%
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
    mutate(lt_rel_max = all(rel_diff < max_rel_diff))

  # highest asw w/o nosie ???
  pam_best <-
    pam_k %>%
    group_by(k) %>%
    filter(!any(is_noise)) %>%
    ungroup() %>%
    filter(k == k[which.max(k_asw)]) %>%
    left_join(filter(pam_ratios, lt_rel_max) %>%
                select(k, cluster, ratio),
              by = c('k', 'cluster')) %>%
    select(k, cluster, medoid, prop, ratio, rid)

  return(list(pam_best = pam_best,
              pam_k = pam_k,
              pam_ratios = pam_ratios,
              pam_1_medoid = pam_1_medoid))

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
