
pam_hap_k <- function(allele_matrix,
                      max_clust = 4,
                      min_ploidy = 1,
                      max_k = max_clust + 2,
                      noise_rel_asw = 0.5,
                      max_noise_prop = 1/3,
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
    mutate(cluster_asw = mean(sil_width),
           cluster_prop = n() / ncol(allele_matrix)) %>%
    group_by(k) %>%
    mutate(k_asw = mean(sil_width)) %>%
    ungroup() %>%
    group_by(k, cluster) %>%
    mutate(is_noise = (cluster_asw < noise_rel_asw * abs(k_asw)) & (cluster_prop < max_noise_prop)) %>%
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
    group_by(k) %>%
    mutate(nc = n()) %>%
    ungroup() %>%
    filter(nc <= max_clust) %>%
    select(k, cluster, nc,  medoid, medoid_state)

  # test if any ratios are plausible
  pam_ratio <-
    pam_medoid %>%
    select(k, cluster, nc) %>%
    left_join(select(pam_sil, k, cluster, prop), by = c('k', 'cluster')) %>%
    select(k, nc, prop, cluster) %>%
    chop(c(prop, cluster)) %>%
    mutate(ratio = map2(nc, prop, function(nc, prop) {
      enum_hap_ratios(nc, min_ploidy = max(nc, min_ploidy), max_ploidy = max_clust) %>%
        rowwise() %>%
        mutate(sim = cosine(ratio, prop)) %>%
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
  # collect all unique hap/ratio combinations

  k_medoids <-
    pam_medoid %>%
    inner_join(select(pam_ratio, k, cluster, ratio),
               by = c('k', 'cluster')) %>%
    # use hashes to collapse identical states
    mutate(medoid_hash = map2_chr(medoid_state, ratio, function(x, y) {
      digest::digest(x*y)
    })) %>%
    group_by(k) %>%
    mutate(k_hash = sort(medoid_hash) %>% digest::digest()) %>%
    group_by(k_hash) %>%
    filter(k == min(k)) %>%
    ungroup() %>%
    select(-ends_with('hash'))


  return(list(k_medoids = k_medoids, sil_data = pam_sil))

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

