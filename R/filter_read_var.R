

filter_read_var <- function(read_allele_mat,
                            max_na_read = 0.25,
                            max_na_var = 0.25,
                            max_maj_af = 0.95,
                            min_p_neighbour = 0.005,
                            max_neighbour_dist = 0.25) {

  v_set <- seq_len(nrow(read_allele_mat))
  r_set <- seq_len(ncol(read_allele_mat))

  while(TRUE) {
    v_mask_1 <-
      filt_vars(read_allele_mat = read_allele_mat[v_set, r_set],
                max_na_var = max_na_var,
                max_maj_af = max_maj_af) %>%
      filter(fail) %>%
      pull(vid) %>%
      { setdiff(seq_along(v_set), .) }

    r_mask_1 <-
      filt_reads(read_allele_mat = read_allele_mat[v_set, r_set],
                 max_na_read = max_na_read) %>%
      filter(fail) %>%
      pull(rid) %>%
      { setdiff(seq_along(r_set), .) }

    v_mask_2 <-
      filt_vars(read_allele_mat = read_allele_mat[v_set, r_set[r_mask_1]],
                max_na_var = max_na_var,
                max_maj_af = max_maj_af) %>%
      filter(fail) %>%
      pull(vid) %>%
      { setdiff(seq_along(v_set), .) }

    r_mask_2 <-
      filt_reads(read_allele_mat = read_allele_mat[v_set[v_mask_1], r_set],
                 max_na_read = max_na_read) %>%
      filter(fail) %>%
      pull(rid) %>%
      { setdiff(seq_along(r_set), .) }

    v_mask <- union(v_mask_1, v_mask_2) %>% sort()
    r_mask <- union(r_mask_1, r_mask_2) %>% sort()

    if (length(v_mask) == length(v_set) && length(r_mask) == length(r_set)) {
      break
    }

    r_set <- r_set[r_mask]
    v_set <- v_set[v_mask]
  }

  return(list(var_set = v_set,
              read_set = r_set))
}

filt_vars <- function(read_allele_mat,
                      max_na_var,
                      max_maj_af) {
  set_colnames(read_allele_mat, str_c('r_', seq_len(ncol(read_allele_mat)))) %>%
    as_tibble() %>%
    mutate(vid = seq_len(nrow(read_allele_mat))) %>%
    pivot_longer(-vid,
                 names_prefix = 'r_',
                 names_to = 'rid',
                 names_ptypes = list(rid = integer()),
                 values_to = 'allele') %>%
    group_by(vid) %>%
    summarise(naf = sum(is.na(allele)) / n(),
              mj_af = mjaf(allele)) %>%
    mutate(fail = naf > max_na_var | mj_af > max_maj_af)
}

mjaf <- function(x) {
  table(x) %>% { . / sum(.) } %>% max()
}

filt_reads <- function(read_allele_mat,
                       max_na_read) {
  set_colnames(read_allele_mat, str_c('r_', seq_len(ncol(read_allele_mat)))) %>%
    as_tibble() %>%
    mutate(vid = seq_len(nrow(read_allele_mat))) %>%
    pivot_longer(-vid,
                 names_prefix = 'r_',
                 names_to = 'rid',
                 names_ptypes = list(rid = integer()),
                 values_to = 'allele') %>%
    group_by(rid) %>%
    summarise(naf = sum(is.na(allele)) / n()) %>%
    mutate(fail = naf > max_na_read)
}

get_read_dist <- function(read_allele_mat) {
  n_var <- nrow(read_allele_mat)
  n_read <- ncol(read_allele_mat)
  read_dist_mat <- matrix(n_var, n_read, n_read)

  for (i in seq_len(n_var)) {
    for (a in unique(read_allele_mat[i, ])) {
      if (! is.na(a)) {
        at <- which(read_allele_mat[i, ] == a)
      } else {
        at <- which(is.na(read_allele_mat[i, ]))
      }
      read_dist_mat[at, at] %<>% {. - 1L}
    }
  }

  return(read_dist_mat)
}
