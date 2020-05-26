
#' @importFrom dplyr first last mutate mutate_all select filter arrange desc
#' @importFrom magrittr '%>%' set_rownames set_colnames
#' @importFrom tibble tibble as_tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @importFrom purrr map pmap map_dbl map_chr map_int map_lgl reduce
#' @importFrom rlang is_integerish is_scalar_double

chimclustr <- function(allele_matrix,
                       var_pos,
                       max_copy_num,
                       max_copy_num,
                       var_max_miss_rate,
                       read_max_miss_rate,
                       read_min_lh,
                       read_max_error_rate) {

  stopifnot(is.matrix(allele_matrix),
            is.integer(allele_matrix),
            is_integerish(var_pos),
            length(var_pos) == nrow(allele_matrix),
            is_scalar_double(var_max_miss_rate),
            is_scalar_double(read_max_miss_rate),
            is_scalar_double(read_min_lh),
            is_scalar_double(read_max_error_rate))

  missingness <- marginal_rate_em(is.na(allele_matrix))
  read_pass <- which(missingness$col_rate < 0.05)
  var_pass <- which(missingness$row_rate < 0.25)

  plot_1_data <-
    list(mat = t(allele_matrix),
         row_annot =
           tibble(rid = colnames(allele_matrix),
                  status = if_else(seq_along(rid) %in% read_pass, 'pass', 'fail')) %>%
           as.data.frame() %>%
           column_to_rownames('rid'),
         col_annot =
           tibble(vid = rownames(allele_matrix),
                  status = if_else(seq_along(vid) %in% var_pass, 'pass', 'fail')) %>%
           as.data.frame() %>%
           column_to_rownames('vid'),
         row_hclust =
           t(allele_matrix) %>%
           as_tibble() %>%
           mutate_all(~ replace_na(., -1L)) %>%
           mutate_all(as.factor) %>%
           cluster::daisy(metric = 'gower') %>%
           hclust(method = 'average'))

  read_allele_flt <- allele_matrix[var_pass, read_pass]
  pam_hap_k <- pam_hap_k(read_allele_flt, max_clust = 4)

  plot_2_data <-
    list(mat = t(read_allele_flt) %>% set_rownames(seq_len(nrow(.))),
         row_annot =
           pam_hap_k$sil_data %>%
           filter(k > 1) %>%
           select(k, cluster, rid) %>%
           unnest(rid) %>%
           pivot_wider(names_from = k, values_from = cluster, names_prefix = 'k_') %>%
           arrange(rid) %>%
           as.data.frame() %>%
           column_to_rownames('rid'),
         row_hclust =
           t(read_allele_flt) %>%
           as_tibble() %>%
           mutate_all(~ replace_na(., -1L)) %>%
           mutate_all(as.factor) %>%
           cluster::daisy(metric = 'gower') %>%
           hclust(method = 'average'))

  plot_2_data$annot_color <-
    plot_2_data$row_annot %>%
    map(unique) %>%
    map(~ set_names(scales::hue_pal()(length(.)), .))

  ## visualise clustering silhouette
  pam_hap_k$sil_data %>%
    filter(k > 1) %>%
    select(k, cluster, rid, sil_width, is_noise, k_asw) %>%
    mutate(lab = str_c('k = ', k,', asw = ', format(k_asw, digits = 2))) %>%
    unnest(c(rid, sil_width)) %>%
    arrange(k, desc(cluster), sil_width) %>%
    mutate(kcr = str_c(k, cluster, rid, sep = '_') %>% as_factor()) %>%
    ggplot(aes(kcr, sil_width, fill = cluster, col = cluster)) +
    coord_flip() +
    geom_col() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    facet_wrap(~lab, scales = 'free_y') +
    ggtitle('Cluster silhouette plot')

  # # speed up EM a bit
  # read_derep <-
  #   tibble(read_id = seq_len(ncol(read_allele_flt))) %>%
  #   mutate(hash = map_chr(read_id, ~ digest::digest(read_allele_flt[, .]))) %>%
  #   group_by(hash) %>%
  #   summarise(rep_id = min(read_id),
  #             read_id = list(read_id),
  #             n = n()) %>%
  #   select(-hash)

  hap_em_search <-
    pam_hap_k$k_medoids %>%
    group_by(k) %>%
    summarise(
      label = str_c('k=', first(k), ', nc=', first(nc), ', ratio=', str_c(ratio, collapse = ':')),
      em_res = list(
        hap_mix_em(allele_mat = unname(read_allele_flt),
                   var_pos = var_pos[var_pass],
                   haps = do.call(cbind, medoid_state),
                   hap_prop = ratio %>% { . / sum(.) },
                   ts_rate = 1e-4,
                   ts_max = 3,
                   max_iter = 10,
                   fixed_hap_prop = TRUE,
                   epsilon = 1e-2))) %>%
    mutate(LH = map_dbl(em_res, ~last(.$LH)),
           label = str_c(label, ', LH=', format(LH)))

  ## final results
  # table w read_id, haplotype, status, missing_rate, error_rate, posterior
}

#' @importFrom rlang is_scalar_character
chimclustr_report <- function(filename,
                              plot_1_data,
                              plot_2_data,
                              sil_data,
                              read_allele_flt,
                              hap_em_search) {

  stopifnot(is_scalar_character(filename))

  filename <- `if`(str_ends(filename, '\\.html'), filename, str_c(filename, '.html'))
  message(str_c("saving report to ", filename))

  rmd_file <- system.file(file.path('Rmd', 'chimclustr_report.Rmd'), package = 'chimclsutr', mustWork = TRUE)
  rmd_env <- list2env(data, envir = new.env())

  rmarkdown::render(input = rmd_file,
                    output_file = filename,
                    envir = environment())
}

