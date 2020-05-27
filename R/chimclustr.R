
#' @importFrom dplyr first last mutate mutate_all select filter arrange desc
#' @importFrom magrittr '%>%' '%<>%' set_rownames set_colnames
#' @importFrom tibble tibble as_tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @importFrom purrr map pmap map_dbl map_chr map_int map_lgl reduce
#' @importFrom rlang is_integerish is_scalar_double
#' @export
chimclustr <- function(allele_matrix,
                       var_pos,
                       max_copy_num,
                       var_max_miss_rate,
                       read_max_miss_rate,
                       read_min_lh,
                       read_max_error_rate,
                       em_max_iter = 25L,
                       em_ts_max = 3L,
                       report_filename = NULL) {

  stopifnot(is.matrix(allele_matrix),
            is.integer(allele_matrix),
            is_integerish(var_pos),
            length(var_pos) == nrow(allele_matrix),
            is_scalar_double(var_max_miss_rate),
            is_scalar_double(read_max_miss_rate),
            is_scalar_double(read_min_lh),
            is_scalar_double(read_max_error_rate),
            is.null(report_filename) || is_scalar_character(report_filename))

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
  pam_hap_k <- pam_hap_k(read_allele_flt,
                         max_clust = max_copy_num,
                         min_ploidy = min_copy_num)

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
                   ts_max = em_ts_max,
                   max_iter = em_max_iter,
                   fixed_hap_prop = TRUE,
                   epsilon = 1e-2))) %>%
    mutate(LH = map_dbl(em_res, ~last(.$LH)),
           label = str_c(label, ', LH=', format(LH)))

  if (!is.null(report_filename)) {
    chimclustr_report(filename = report_filename,
                      plot_1_data = plot_1_data,
                      plot_2_data = plot_2_data,
                      sil_data = pam_hap_k$sil_data,
                      read_allele_flt = read_allele_flt,
                      hap_em_search = hap_em_search)
  }

  ## final results
  # table w read_id, haplotype, status, missing_rate, error_rate, posterior
  result <-
    hap_em_search %>%
    with(em_res[[which.max(LH)]]) %>%
    with(last(read_hap_post) %>%
           group_by(read_id) %>%
           slice(which.max(post_lh)) %>%
           ungroup() %>%
           arrange(read_id) %>%
           mutate(error_rate = last(read_error_rate))) %>%
    full_join(tibble(read_name = colnames(allele_matrix),
                     pass = seq_along(read_name) %in% read_pass,
                     read_id = if_else(pass, cumsum(pass), NA_integer_),
                     missing_rate = missingness$col_rate) %>%
                select(-pass),
              by = 'read_id') %>%
    select(read_name, haplotype, posterior = post_lh, error_rate, missing_rate) %>%
    mutate(status = if_else(posterior >= read_min_lh & error_rate <= read_max_error_rate & haplotype != 'chimera',
                            'pass', 'fail', missing = 'fail')) %>%
    arrange(read_name)

  return(result)
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
  if (!file.exists(filename)) {
    # creating file forces normalizePath to return full path when a relative filename is used
    invisible(file.create(filename))
  }
  filename <- normalizePath(filename)

  rmd_file <- system.file(file.path('Rmd', 'chimclustr_report.Rmd'), package = 'chimclustr', mustWork = TRUE)
  rmd_file_tmp <- file.path(dirname(filename), str_c('.', basename(rmd_file)))
  invisible(file.copy(rmd_file, rmd_file_tmp, overwrite = TRUE))

  rmarkdown::render(input = rmd_file_tmp,
                    output_file = filename,
                    envir = environment())
  message(str_c("Chimclustr report saved to ", filename))
  invisible(file.remove(rmd_file_tmp))
}

