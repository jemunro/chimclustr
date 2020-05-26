
#' @importFrom dplyr first last mutate mutate_all select filter arrange desc
#' @importFrom magrittr '%>%' set_rownames set_colnames
#' @importFrom tibble tibble as_tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @importFrom purrr map pmap map_dbl map_chr map_int map_lgl reduce

chimclustr <- function(allele_matrix,
                       var_pos,
                       max_copy_num,
                       max_copy_num,
                       var_max_miss_rate,
                       read_max_miss_rate,
                       read_min_lh,
                       read_max_error_rate) {

  missingness <- marginal_rate_em(is.na(allele_matrix))
  read_pass <- which(missingness$col_rate < 0.05)
  var_pass <- which(missingness$row_rate < 0.25)

  row_annot <-
    tibble(rid = colnames(allele_matrix),
           status = if_else(seq_along(rid) %in% read_pass, 'pass', 'fail')) %>%
    as.data.frame() %>%
    column_to_rownames('rid')

  col_annot <-
    tibble(vid = rownames(allele_matrix),
           status = if_else(seq_along(vid) %in% var_pass, 'pass', 'fail')) %>%
    as.data.frame() %>%
    column_to_rownames('vid')

  read_dissim_1 <-
    t(allele_matrix) %>%
    as_tibble() %>%
    mutate_all(~ replace_na(., -1L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')

  read_hclust <- hclust(read_dissim_1, method = 'average')
  ## visualise pre filter
  pheatmap::pheatmap(t(allele_matrix),
                     cluster_rows = read_hclust,
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_row = row_annot,
                     annotation_col = col_annot,
                     annotation_colors = list(status = c(pass = 'chartreuse4', fail = 'darkorange4')),
                     color = viridisLite::cividis(2),
                     na_col = 'gray50',
                     legend = FALSE)

  read_allele_flt <- allele_matrix[var_pass, read_pass]




  pam_hap_k <- pam_hap_k(read_allele_flt, max_clust = 4)

  ## visualise clustering
  read_dissim_2 <-
    t(read_allele_flt) %>%
    as_tibble() %>%
    mutate_all(~ replace_na(., -1L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')
  read_hclust_2 <- hclust(read_dissim_2, method = 'average')

  row_clust_annot <-
    pam_hap_k$sil_data %>%
    select(k, cluster, rid) %>%
    unnest(rid) %>%
    pivot_wider(names_from = k, values_from = cluster, names_prefix = 'k_') %>%
    arrange(rid) %>%
    as.data.frame() %>%
    column_to_rownames('rid')

  row_clust_col <-
    row_clust_annot %>%
    map(unique) %>%
    map(~ set_names(scales::hue_pal()(length(.)), .))

  pheatmap::pheatmap(t(read_allele_flt) %>% set_rownames(seq_len(nrow(.))),
                     cluster_rows = read_hclust_2,
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     color = viridisLite::cividis(2),
                     annotation_row = row_clust_annot,
                     annotation_colors = row_clust_col,
                     na_col = 'gray50',
                     legend = FALSE)

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


  pam_k_asw <-
    pam_hap_k$sil_data %>%
    filter(k > 1) %>%
    select(k, k_asw, k_asw_filt) %>%
    distinct()

  pam_hap_k$sil_data %>%
    filter(k > 1) %>%
    select(k, cluster, sil_width, is_noise, k_asw, k_asw_filt) %>%
    unnest(sil_width) %>%
    ggplot(aes(cluster, sil_width)) +
    geom_boxplot(aes(fill = is_noise), col = 'gray25') +
    geom_hline(aes(linetype = 'max asw w/o noise', yintercept = max(pam_k_asw$k_asw_filt)), col = 'red') +
    geom_hline(data = pam_k_asw, aes(yintercept = k_asw, linetype = 'asw')) +
    geom_hline(data = pam_k_asw, aes(yintercept = k_asw/2, linetype = 'noise threshold'), col='darkblue') +
    geom_hline(data = pam_k_asw, aes(yintercept = k_asw_filt, linetype = 'asw w/o noise')) +
    scale_fill_manual(name = 'noise', values = c(`TRUE` = 'gray50', `FALSE` = 'gray80')) +
    scale_linetype_manual(name = NULL, values = c(1,2,1,1),
                          guide = guide_legend(override.aes = list(color = c("black", "black", "red", 'darkblue')))) +
    facet_grid(~k, scales = 'free_x', space = 'free_x', labeller = 'label_both') +
    ggtitle('Cluster silhouette width boxplot')

  # speed up EM a bit
  read_derep <-
    tibble(read_id = seq_len(ncol(read_allele_flt))) %>%
    mutate(hash = map_chr(read_id, ~ digest::digest(read_allele_flt[, .]))) %>%
    group_by(hash) %>%
    summarise(rep_id = min(read_id),
              read_id = list(read_id),
              n = n()) %>%
    select(-hash)

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

  em_plots <-
    hap_em_search %>%
    pmap(function(k, label, em_res, LH, ...) {
      p_lh <-
        em_res %>%
        select(iter, LH) %>%
        ggplot(aes(iter, LH)) +
        geom_line() +
        geom_point(size = 0.5) +
        ylab('log likelihood') +
        xlab('iteration')

      p_ve <-
        em_res %>%
        select(iter, rate = var_error_rate) %>%
        mutate(rate = map(rate, ~ tibble(rate = ., id = seq_along(.)))) %>%
        unnest(rate) %>%
        mutate(id = as_factor(id)) %>%
        ggplot(aes(iter, rate, group = id)) +
        geom_line(alpha = 0.25) +
        ylab('variant error rate') +
        xlab('iteration')

      p_re <-
        em_res %>%
        select(iter, rate = read_error_rate) %>%
        mutate(rate = map(rate, ~ tibble(rate = ., id = seq_along(.)))) %>%
        unnest(rate) %>%
        mutate(id = as_factor(id)) %>%
        ggplot(aes(iter, rate, group = id)) +
        geom_line(alpha = 0.10) +
        ylab('read error rate') +
        xlab('iteration')

      p_post <-
        last(em_res$read_hap_post) %>%
        group_by(read_id) %>%
        arrange(desc(post_lh)) %>%
        slice(1) %>%
        ungroup() %>%
        select(read_id, haplotype) %>%
        left_join(select(em_res, iter, read_hap_post) %>% unnest(read_hap_post),
                  by = c('read_id', 'haplotype')) %>%
        ggplot(aes(iter, post_lh, group = read_id)) +
        geom_line(alpha = 0.10) +
        ylab('posterior probability') +
        xlab('iteration') +
        facet_wrap(~haplotype, ncol = 1)

      p1 <- plot_grid(p_lh, p_ve, p_re, align = 'v', ncol = 1)
      p2 <- plot_grid(p1, p_post, ncol = 2)
      plot_grid(ggdraw() + draw_text(str_c('Haplotype Mixture EM Inspection (', label, ')')), p2,
                ncol = 1, rel_heights = c(1,15))
    })


  phase_plots <-
    hap_em_search %>%
    pmap(function(k, label, em_res, LH, ...) {
      read_data <-
        last(em_res$read_hap_post) %>%
        group_by(read_id) %>%
        arrange(desc(post_lh)) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(read_id) %>%
        mutate(error_rate = last(em_res$read_error_rate)) %>%
        select(read_id, haplotype, post_lh, error_rate) %>%
        mutate(status = if_else(error_rate < 0.05 & post_lh > 0.99 & haplotype != 'chimera', 'pass', 'fail')) %>%
        arrange(desc(haplotype), status, desc(error_rate), post_lh) %>%
        mutate(read_id = as_factor(as.character(read_id)))

      p1 <-
        unname(t(read_allele_flt)) %>%
        as_tibble() %>%
        mutate(read_id = seq_len(n())) %>%
        pivot_longer(starts_with('V'),
                     names_to = 'variant',
                     names_prefix = 'V',
                     names_ptypes = list(variant = integer()),
                     values_to = 'allele') %>%
        mutate(read_id = factor(read_id, levels(read_data$read_id)),
               allele = as.factor(allele),
               variant = as_factor(variant)) %>%
        left_join(select(read_data, haplotype, status, read_id), 'read_id') %>%
        ggplot(aes(variant, read_id, fill = allele)) +
        geom_tile() +
        scale_fill_manual(values = setNames(viridisLite::cividis(2), 1:2),
                          na.value = 'gray50') +
        facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
        theme(axis.ticks = element_blank(), axis.text = element_blank(),
              axis.title.y = element_blank(), legend.position = 'top',
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

      p2 <-
        read_data %>%
        ggplot(aes(x='x', read_id, fill = status)) +
        geom_tile() +
        theme() +
        facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
        theme(legend.position = 'top', strip.background = element_blank(), strip.text = element_blank(),
              axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
        xlab('status')

      p3 <-
        read_data %>%
        ggplot(aes(x='x', read_id, fill = error_rate)) +
        geom_tile() +
        theme() +
        scale_fill_viridis_c(limits = c(0,1)) +
        facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
        theme(legend.position = 'top', strip.background = element_blank(), strip.text = element_blank(),
              axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
        xlab('error')

      p4 <-
        read_data %>%
        rename(posterior = post_lh) %>%
        ggplot(aes(x='x', read_id, fill = posterior)) +
        geom_tile() +
        theme() +
        scale_fill_viridis_c(option = 'plasma', limits = c(0, 1)) +
        facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
        theme(legend.position = 'top', strip.background = element_blank(), strip.text = element_blank(),
              axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
        xlab('posterior')

      rm_legend <- theme(legend.position = 'none')

      title <- ggdraw() + draw_text(str_c('Haplotype Mixture EM Results (', label, ')'))
      pt <- plot_grid(get_legend(p1), get_legend(p2), get_legend(p3), get_legend(p4), nrow = 1)
      pb <- plot_grid(p1 + rm_legend, p2 + rm_legend, p3 + rm_legend, p4 + rm_legend, nrow = 1, rel_widths = c(10,1,1,1), align = 'h', axis = 'bt')
      plot_grid(title, pt, pb, ncol = 1, rel_heights = c(1, 1, 10))
    })
}



