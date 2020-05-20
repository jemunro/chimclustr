# obs. vars. (y) = read variant alleles
# lat. vars. (z) = haplotype at each variant
# parameters     = p_err, p_chim, p_hap_i

# require(tidyverse)
# require(magrittr)
# require(furrr)
# options(future.fork.enable = TRUE)
# options(future.globals.maxSize = Inf)
# plan(multiprocess, workers = 8)


run_test <- function() {

  data <- test_data()
  read_allele_mat <- test_matrix(data) %>% unname()
  max_ploidy <- 4L
  max_iter <- 10L

  flt <- filter_read_var(read_allele_mat,
                         max_na_read = 0.25,
                         max_na_var = 0.25,
                         max_maj_af = 0.90)

  read_allele_flt <- read_allele_mat[flt$var_set, flt$read_set, drop = FALSE]
  var_pos <- data$pos[flt$var_set]

  res <- chimclustr(read_allele_mat = read_allele_flt,
                    var_pos = var_pos,
                    max_ploidy = max_ploidy,
                    max_iter = max_iter)

  saveRDS(res, 'res.rds')

  # plots
  # 1) pre-filter heatmap ( + read/var NA hist and cutoff)
  # 2) chimclust results
  # 3) post-filter heatmap, split by phase

  # plotting
  # read_sample <- sample(ro, 500) %>% as.character()
  #
  # read_heatmap <-
  #   read_dissim_mat %>%
  #   as.data.frame() %>%
  #   rownames_to_column('row_id') %>%
  #   as_tibble() %>%
  #   gather(-row_id, key = 'col_id', value = 'sim') %>%
  #   mutate(sim = 1 - sim) %>%
  #   filter(row_id %in% read_sample, col_id %in% read_sample) %>%
  #   mutate(col_id = factor(col_id, levels = as.character(ro)),
  #          row_id = factor(row_id, levels = as.character(rev(ro)))) %>%
  #   ggplot() +
  #   geom_tile(aes(row_id, col_id, fill = sim)) +
  #   scale_fill_viridis_c(option = 'E') +
  #   theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
  #         axis.ticks.x = element_blank(), legend.position = 'top') +
  #   ylab('read') +
  #   xlab('read')
  # # read_heatmap
  #
  # ra_plot <-
  #   read_allele_df %>%
  #   { set_colnames(., 1:ncol(.)) } %>%
  #   as_tibble() %>%
  #   mutate_all(as.integer) %>%
  #   mutate(read_id = as.character(seq_len(n()))) %>%
  #   filter(read_id %in% read_sample) %>%
  #   gather(-read_id, key = 'var_id', value = 'allele') %>%
  #   mutate(read_id = factor(read_id, levels = as.character(ro)),
  #          var_id = factor(var_id, levels = as.character(rev(vo)))) %>%
  #   mutate(allele = as.character(allele)) %>%
  #   ggplot() +
  #   geom_tile(aes(var_id, read_id, fill = allele)) +
  #   theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  #   guides(fill = F) +
  #   xlab('variant') +
  #   scale_fill_manual(values =  c(`1` = "#FFEA46FF", `2`= "#00204DFF"))
  #
  # for( i in (seq_along(model_search$em_res))) {
  #   hap_plot <-
  #     model_search$em_res[[i]]$read_hap_lh_smry %>%
  #     mutate(read_id = as.character(read_id)) %>%
  #     filter(read_id %in% read_sample) %>%
  #     mutate(read_id = factor(read_id, ro)) %>%
  #     ggplot() +
  #     geom_tile(aes(x=1, y=read_id, fill = haplotype)) +
  #     theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
  #           axis.title.y = element_blank(), legend.position = 'top',
  #           axis.ticks.x = element_blank()) +
  #     scale_fill_brewer(palette = 'Dark2') +
  #     xlab('haplotype')
  #
  #   hap_lh_plot <-
  #     model_search$em_res[[i]]$read_hap_lh %>%
  #     mutate(read_id = as.character(read_id)) %>%
  #     filter(read_id %in% read_sample) %>%
  #     mutate(read_id = factor(read_id, ro)) %>%
  #     ggplot() +
  #     geom_tile(aes(x=haplotype, y=read_id, fill = raw_llh)) +
  #     theme(axis.text.y = element_blank(),
  #           axis.title.y = element_blank(),
  #           axis.ticks.x = element_blank(),
  #           axis.text.x = element_text(angle = 90)) +
  #     scale_fill_viridis_c() +
  #     guides(fill = guide_colorbar(title = 'log likelihood')) +
  #     ylab('read') +
  #     xlab('haplotype')
  #
  #
  #   plot_grid(read_heatmap, ra_plot, hap_plot, hap_lh_plot,
  #             rel_widths = c(15, 5, 3, 6),
  #             nrow = 1, align = 'h', axis = 'tblr')
  #
  #   ggsave(str_c(i, '.png'), width = 11, height = 8)
  # }

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
    select(vid, read_data) %>%
    unnest(read_data) %>%
    spread(rid, allele) %>%
    arrange(vid) %>%
    { as.matrix(select(., -vid)) }
}


chimclustr <- function(read_allele_mat,
                       var_pos,
                       max_ploidy = 4L,
                       max_iter = 10) {

  stopifnot(is.matrix(read_allele_mat),
            is.integer(read_allele_mat),
            all(is.na(read_allele_mat) | read_allele_mat > 0L),
            is_integerish(max_ploidy),
            is_integerish(max_iter),
            max_iter > 0,
            max_ploidy >= 1,
            is_integerish(var_pos),
            !any(is.na(var_pos)),
            length(var_pos) == nrow(read_allele_mat))

  read_allele_mat <- unname(read_allele_mat)

  read_derep <-
    tibble(read_index = seq_len(ncol(read_allele_mat)),
           hash = map_chr(read_index, ~ digest(read_allele_mat[, .]))) %>%
    group_by(hash) %>%
    summarise(rep_index = min(read_index),
              read_index = list(read_index),
              n = n()) %>%
    select(-hash) %>%
    arrange(rep_index)

  read_dissim <-
    t(read_allele_mat) %>%
    as_tibble() %>%
    mutate_all(function(x) replace_na(x, 0L)) %>%
    mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower')

  # this method can't choose k = 1
  # will need to compare to k = 1 in all cases
  pam_best_k <- fpc::pamk(read_dissim, krange = seq_len(max_ploidy))

  model_search <-
    tibble(k = 1,
           medoid = list(as.integer(cluster::pam(read_dissim, k=1)$medoids)),
           prop = list(1)) %>%
    add_row(k = pam_best_k$nc,
            medoid = list(as.integer(pam_best_k$pamobject$medoids)),
            prop = list(c(table(pam_best_k$pamobject$clustering)) %>% { . / sum(.) } )) %>%
    mutate(ratios = map(prop, ~ hap_ratios(., max_ploidy = max_ploidy))) %>%
    select(-prop) %>%
    unnest(ratios) %>%
    mutate(em_res = pmap(., function(medoid, ratio, ...) {
      message('ratio: ', str_c(ratio, collapse = ':'))
      chimclustr:::hap_mix_em(read_allele_mat = read_allele_mat[, read_derep$rep_index],
                              var_pos = var_pos,
                              haps = read_allele_mat[, medoid, drop = FALSE],
                              hap_prop = ratio / sum(ratio),
                              ts_rate = 1e-4,
                              ts_max = 2,
                              var_error_rate = 1e-2,
                              read_error_rate = 1e-2,
                              max_iter = max_iter,
                              epsilon = 0.1,
                              fixed_hap_prop = TRUE,
                              read_weight = read_derep$n)
    })) %>%
    mutate(LH = map_dbl(em_res, ~ last(.$search$LH)))

  return(model_search)
}




