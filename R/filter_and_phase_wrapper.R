#!/usr/bin/env Rscript

stopifnot(require(rlang),
          require(tidyverse),
          require(magrittr),
          require(docopt),
          require(cowplot),
          require(chimclustr))
"
Usage:
  filter_and_phase.R <read_vars_tsv> <var_pos> <out_prefix> [options]

Options:
  read_vars_tsv               Tsv file of qname, variant_id, genotype.
  var_pos                     File of variant positions, one per line.
  out_prefix                  Prefix of output files.
  --max-copy-num=<f>          Maximum copy number of target
  --min-copy-num=<f>          Minimum copy number of target
  --var-max-maj-af=<f>        Maximum major allele frequency for pre-filter [default: 0.95].
  --var-max-miss-rate=<f>     Maximum missing rate per variant frequency for pre-filte[default: 0.25].
  --read-max-miss-rate=<f>    Maximum rate of missing alleles per read for pre-filter [default: 0.20].
  --read-min-lh=<f>           Minimum posterior likelihood of a read assigned to phase [default: 0.75].
  --read-max-error-rate=<f>   Maximum expected error rate of a read assigned to phase [default: 0.25].
" -> doc

opts <- 
  docopt(doc, c('/stornext/HPCScratch/home/munro.j/runs/pba/test-gt/work/73/e42f60e2a6a6a096899456c3f08723/SM-AR-003.haplotag.tsv.gz',
                '/stornext/HPCScratch/home/munro.j/runs/pba/test-gt/progress/gatk_jdc_1/snp_pos.gz', 
                'out_tsv'))

###################### Chimclustr to do #########################
# Include all sites? - better estimates of rates, more consistent across models
#   - probably not too computationally expensive, as don't need to consider chimeras at most sites
# Start with better read error rate estimate – this takes a long time to converge currently (implicitly var error rate will be zero)
# Ultimately missingness and errors could be modelled separately
# Silhouette plots for tested clustering
# Heuristic noise detection and exclusion:
#   min_prop (=0.15) and min_rel_sil (=0.5)
#   Up to k + 2
# Is it useful to cap error rate? Does this also slow down convergence?


# set and check params
max_copy_num <- as.numeric(opts$`max-copy-num`)
max_copy_num <- as.numeric(opts$`min-copy-num`)
var_max_maj_af <- as.numeric(opts$`var-max-maj-af`)
var_max_miss_rate <- as.numeric(opts$`var-max-miss-rate`)
read_max_miss_rate <- as.numeric(opts$`read-max-miss-rate`)
read_min_lh <- as.numeric(opts$`read-min-lh`)
read_max_error_rate <- as.numeric(opts$`read-max-error-rate`)

stopifnot(file.exists(opts$read_vars_tsv),
          file.exists(opts$var_pos))

var_pos <- scan(opts$var_pos, what = integer(), quiet = TRUE)

read_allele_mat <-
  read_tsv(opts$read_vars_tsv,
         col_types = cols(PS = col_integer(), HP = col_integer())) %>%
  mutate(vid = as.integer(PS + 1)) %>%
  select(rid = name, vid, allele = HP) %>% 
  pivot_wider(names_from = rid, values_from = allele) %>% 
  arrange(vid) %>% 
  as.data.frame() %>% 
  column_to_rownames('vid') %>% 
  as.matrix()

## filter based on missingess
missingness <- marginal_rate_em(is.na(read_allele_mat))
read_pass <- which(missingness$col_rate < 0.05)
var_pass <- which(missingness$row_rate < 0.25)

row_annot <- 
  tibble(rid = colnames(read_allele_mat),
         status = if_else(seq_along(rid) %in% read_pass, 'pass', 'fail')) %>% 
  as.data.frame() %>% 
  column_to_rownames('rid')

col_annot <-
  tibble(vid = rownames(read_allele_mat),
         status = if_else(seq_along(vid) %in% var_pass, 'pass', 'fail')) %>% 
  as.data.frame() %>% 
  column_to_rownames('vid')

read_dissim_1 <-
  t(read_allele_mat) %>%
  as_tibble() %>%
  mutate_all(~ replace_na(., -1L)) %>% 
  mutate_all(as.factor) %>%
  cluster::daisy(metric = 'gower')

read_hclust <- hclust(read_dissim_1, method = 'average')
## visualise pre filter
pheatmap::pheatmap(t(read_allele_mat),
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

read_allele_flt <- read_allele_mat[var_pass, read_pass]




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
  summarise(em_res = list(
    hap_mix_em(allele_mat = unname(read_allele_flt),
                var_pos = var_pos[var_pass],
                haps = do.call(cbind, medoid_state),
                hap_prop = ratio %>% { . / sum(.) },
                ts_rate = 1e-4,
                ts_max = 2,
                max_iter = 10,
                fixed_hap_prop = TRUE,
                epsilon = 1e-2)
  )) %>% 
  mutate(LH = map_dbl(em_res, ~last(.$search$LH)))

plots <- 
  hap_em_search %>% 
  arrange(desc(LH)) %>% 
  pmap(function(k, em_res, LH) {
    
    read_data <-
      em_res$read_hap_lh_smry %>% 
      mutate(error_rate = last(em_res$search$read_error_rate)) %>% 
      select(read_id, haplotype, post_lh, error_rate) %>% 
      mutate(status = if_else(error_rate < 0.10 & post_lh > 0.95 & haplotype != 'chimera', 'pass', 'fail')) %>% 
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
    
    p2a <-
      read_data %>% 
      ggplot(aes(x='x', read_id, fill = status)) + 
      geom_tile() +
      theme() +
      facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
      theme(legend.position = 'top', strip.background = element_blank(), strip.text = element_blank(),
            axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
      xlab('status')
    
    p2 <-
      read_data %>% 
      ggplot(aes(x='x', read_id, fill = error_rate)) + 
      geom_tile() +
      theme() +
      scale_fill_viridis_c(limits = c(0,1)) +
      facet_grid(haplotype ~ ., space='free_y', scales = 'free_y', switch = 'both') +
      theme(legend.position = 'top', strip.background = element_blank(), strip.text = element_blank(),
            axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
      xlab('error_rate')
    
    p3 <-
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
    
    pt <- plot_grid(get_legend(p1), get_legend(p2a), get_legend(p2), get_legend(p3), nrow = 1)
    pb <- plot_grid(p1 + rm_legend, p2a+ rm_legend, p2+ rm_legend, p3+ rm_legend, nrow = 1, rel_widths = c(10,1,1,1), align = 'h', axis = 'bt')
    plot_grid(pt, pb, ncol = 1, rel_heights = c(1, 10))
    # return(plot_grid(pt, pb, ncol = 1, rel_heights = c(1, 10)))
  })
