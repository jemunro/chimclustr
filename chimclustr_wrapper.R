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
  --max-copy-num=<f>          Maximum copy number of target [default: 2].
  --min-copy-num=<f>          Minimum copy number of target [default: 1].
  --var-max-miss-rate=<f>     Maximum missing rate per variant frequency for pre-filte [default: 0.25].
  --read-max-miss-rate=<f>    Maximum rate of missing alleles per read for pre-filter [default: 0.05].
  --read-min-lh=<f>           Minimum posterior likelihood of a read assigned to phase [default: 0.99].
  --read-max-error-rate=<f>   Maximum expected error rate of a read assigned to phase [default: 0.10].
  --em-max-iter=<f>           Maximum number of EM iterations [default: 10].
  --em-ts-max=<f>             Maximum number of template switches considered in em [default: 2].
" -> doc

opts <-
  docopt(doc, c('/stornext/HPCScratch/home/munro.j/runs/pba/test-gt/work/73/e42f60e2a6a6a096899456c3f08723/SM-AR-003.haplotag.tsv.gz',
                '/stornext/HPCScratch/home/munro.j/runs/pba/test-gt/progress/gatk_jdc_1/snp_pos.gz',
                'output'))

max_copy_num <- as.integer(opts$`max-copy-num`)
min_copy_num <- as.integer(opts$`min-copy-num`)
var_max_miss_rate <- as.numeric(opts$`var-max-miss-rate`)
read_max_miss_rate <- as.numeric(opts$`read-max-miss-rate`)
read_min_lh <- as.numeric(opts$`read-min-lh`)
read_max_error_rate <- as.numeric(opts$`read-max-error-rate`)
em_max_iter <- as.integer(opts$`em-max-iter`)
em_ts_max <- as.integer(opts$`em-ts-max`)
out_prefix = opts$out_prefix

stopifnot(file.exists(opts$read_vars_tsv),
          file.exists(opts$var_pos))

var_pos <- scan(opts$var_pos, what = integer(), quiet = TRUE)

allele_matrix <-
  read_tsv(opts$read_vars_tsv,
           col_types = cols(PS = col_integer(), HP = col_integer())) %>%
  mutate(vid = as.integer(PS + 1)) %>%
  select(rid = name, vid, allele = HP) %>%
  pivot_wider(names_from = rid, values_from = allele) %>%
  arrange(vid) %>%
  as.data.frame() %>%
  column_to_rownames('vid') %>%
  as.matrix()

result <- chimclustr(allele_matrix = allele_matrix,
                     var_pos = var_pos,
                     max_copy_num = max_copy_num,
                     var_max_miss_rate = var_max_miss_rate,
                     read_max_miss_rate = read_max_miss_rate,
                     read_min_lh = read_min_lh,
                     read_max_error_rate = read_max_error_rate,
                     em_max_iter = em_max_iter,
                     em_ts_max = em_ts_max,
                     report_filename = out_prefix)

write_tsv(result, str_c(out_prefix, '.tsv.gz'))

