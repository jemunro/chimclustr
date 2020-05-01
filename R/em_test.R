

gm_sample <- function(n = 1000, mu = c(0, 2), prob = c(0.4, 0.6)) {
  state <- sample(mu, n, replace = T, prob = prob)
  tibble(mu = state,
         x = rnorm(n, state))
}

gmm_em <- function(x, k = 2, max_iter = 1000L, weighted = TRUE) {
  # init with kmeans
  cl <- kmeans(x, centers = k)
  state <- cl$cluster
  mu <- c(cl$centers)
  sigma <- map_dbl(seq_len(k), function(i) sd(x[state == i]))
  prop <- map_dbl(seq_len(k), function(i) sum(state == i)/ length(x))

  # EM phase
  n_iter <- 0
  res <- tibble(iter = n_iter,
                mu = list(mu),
                sigma = list(sigma),
                prop = list(prop),
                cluster = list(seq_len(k)))
  while (n_iter < max_iter) {
    n_iter <- n_iter + 1
    # E
    ps <- expectation(x, prop, mu, sigma)
    # M
    mx <- maximise(x, ps, weighted = weighted)
    mu <- mx$mu
    sigma <- mx$sigma
    prop <- colSums(ps) / length(x)
    res <- add_row(res,
                   iter = n_iter,
                   mu = list(mu),
                   sigma = list(sigma),
                   prop = list(prop),
                   cluster = list(seq_len(k)))
  }
  return(unnest_legacy(res) %>%
           mutate(cluster = as.factor(cluster)))
}

expectation <- function(x, prop, mu, sigma) {
  lh <- vapply(seq_along(mu),
               function(i) { prop[i] * dnorm(x, mu[i], sigma[i], log = F)},
               double(length(x)))
  lh / rowSums(lh)
}

maximise <- function(x, ps, weighted = TRUE) {
  # probability weighted params - why
  if (weighted) {
    mu <- map_dbl(seq_len(ncol(ps)),
                  function(i) weighted.mean(x, ps[, i]))
    sigma <- map_dbl(seq_len(ncol(ps)),
                     function(i) radiant.data::weighted.sd(x, ps[, i] / sum(ps[, i])))
  } else {
    state <- map_int(seq_len(nrow(ps)), ~ which.max(ps[., ]))
    mu <- map_dbl(seq_len(ncol(ps)), function(i) mean(x[state == i]))
    sigma <- map_dbl(seq_len(ncol(ps)), function(i) sd(x[state == i]))
  }
  return(list(mu = mu, sigma = sigma))
}


sm <- gm_sample()
x <- sm$x
p <- table(sm$mu) / length(x)


x <- c(rnorm(500, mean = 0), rnorm(500, mean = 2))
res_w <- gmm_em(x)

res_w %>%
  gather(-iter, -cluster, key = 'param', value = 'value') %>%
  ggplot(aes(iter, value, col = cluster)) +
  geom_line() +
  facet_wrap(~ param, scales = 'free', nrow = 1)

# unweighted version does not work
# res_uw <- gmm_em(x, weighted = F)
# res_uw %>%
#   gather(-iter, -cluster, key = 'param', value = 'value') %>%
#   ggplot(aes(iter, value, col = cluster)) +
#   geom_line() +
#   facet_wrap(~ param, scales = 'free', nrow = 1)

