# Whiten/decorrelate likely-correlated environmental variables. This makes
# Euclidean distance equal Mahalanobis distance in the original variables, so
# correlated variables (e.g. temperature bioclims) count once instead of
# once each.
whiten = function(env, q) {
  varying = apply(env, 2, stats::sd) > 0
  pca = stats::prcomp(env[, varying, drop = FALSE], center = TRUE, scale. = TRUE)
  explained = cumsum(pca$sdev^2) / sum(pca$sdev^2)
  usable = sum(pca$sdev > .Machine$double.eps * pca$sdev[1])
  k = min(which(explained >= q)[1], usable)
  sweep(pca$x[, seq_len(k), drop = FALSE], 2, pca$sdev[seq_len(k)], "/")
}

# Gaussian similarity, ~1 within a bandwidth of each other and ~0 far apart. The
# median nearest-neighbour distance is a natural yardstick for "same
# environment".
similarity = function(z, bw_mult) {
  d = as.matrix(stats::dist(z))
  nearest = d
  diag(nearest) = Inf
  sigma = bw_mult * stats::median(apply(nearest, 1, min))
  if (!is.finite(sigma) || sigma <= 0) sigma = max(d) / 10
  exp(-d^2 / (2 * sigma^2))
}


# Greedy backward elimination on the species-additive coverage objective.
backward_order = function(kernel, presence, min_sites, weights) {
  n = nrow(kernel)
  held = presence > 0
  n_held = colSums(presence)
  # mean over a species' own sites, so every species weighs the same
  species_weight = 1 / n_held

  best = matrix(0, n, ncol(presence))
  second = matrix(0, n, ncol(presence))
  best_site = matrix(0L, n, ncol(presence))
  second_site = matrix(0L, n, ncol(presence))
  for (s in seq_len(ncol(presence))) {
    holders = which(held[, s])
    top = top_two(kernel[holders, holders, drop = FALSE], holders)
    best[holders, s] = top$best
    best_site[holders, s] = top$best_site
    second[holders, s] = top$second
    second_site[holders, s] = top$second_site
  }

  gap = (best - second) * rep(species_weight, each = n) * held
  loss = add_by_site(numeric(n), as.vector(best_site), as.vector(gap))

  retained = rep(TRUE, n)
  order = integer(n)
  loss_at_removal = numeric(n)
  floor_released_at = NA_integer_

  for (step in seq_len(n)) {
    at_floor = n_held <= min_sites
    blocked = if (any(at_floor)) as.vector(presence %*% at_floor) > 0 else rep(FALSE, n)
    droppable = retained & !blocked
    if (!any(droppable)) {
      # The order must still be total, so give the floor up and note where
      droppable = retained
      if (is.na(floor_released_at)) floor_released_at = n - step + 1L
    }
    p = which.min(ifelse(droppable, loss * weights, Inf))
    order[step] = p
    loss_at_removal[step] = loss[p]
    retained[p] = FALSE

    for (s in which(held[p, ])) {
      rows = which(held[, s] & (best_site[, s] == p | second_site[, s] == p))
      if (length(rows) == 0) next
      loss = add_by_site(loss, best_site[rows, s], -gap[rows, s])
      holders = which(held[, s] & retained)
      top = top_two(kernel[rows, holders, drop = FALSE], holders)
      best_site[rows, s] = top$best_site
      second_site[rows, s] = top$second_site
      gap[rows, s] = (top$best - top$second) * species_weight[s]
      loss = add_by_site(loss, best_site[rows, s], gap[rows, s])
    }
    n_held = n_held - presence[p, ]
    loss[p] = 0
  }
  list(order = order, loss = loss_at_removal, floor_released_at = floor_released_at)
}

# Best and second-best column per row, as site ids (`cols`) and similarities.
# A missing runner-up (or no column at all) has similarity 0 and site 0.
top_two = function(sim, cols) {
  n = nrow(sim)
  none = list(best = numeric(n), best_site = integer(n),
               second = numeric(n), second_site = integer(n))
  if (length(cols) == 0) return(none)
  rows = seq_len(n)
  j1 = max.col(sim, ties.method = "first")
  best = sim[cbind(rows, j1)]
  if (ncol(sim) == 1) {
    return(list(best = best, best_site = cols[j1], second = numeric(n), second_site = integer(n)))
  }
  sim[cbind(rows, j1)] = -Inf
  j2 = max.col(sim, ties.method = "first")
  list(best = best, best_site = cols[j1],
       second = sim[cbind(rows, j2)], second_site = cols[j2])
}

# Add `value` into the elements of `total` named by site index; 0 means no site.
add_by_site = function(total, site, value) {
  used = site > 0
  if (any(used)) {
    sums = rowsum(value[used], site[used], reorder = FALSE)
    total[as.integer(rownames(sums))] = total[as.integer(rownames(sums))] + sums[, 1]
  }
  total
}
