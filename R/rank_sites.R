#' Rank sites through backward elimination of redundancy
#'
#' This function uses a backward elimination algorithm to order a set of sites
#' by their environmental redundancy. A site is redundant if it adds little
#' environmental range, summed across species. Therefore, we define the most
#' redundant site as the site that, in whitened environmental space, has the
#' smallest distance to its nearest environmental neighbour, summed across
#' species (and evaluated independently for each species). This metric rewards
#' sites that are both environmentally distinct and speciose.
#'
#' Sites are removed one at a time, always dropping the most redundant site.
#' This gives a nested ranking, whereby the top `k` sites are near-opimal for
#' any given `k`.
#'
#' @param species Sites x species matrix of presence or abundance (we only
#'   consider presence, i.e. all values > 0 are equivalent). Species absent
#'   from all sites are ignored. Row names are the site ids.
#' @param env Sites x environment numeric matrix without `NA`, with the same
#'   row names as `species` (in any order). These need not be normalised or
#'   decorrelated, as all computation happens on a whitened subset of these
#'   variables (i.e. a PCA space explainging a proportion, `q`, of the overall
#'   variance)
#' @param weights Optional vector of site preferences, named by site, default 1. The
#'   redundancy of sites are multiplied by these weights, e.g. a site with
#'   weight 3 can be three times as redundant as a site with weight 1 before it
#'   is removed. Can be given for a subset of sites (e.g. `c(site2=3,
#'   site6=4)`); other values will be filled as 1.
#' @param min_sites Avoid droping a site if that would leave a species fewer
#'   than this many sites. Beyond some point, this will be impossible, so once
#'   no site can be dropped without this constrait it is ignored.
#' @param q Proportion of environmental variance to keep after whitening. 
#' @param bw_mult Kernel bandwidth as a multiple of the median
#'   nearest-neighbour distance. Sites closer than this are treated as the same
#'   environment.
#'
#' @return A data frame with one row per site: `site`, `rank` (1 = least
#'   redundant), `loss`, and `weight`.
#' @export
#' @importFrom stats setNames
rank_sites = function(species, env, weights = NULL, min_sites = 0, q = 0.9, bw_mult = 1) {
  species = as.matrix(species)
  env = as.matrix(env)
  ids = rownames(species)
  check_inputs(species, env, min_sites, q)
  env = env[ids, , drop = FALSE]
  weights = site_weights(weights, ids)

  species = species[, colSums(species > 0) > 0, drop = FALSE]
  if (ncol(species) == 0) stop("no species are held by any site", call. = FALSE)
  presence = (species > 0) * 1

  kernel = similarity(whiten(env, q), bw_mult)
  res = backward_order(kernel, presence, min_sites, weights)
  if (!is.na(res$floor_released_at)) {
    warning("min_sites cannot be kept below k = ", res$floor_released_at,
            "; ranks beyond that ignore it", call. = FALSE)
  }

  # last dropped = collected first
  keep = rev(res$order)
  data.frame(site = ids[keep], rank = seq_along(ids), loss = rev(res$loss),
             weight = weights[keep], row.names = NULL, stringsAsFactors = FALSE)
}

check_inputs = function(species, env, min_sites, q) {
  ids = rownames(species)
  if (is.null(ids) || anyDuplicated(ids)) {
    stop("`species` needs unique row names", call. = FALSE)
  }
  if (!setequal(ids, rownames(env)) || anyDuplicated(rownames(env))) {
    stop("`species` and `env` must have the same sites (row names)", call. = FALSE)
  }
  if (!is.numeric(env) || anyNA(env)) {
    stop("`env` must be numeric without missing values", call. = FALSE)
  }
  if (anyNA(species)) stop("`species` must not contain missing values", call. = FALSE)
  if (!(is.numeric(min_sites) && length(min_sites) == 1 && min_sites >= 0)) {
    stop("`min_sites` must be a single number >= 0", call. = FALSE)
  }
  if (!(is.numeric(q) && length(q) == 1 && q > 0 && q <= 1)) {
    stop("`q` must be a single number in (0, 1]", call. = FALSE)
  }
}

# Sites without a weight keep the default of 1
site_weights = function(weights, ids) {
  out = setNames(rep(1, length(ids)), ids)
  if (is.null(weights)) return(out)
  if (is.null(names(weights))) stop("`weights` must be named by site", call. = FALSE)
  unknown = setdiff(names(weights), ids)
  if (length(unknown)) {
    stop("`weights` names unknown sites, e.g. ", unknown[[1]], call. = FALSE)
  }
  out[names(weights)] = weights
  if (!all(is.finite(out)) || any(out <= 0)) {
    stop("`weights` must be finite and positive", call. = FALSE)
  }
  out
}
