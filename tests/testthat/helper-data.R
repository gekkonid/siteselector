# Small random problem with named sites, used by several tests
toy_problem = function(n_sites = 12, n_species = 4, n_env = 3, seed = 1) {
  set.seed(seed)
  ids = sprintf("s%02d", seq_len(n_sites))
  env = matrix(rnorm(n_sites * n_env), n_sites, dimnames = list(ids, NULL))
  species = matrix(rbinom(n_sites * n_species, 1, 0.6), n_sites, dimnames = list(ids, NULL))
  # every species needs a holder, or it is silently dropped
  species[cbind(sample(n_sites, n_species), seq_len(n_species))] = 1
  list(species = species, env = env)
}
