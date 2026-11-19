## Zero-truncated negative binomial distribution
## d, p, q and r functions for fitdistrplus (used for culling events per LGA).
## Return NaN for invalid parameters so the fitting algorithm moves away from them.

library(truncdist)

dtrunc_nbinom <- function(x, size, mu) {
  tryCatch({
    if (mu < 1 | size < 0) stop()
    dtrunc(x, "nbinom", a = 0, b = Inf, size = size, mu = mu)
  }, error = function(e) rep(NaN, length(x)))
}

ptrunc_nbinom <- function(q, size, mu) {
  tryCatch({
    if (mu < 1 | size < 0) stop()
    ptrunc(q, "nbinom", a = 0, b = Inf, size = size, mu = mu)
  }, error = function(e) rep(NaN, length(q)))
}

qtrunc_nbinom <- function(p, size, mu) {
  tryCatch({
    if (mu < 1 | size < 0) stop()
    qtrunc(p, "nbinom", a = 0, b = Inf, size = size, mu = mu)
  }, error = function(e) rep(NaN, length(p)))
}

rtrunc_nbinom <- function(n, size, mu) {
  rtrunc(n, "nbinom", a = 0, b = Inf, size = size, mu = mu)
}
