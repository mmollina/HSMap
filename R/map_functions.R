#' Mapping functions (centiMorgan scale) and their inverses
#'
#' @description
#' Compute recombination fractions \eqn{r \in [0, 0.5)} from genetic distances
#' (in centiMorgans), and the inverse transforms (distance in cM from \eqn{r}).
#' Implementations are numerically stable, vectorized, and apply conservative
#' domain handling:
#'
#' - For functions that take \eqn{r} as input (`*_inv()`): values \eqn{r \le 0}
#'   are set to \eqn{0}; values \eqn{r \ge 0.5} are set to \eqn{0.5 - 1e-14}
#'   to avoid infinite distances.
#' - For functions that take distance `d` (in cM) as input: values \eqn{d < 0}
#'   are set to \eqn{0}.
#'
#' The three families covered here:
#'
#' \itemize{
#'   \item \strong{Haldane}: no interference.
#'     \deqn{r = \tfrac{1}{2}\left(1 - e^{-d/50}\right), \quad
#'           d = -50 \log(1 - 2r)}
#'
#'   \item \strong{Kosambi}: positive interference.
#'     \deqn{r = \frac{e^{d/25} - 1}{2\left(e^{d/25} + 1\right)}
#'           \;=\; \tfrac{1}{2}\tanh\!\left(\frac{d}{50}\right), \quad
#'           d = 25 \log\!\left(\frac{1 + 2r}{1 - 2r}\right)}
#'
#'   \item \strong{Morgan (linear)}: distance proportional to recombination.
#'     \deqn{r = d/100, \quad d = 100\,r}
#'     (This “Morgan” option is the linear/identity mapping in Morgans,
#'     provided for convenience alongside Haldane/Kosambi.)
#' }
#'
#' @section Units:
#' All distances `d` are in \strong{centiMorgans (cM)}; recombination fractions
#' `r` are unitless. If you prefer to work in Morgans, divide cM by 100.
#'
#' @param r Numeric vector of recombination fractions (expected in \eqn{[0, 0.5)}).
#' @param d Numeric vector of genetic distances in centiMorgans (cM).
#'
#' @return A numeric vector of the same length as the input, with `NA` preserved.
#'
#' @examples
#' ## Haldane:
#' inv_haldane(c(0.01, 0.1, 0.25))     # r -> d (cM)
#' haldane(     c(1, 10, 25))          # d (cM) -> r
#'
#' ## Kosambi:
#' inv_kosambi(c(0.01, 0.1, 0.25))
#' kosambi(     c(1, 10, 25))
#'
#' ## Morgan (linear):
#' inv_morgan(c(0.01, 0.1, 0.25))
#' morgan(     c(1, 10, 25))
#'
#' ## Composition checks (within numerical tolerance):
#' all.equal(haldane(inv_haldane(0.1)), 0.1, tol = 1e-12)
#' all.equal(kosambi(inv_kosambi(0.1)), 0.1, tol = 1e-12)
#' all.equal(morgan(inv_morgan(0.1)),   0.1, tol = 1e-12)
#'
#' @name map_functions_cm
NULL

# ---- Haldane (cM) ------------------------------------------------------------

#' @rdname map_functions_cm
#' @export
inv_haldane <- function(r) {
  r <- as.numeric(r)
  if (!length(r)) return(r)
  # domain guard on r
  r[is.na(r)] <- NA_real_
  r[r < 0]    <- 0
  r[r >= 0.5] <- 0.5 - 1e-14
  # d(cM) = -50 * log(1 - 2r)
  # use log1p for better stability when r ~ 0
  res <- -50 * log1p(-2 * r)
  res
}

#' @rdname map_functions_cm
#' @export
haldane <- function(d) {
  d <- as.numeric(d)
  if (!length(d)) return(d)
  # domain guard on d
  d[is.na(d)] <- NA_real_
  d[d < 0]    <- 0
  # r = 0.5 * (1 - exp(-d/50))
  # use expm1 for better stability when d ~ 0
  res <- 0.5 * (-expm1(-d / 50))
  # clamp to [0, 0.5)
  res[res < 0]    <- 0
  res[res >= 0.5] <- 0.5 - 1e-14
  res
}

# ---- Kosambi (cM) ------------------------------------------------------------

#' @rdname map_functions_cm
#' @export
inv_kosambi <- function(r) {
  r <- as.numeric(r)
  if (!length(r)) return(r)
  r[is.na(r)] <- NA_real_
  r[r < 0]    <- 0
  r[r >= 0.5] <- 0.5 - 1e-14
  # d(cM) = 25 * log((1 + 2r) / (1 - 2r))
  # use log1p for stability
  num <- log1p(2 * r)          # log(1 + 2r)
  den <- log1p(-2 * r)         # log(1 - 2r)
  res <- 25 * (num - den)
  res
}

#' @rdname map_functions_cm
#' @export
kosambi <- function(d) {
  d <- as.numeric(d)
  if (!length(d)) return(d)
  d[is.na(d)] <- NA_real_
  d[d < 0]    <- 0
  # r = 0.5 * tanh(d/50)
  # tanh(x) = (exp(2x) - 1) / (exp(2x) + 1); use expm1 for stability
  x   <- d / 50
  e2x <- exp(2 * x)
  res <- 0.5 * (e2x - 1) / (e2x + 1)
  # clamp to [0, 0.5)
  res[res < 0]    <- 0
  res[res >= 0.5] <- 0.5 - 1e-14
  res
}

# ---- Morgan (linear; cM) -----------------------------------------------------

#' @rdname map_functions_cm
#' @export
inv_morgan <- function(r) {
  r <- as.numeric(r)
  if (!length(r)) return(r)
  r[is.na(r)] <- NA_real_
  r[r < 0]    <- 0
  r[r >= 0.5] <- 0.5 - 1e-14
  # linear: d(cM) = 100 * r
  100 * r
}

#' @rdname map_functions_cm
#' @export
morgan <- function(d) {
  d <- as.numeric(d)
  if (!length(d)) return(d)
  d[is.na(d)] <- NA_real_
  d[d < 0]    <- 0
  # linear: r = d / 100
  res <- d / 100
  res[res < 0]    <- 0
  res[res >= 0.5] <- 0.5 - 1e-14
  res
}
