suppressMessages(library(jsonlite))
ref <- fromJSON("ref.json", simplifyVector = FALSE)
f_ci <- function(A, B, imp, tol) {
  if (imp == "decrease") { A <- -A; B <- -B }
  m <- length(A); n <- length(B)
  Q <- outer(A, B, function(a, b) (b > a) + 0.5 * (b == a)); NAP <- mean(Q)
  h <- (m + n) / 2 - 1; z <- qnorm(.975)
  f <- function(x) m * n * (NAP - x)^2 * (2 - x) * (1 + x) - z^2 * x * (1 - x) * (2 + h + (1 + 2 * h) * x * (1 - x))
  c(if (NAP > 0) uniroot(f, c(0, NAP), tol = tol)$root else 0, if (NAP < 1) uniroot(f, c(NAP, 1), tol = tol)$root else 1)
}
out <- t(sapply(ref$pairs, function(p) f_ci(unlist(p$A), unlist(p$B), p$improvement, 1e-13)))
write_json(out, "ci_tight.json", digits = NA)
cat("ok\n")
