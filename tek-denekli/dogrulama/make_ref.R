# V-Able SCED motoru için referans değerler: SingleCaseES, scdhlm, scan
suppressMessages({library(SingleCaseES); library(scdhlm); library(jsonlite); library(nlme)})
set.seed(20261108)

## 1) Rastgele A-B çiftleri -------------------------------------------------
pairs <- vector("list", 300)
for (k in seq_along(pairs)) {
  m <- sample(3:12, 1); n <- sample(3:15, 1)
  type <- sample(c("count", "cont", "ties", "zeros"), 1)
  A <- switch(type,
    count = rpois(m, runif(1, 2, 10)),
    cont  = round(rnorm(m, 20, 5), 2),
    ties  = sample(1:4, m, TRUE),
    zeros = rpois(m, 0.6))
  shift <- runif(1, -3, 8)
  B <- switch(type,
    count = rpois(n, max(0.5, mean(A) + shift)),
    cont  = round(rnorm(n, 20 + shift, 5), 2),
    ties  = sample(2:6, n, TRUE),
    zeros = rpois(n, runif(1, 0.2, 4)))
  if (var(A) == 0) A[1] <- A[1] + 1
  imp <- sample(c("increase", "decrease"), 1)
  es <- function(f, ...) tryCatch(f(A_data = A, B_data = B, improvement = imp, ...), error = function(e) NULL)
  r <- list(A = A, B = B, improvement = imp,
    NAP = es(NAP), Tau = es(Tau), TauU = es(Tau_U), TauBC = es(Tau_BC),
    PND = es(PND), PEM = es(PEM), PAND = es(PAND), IRD = es(IRD),
    SMDb = es(SMD, std_dev = "baseline"), SMDp = es(SMD, std_dev = "pool"))
  if (all(c(A, B) >= 0)) {
    r$LRRi <- es(LRRi, scale = "count", pct_change = TRUE)
    r$LRRd <- es(LRRd, scale = "count", pct_change = TRUE)
  }
  pairs[[k]] <- r
}

## 2) BC-SMD: scdhlm veri setleri ve sentetik çoklu başlama düzeyi ----------
std <- function(d, case, session, outcome, trtvar, baseline) {
  data.frame(case = as.character(d[[case]]), session = as.numeric(d[[session]]),
             outcome = as.numeric(d[[outcome]]), trt = as.integer(as.character(d[[trtvar]]) != baseline))
}
sets <- list(
  Laski = std(Laski, "case", "time", "outcome", "treatment", "baseline"),
  Saddler = std(subset(Saddler, measure == "writing quality"), "case", "time", "outcome", "treatment", "baseline"),
  Schutte = std(Schutte, "case", "week", "fatigue", "treatment", "baseline"),
  BartonArwood = std(BartonArwood, "case", "session", "outcome", "condition", "A"),
  Rodriguez = std(Rodriguez, "case", "session", "outcome", "condition", "A"),
  AlberMorgan = std(AlberMorgan, "case", "session", "outcome", "condition", "baseline"),
  Datchuk = std(Datchuk, "case", "session", "outcome", "condition", "baseline"),
  CaseHarrisGraham = std(CaseHarrisGraham, "case", "session", "outcome", "condition", "baseline")
)
for (s in 1:12) {
  J <- sample(3:6, 1); rows <- list()
  for (j in 1:J) {
    nT <- sample(14:28, 1); start <- sample(4:(nT - 4), 1)
    u <- rnorm(1, 0, runif(1, 0.1, 1.5)); phi <- runif(1, -0.3, 0.6)
    e <- as.numeric(arima.sim(list(ar = phi), nT)) * runif(1, 0.5, 2)
    trt <- as.integer(1:nT > start)
    rows[[j]] <- data.frame(case = paste0("V", j), session = 1:nT,
      outcome = round(10 + u + runif(1, 0.5, 3) * trt + e, 2), trt = trt)
  }
  sets[[paste0("sim", s)]] <- do.call(rbind, rows)
}
# bir veri setinde oturum aralıkları (boşluklu oturumlar)
g <- sets$sim1; g <- g[!(g$session %in% c(5, 9, 13)), ]; sets$sim_gaps <- g

bcsmd <- lapply(names(sets), function(nm) {
  d <- sets[[nm]]
  d <- d[order(d$case, d$session), ]
  bad <- tapply(d$trt, d$case, function(x) length(unique(x)) < 2)
  if (any(bad)) d <- d[!d$case %in% names(bad)[bad], ]
  res <- tryCatch({
    fit <- calc_BCSMD(design = "MBP", case = case, phase = trt, session = session, outcome = outcome,
                      treatment_name = "1", FE_base = 0, RE_base = 0, FE_trt = 0, data = d,
                      summary = FALSE)
    m <- fit$model; vc <- as.numeric(nlme::VarCorr(m)[, "Variance"])
    list(g = fit$g_AB, SE = fit$SE_g_AB, nu = fit$nu, CI_L = fit$CI_L, CI_U = fit$CI_U, phi = fit$phi, rho = fit$rho,
         tau2 = vc[1], sigma2 = vc[2], beta = unname(nlme::fixef(m)), logLik = as.numeric(logLik(m)))
  }, error = function(e) list(error = conditionMessage(e)))
  list(name = nm, data = d, ref = res)
})

write_json(list(pairs = pairs, bcsmd = bcsmd), "ref.json",
           auto_unbox = TRUE, digits = NA, null = "null")
cat("ok\n")
