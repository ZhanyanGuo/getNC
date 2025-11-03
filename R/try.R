suppressPackageStartupMessages({
  library(getNC)
})

say <- function(...) cat(sprintf(...), "\n")

say("Step 0: Fit glasso using the package's default dataset (obj = NULL)...")
fit <- fit_glasso(
  obj       = NULL,   # <-- uses your packaged pbmc_small (or whatever you bundled)
  nfeatures = 100,    # small for a quick run
  rho       = 0.2
)
say("  - features: %d", length(fit$features))
say("  - omega dim: %d x %d", nrow(fit$omega), ncol(fit$omega))
say("  - sigma dim: %d x %d", nrow(fit$sigma), ncol(fit$sigma))

say("Step 2: Conditional prediction (pick a few genes from the fit)...")
genes   <- fit$features
target  <- genes[5]
knocked <- genes[c(3, 7, 9)]
say("  - target:  %s", target)
say("  - knocked: %s", paste(knocked, collapse = ", "))

res <- predict_knockout_from_fit(
  fit,
  target  = target,
  knocked = knocked,
)
say("  - E[target | knocked=0] = %.4f", res$mean_cond[1])
say("  - Var[target | knocked=0] = %.4f", as.numeric(res$cov_cond[1,1]))

# 3D density plot for partner knockouts (uses sweep_partner_knockouts)
suppressPackageStartupMessages({
  library(plotly)   # for 3D plot
})

say("Step 3: Plotting partner knockout densities for target under {knocked ∪ g}...")
plots <- plot_partner_knockout_densities_dual(
  fit,
  target  = target,
  knocked = knocked,
  k = 20,                     # how many partners per plot
  use_original_units = TRUE   # or FALSE to work on correlation scale
)

# View them (in RStudio/VSCode/interactive R, just print):
plots$mean_plot
plots$var_plot


say("Done.")