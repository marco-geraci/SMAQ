# =========================================================
# Statistical Methods and Applications Using Quantiles
# Master package loader
# Marco Geraci
# =========================================================

pkgs <- c(
  "AER", "EnvStats", "KernSmooth", "MASS", "Qtools",
  "bbmle", "boot", "cdfquantreg", "conquer", "datasets",
  "discreteQ", "evd", "fitdistrplus", "ggplot2", "gridExtra",
  "lpSolve", "lqmm", "microbenchmark", "mplot",
  "nlme", "nlqmm", "numDeriv", "pak",
  "quantreg", "reldist", "reshape2",
  "stats", "unitquantreg", "xtable"
)

invisible(lapply(pkgs, library, character.only = TRUE))


###################################################
# Set ggplot global options

# lightgrey (#D3D3D3), gainsboro (#DCDCDC), whitesmoke (#F5F5F5)
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(margin = margin(b = 10), vjust = 1),

      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.x.top = element_text(margin = margin(b = 10)),
      axis.text.x.top = element_text(margin = margin(b = 4)),

      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.y.right = element_text(margin = margin(l = 10)),
      axis.text.y.right = element_text(margin = margin(l = 4)),

      panel.grid.major = element_line(color = "gainsboro"),
      panel.grid.minor = element_line(color = "gainsboro")
    )
)

# pretty rounding format
fround <- function(x, digits = 2){
    format(round(x = x, digits = digits), nsmall = digits)
}
