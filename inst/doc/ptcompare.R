## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.6, dpi = 120,
                      out.width = "100%", fig.align = "center")
has_mclust <- requireNamespace("mclust", quietly = TRUE)

## -----------------------------------------------------------------------------
library(ptcompare)

## -----------------------------------------------------------------------------
set.seed(11)
n <- 1000L
x <- matrix(rnorm(2 * n), ncol = 2,
            dimnames = list(NULL, c("marker_A", "marker_B")))
y <- matrix(rnorm(2 * n), ncol = 2,
            dimnames = list(NULL, c("marker_A", "marker_B")))

extra <- runif(n) < 0.15
y[extra, ] <- cbind(rnorm(sum(extra), 1.6, 0.3),
                    rnorm(sum(extra), 1.6, 0.3))
sum(extra)

## ----fig.cap = "The two samples. Group 2 carries an extra cluster near (1.6, 1.6)."----
plot(rbind(x, y), type = "n", xlab = "marker_A", ylab = "marker_B")
points(x, pch = 0, col = "#123E7C", cex = 0.6)
points(y, pch = 2, col = "#B65D35", cex = 0.6)
legend("topleft", c("group 1", "group 2"), pch = c(0, 2),
       col = c("#123E7C", "#B65D35"), bty = "n")

## -----------------------------------------------------------------------------
set.seed(12)
fit <- compare_pt(x, y)
fit

## -----------------------------------------------------------------------------
regions(fit, prob = 0.9, min_n = 40)

## ----fig.height = 5.4, fig.cap = "The four most probable regions. The box is the region, the dashed line its threshold."----
plot(fit, type = "regions", top = 4, min_n = 40)

## ----fig.cap = "The representative tree, drawn to depth 3."-------------------
plot(fit, max_depth = 3)

## ----fig.height = 5.4, fig.cap = "Regions by depth, shaded by the probability of a difference."----
plot(fit, type = "levels", max_depth = 3)

## ----eval = has_mclust--------------------------------------------------------
data("wdbc", package = "mclust")
features <- setdiff(names(wdbc)[sapply(wdbc, is.numeric)], "ID")
features
table(wdbc$Diagnosis)

## ----eval = has_mclust--------------------------------------------------------
set.seed(2)
fit_wdbc <- compare_pt(wdbc[, features], group = wdbc$Diagnosis)
summary(fit_wdbc, top = 3)

## ----eval = has_mclust--------------------------------------------------------
found <- regions(fit_wdbc, prob = 0.9, min_n = 25)
found[, c("depth", "variable", "n1", "n2", "prob_difference", "effect")]

## ----eval = has_mclust, fig.cap = "Left: the division at the root. Right: a division within the large nuclei."----
plot(fit_wdbc, type = "regions", top = 2, min_n = 25)

## ----eval = has_mclust--------------------------------------------------------
local <- found[found$depth > 0, ][1, ]
local$condition

## ----eval = has_mclust--------------------------------------------------------
cells <- region_data(fit_wdbc, local$node_id)
tapply(cells$Nconcave_mean, cells$group, median)
tapply(cells$Nconcave_mean <= local$threshold, cells$group, mean)

## ----eval = has_mclust--------------------------------------------------------
benign <- wdbc[wdbc$Diagnosis == "B", features]
set.seed(3)
half <- sample(c(TRUE, FALSE), nrow(benign), replace = TRUE)

set.seed(4)
fit_null <- compare_pt(benign[half, ], benign[!half, ])
fit_null
nrow(regions(fit_null, prob = 0.9))

