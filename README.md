# ptcompare: interpretable multivariate two-sample comparison

`ptcompare` compares the distributions of two samples of continuous
multivariate observations. It returns the posterior probability that the two
distributions differ, and describes where they differ as regions defined by
thresholds on the variables. The method is the hidden Markov Pólya tree of
Awaya and Ma (2024). Adaptive partition trees are explored by sequential Monte
Carlo (SMC), and the latent states given a tree are computed by message
passing.

## Installation

From CRAN:

```r
install.packages("ptcompare")
```

The development version from GitHub needs a C++ compiler (Rtools on Windows):

```r
# install.packages("remotes")
remotes::install_github("nawaya040/ptcompare")
```

## Usage

```r
library(ptcompare)

set.seed(1)
x <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("marker_A", "marker_B")))
y <- matrix(rnorm(400, mean = 0.8), ncol = 2, dimnames = list(NULL, c("marker_A", "marker_B")))

fit <- compare_pt(x, y)
fit                            # posterior probability of an overall difference
summary(fit)                   # settings, computation report, divisions with the largest differences
regions(fit, prob = 0.9)       # table of regions where the groups differ

plot(fit)                      # representative tree
plot(fit, type = "regions")    # regions with the largest effects, drawn over the data
plot(fit, type = "levels")     # two variables only: regions at each depth shaded by probability of a difference
```

## Vignette

```r
vignette("ptcompare")
```

The vignette has three parts: a simulation in two variables where a difference
placed in a small part of the plane is recovered, an analysis of 30
measurements of cell nuclei from `mclust::wdbc` where the groups are two
diagnoses and the difference is not explained by nucleus size alone, and a
pair of samples drawn from the same distribution.

## Comparing your own data

`x` and `y` are numeric matrices, or data frames with only numeric columns,
with the same columns in the same order.

```r
library(ptcompare)
set.seed(123)
fit <- compare_pt(x, y)

# Data with both groups in one data frame and a group column
fit <- compare_pt(data, group = condition)
```

The two samples are treated as independent. `NA` and `Inf` values give an
error. Columns that are constant over both groups combined are dropped with a
warning. When the group variable is a factor, its levels give group 1 and
group 2 in order; otherwise the order of first appearance is used.

The computational budget and the statistical model are set by separate
arguments.

```r
fit <- compare_pt(
  x, y,
  control = pt_control(particles = 512, effect_draws = 4096),
  prior = pt_prior(max_depth = 8, split_grid = 15, min_node = 20)
)
diagnostics(fit)
```

`control` sets how the posterior is computed. The number of particles
determines the accuracy of the SMC approximation. `prior` sets the model: the
maximum tree depth, the number of candidate split points and the prior
probabilities define the posterior being computed.

The default `pt_prior(midpoint_spike = TRUE)` is the spike-and-slab prior on
split locations in Appendix A of the paper: once a region is divided at its
midpoint, all of its descendants are also divided at their midpoints. With
`midpoint_spike = FALSE`, equation (2) of the paper is used at every node.

Resampling follows Section 3.1 of the paper. Ancestors are drawn with
probability proportional to the weights raised to the power κ, and the weights
are then corrected by W/a. The default is `pt_control(resample_power = 0.5)`;
`resample_power = 1` gives standard resampling.

## Transformation of the variables

The model is defined on the unit cube, so each variable is first mapped to
[0, 1]. Without a `support` argument, `compare_pt(x, y)` uses a rank
transformation computed from both groups together. When the known range of
each variable is given as `support = cbind(lower, upper)`, a linear
transformation is the default. The transformation used is stored in the result
and shown by `summary()`.

```r
fit_rank  <- compare_pt(x, y)                       # default: transform = "rank"
fit_paper <- compare_pt(x, y, transform = "range")  # linear transformation as in the paper
```

`transform = "range"` maps each variable linearly onto the unit interval, as in
the paper. When `support` is omitted, the common range of the two groups is
widened by 1% on each side. `transform = "rank"` maps the midranks of the
pooled sample onto [0, 1]; the groups are not ranked separately. Ties are kept
as ties and no jitter is added. Variable names, thresholds and region
conditions are reported in the original units through the stored inverse
transformation.

## Reading the results

`fit$prob_difference` is the posterior probability, under the model, that the
two distributions differ anywhere in the data space. It is averaged over the
trees in the SMC sample. A value near 1 means the data support a difference; a
value near 0 means they support equal distributions.

Each row of `regions()` is an internal node of the representative tree, which
is the sampled tree with the highest posterior score. Inside the region given
by `condition`, the observations are divided into those with
`variable <= threshold` (the left part) and the rest. The row describes whether
the two groups fall into the left part in different proportions:

- `prob_difference`: posterior probability that the proportions differ at this
  division.
- `effect`: posterior mean of the log odds ratio of falling into the left part,
  group 1 versus group 2. Positive values mean group 1 is more likely to be in
  the left part.
- `effect_abs`: posterior mean of the absolute log odds ratio, the size of the
  difference regardless of direction.
- `lower`, `upper`: 95% credible interval for the log odds ratio. The posterior
  includes a point mass at no difference, so a limit can be exactly 0.

By default `regions()` sorts the rows by `prob_difference`. Use
`rank_by = "effect"` to sort by effect size; `plot(fit, type = "regions")` uses
the same argument.

`region_rows()` and `region_data()` return the observations in one region.

```r
top <- regions(fit, prob = 0.9)$node_id[1]   # region with the largest difference
region_data(fit, top)                        # its observations, with group labels
mydata[region_rows(fit, top), ]              # the same rows of the original data frame
```

`region_rows()` returns one `TRUE`/`FALSE` per observation, in the order passed
to `compare_pt()`: the rows of `x` followed by the rows of `y`, or the row
order of the original data frame when `group =` is used. Both functions need
the data stored in the result, which is the default
(`pt_control(keep_data = TRUE)`).

Tree plots are cut at a display depth, and `+` marks a node with further
branches below it. The cut affects only the plot. `plot(fit, type = "regions")`
draws the regions with the largest differences (by default those with
`prob_difference` of at least 0.5 and at least 50 observations of each group,
up to 9 regions). With two or more variables, each region is drawn as a box on
a scatter plot, with the observations
inside the box marked by group. With one variable, the distributions of the two
groups within the region are drawn. `node =` selects specific regions; a single
region is drawn as one large panel. `type = "levels"` (two variables only)
shades the regions at each depth of the tree by their probability of a
difference.

The model, the SMC target distribution, the overall null probability and the
effect sizes are defined in the statistical specification.

```r
file.show(system.file("STATISTICAL-SPECIFICATION.md", package = "ptcompare"))
```

## Using the repository

The following files are in the source repository and are not part of the
installed package.

Open **ptcompare.Rproj** and run the following once in the Console. It
installs the package into `.local-library` inside the project. After changing
the source, run Session → Restart R in RStudio before installing again.

```r
source("scripts/install-local.R")
source("examples/01-local-shift.R")
```

The script leaves `fit`, `region_table`, `x` and `y` in the Environment and
draws the tree, the comparison of local regions and the computational
diagnostics in the Plots pane. Use the back arrow in the Plots pane to see
earlier plots.

| File | Contents |
| --- | --- |
| `examples/00-tryout.R` | Basic use, reading the results, an example with `iris`, and changing the settings (about 10 seconds in total) |
| `examples/00-tryout-dimensions.R` | The same steps in one, two and three variables, to compare how the output looks in each case |
| `examples/01-local-shift.R` | Two variables, with a difference only in a known subpopulation |
| `examples/02-dependence-and-null.R` | Equal marginal distributions with different correlation (two and three variables), and two samples from the same distribution |
| `examples/03-stability.R` | The same data and model fitted with different numbers of particles and random seeds |
| `examples/04-high-dimensional.R` | Twenty variables, fitted with the default rank transformation and with the linear transformation of the paper |

The third example fits the model several times and takes longer than a single
fit. Run time depends on the data and the model settings.

## Verification and reproduction

These files are also only in the repository.

- `tests/testthat/`: checks of fixed and adaptive trees against an independent
  enumeration in R, and tests of the API.
- `evaluation/REPORT.md`: record of the checks made during development
  (repeated simulations, comparison methods, computational cost).
- `evaluation/results/`: simulation results as CSV files with seeds, and the
  computing environment.

```r
source("scripts/load-local.R")
testthat::test_dir("tests/testthat")

# Rerun the evaluation. The study and rank-study scripts resume from existing CSV files.
source("evaluation/run-study.R")
source("evaluation/run-rank-study.R")
source("evaluation/run-verification.R")
source("evaluation/summarize.R")
```

The comparison with other methods also uses the MRS package. Memory use in a
separate process on Windows is measured by `evaluation/run-benchmarks.ps1`, and
the benchmark of the computational kernel is `evaluation/kernel-benchmark.cpp`.

## Notes

- With `transform = "rank"`, the posterior is conditional on the estimated
  transformation; its uncertainty is not integrated over. The transformation
  also changes the split prior and the terminal distributions in the original
  coordinates.
- The local quantities from `regions()` are conditional on the representative
  tree. Regions can overlap through parent–child relations, and selecting
  regions by a probability threshold carries no frequentist multiple-testing
  guarantee.

## References

Awaya, N. and Ma, L. (2024). Hidden Markov Pólya trees for high-dimensional
distributions. *Journal of the American Statistical Association*, 119(545),
189–201. [doi:10.1080/01621459.2022.2105223](https://doi.org/10.1080/01621459.2022.2105223)

The state model follows the MRS model of Soriano, J. and Ma, L. (2017).
Probabilistic multi-resolution scanning for two-sample differences. *Journal of
the Royal Statistical Society: Series B*, 79(2), 547–572.
[doi:10.1111/rssb.12180](https://doi.org/10.1111/rssb.12180)
