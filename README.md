Estimating net transition probabilities from cross-sectional data
================

## Reference

van de Kassteele J, Hoogenveen RT, Engelfriet PM, van Baal PHM,
Boshuizen HC (2012). Estimating net transition probabilities from
cross-sectional data with application to risk factors in chronic disease
modeling. *Statistics in Medicine* **31**(6), 533-543.
<http://onlinelibrary.wiley.com/doi/10.1002/sim.4423/abstract>

## Abstract

A problem occurring in chronic disease modeling is the estimation of
transition probabilities of moving from one state of a categorical risk
factor to another. Transitions could be obtained from a cohort study,
but often such data may not be available. However, under the assumption
that transitions remain stable over time, age specific cross-sectional
prevalence data could be used instead. Problems that then arise are
parameter identifiability and the fact that age dependent
cross-sectional data are often noisy or are given in age intervals. In
this paper we propose a method to estimate so-called net annual
transition probabilities from cross-sectional data, including their
uncertainties. Net transitions only describe the net inflow or outflow
into a certain risk factor state at a certain age. Our approach consists
of two steps: first, smooth the data using multinomial P-splines,
second, from these data estimate net transition probabilities. This
second step can be formulated as a transportation problem, which is
solved using the simplex algorithm from linear programming theory. A
sensible specification of the cost matrix is crucial to get meaningful
results. Uncertainties are assessed by parametric bootstrapping. We
illustrate our method using data on body mass index. We conclude that
this method provides a flexible way of estimating net transitions and
that the use of net transitions has implications for model dynamics, for
example when modeling interventions.

## Worked example

### Initial stuff

First we need to load some packages. If you don’t have them, install
them using `install.packages`.

``` r
# Load packages
library(tidyverse)
library(mgcv)
library(Matrix)
library(lpSolve)
```

Next, we load the data. Here we use the `bmi data.csv` on the Dutch
males ages 0-85 in 2006 and 2007 (POLS survey). The data are available
on this GitHub repository.

``` r
# Read BMI data of Dutch males aged 0-85 in 2006/2007
bmi.data <- read.csv(file = "bmi data.csv")

# Show first records of bmi.data
head(bmi.data)
```

``` 
  age normal overweight obese
1   0     32          1     1
2   1     47          0     2
3   2     52          0     2
4   3     55          0     1
5   4     46          0     0
6   5     68          0     0
```

### Data preparation

We see that `bmi.data` has four columns: `age` and the three BMI
categories: `normal`, `overweight`, `obese`. The numbers are counts. The
data are organized in wide format. For our calculations we need them in
long format (i.e. tidy format). More specifically, we need to have the
data at the individual level, where each person gets one record with
`age` and its BMI category. So 32 records with `age = 0` and `BMI_cat =
normal`, 1 record with `age = 0` and `BMI_cat = overweight`, etc.

We use the `gather` function to reshape the data from wide to long
format. However, the counts are still aggregated then. Therefore we
apply the `uncount` function to expand the records to the individual
level.

``` r
bmi.data <- bmi.data %>%
  # Reshape bmi.data into long format
  gather(key = "BMI_cat", value = "count", normal, overweight, obese) %>% 
  # Duplicate the rows according to the counts
  uncount(weights = count)

# Check
bmi.data %>% count(age, BMI_cat)
```

    # A tibble: 239 x 3
         age BMI_cat        n
       <int> <chr>      <int>
     1     0 normal        32
     2     0 obese          1
     3     0 overweight     1
     4     1 normal        47
     5     1 obese          2
     6     2 normal        52
     7     2 obese          2
     8     3 normal        55
     9     3 obese          1
    10     4 normal        46
    # … with 229 more rows

We notice that `obese` is before `overweight`, according to the
alphabet. We turn `BMI_cat` into a categorical variable (factor) having
the levels in the correct order.

``` r
# Make BMI_cat a factor with levels normal, overweight and obese
bmi.data <- bmi.data %>%
  mutate(
    BMI_cat = BMI_cat %>% factor(levels = c("normal", "overweight", "obese")))

# Check
bmi.data %>% count(age, BMI_cat)
```

    # A tibble: 239 x 3
         age BMI_cat        n
       <int> <fct>      <int>
     1     0 normal        32
     2     0 overweight     1
     3     0 obese          1
     4     1 normal        47
     5     1 obese          2
     6     2 normal        52
     7     2 obese          2
     8     3 normal        55
     9     3 obese          1
    10     4 normal        46
    # … with 229 more rows

bmi.data %\>% ggplot( mapping = aes( x = age, y = ..prop.., colour =
BMI\_cat)) + geom\_point()

### Fit multinomial P-splines

In order to estimate the net transition probabilties, we need to have a
smooth estimate of the age-specific prevalences for each BMI category.
We use the `gam` function from the `mgcv` package with a multinomial
likelihood. We need to specify the number of categories, counting from
0. As we have three categories, we need have `family = multinom(K = 2)`.

Furthermore, we need to translate the factor levels of `BMI_cat` into 0,
1 or 2.

``` r
# Create new variable BMI_int, 0 = normal, 1 = overweight, 2 = obese
bmi.data <- bmi.data %>%
  mutate(BMI_int = BMI_cat %>% as.integer %>% "-"(1))
```

Each category gets its own P-spline (`bs = "ps"`). Because the third
category is one minus the sum of the others, we only need to specify two
smooths of `age`. We must choose `k` large enough to allow enough
flexibility, the penalty does the rest. Because of the `method = REML`
smoothing parameter estimation, the result is somewhat different from
the fit in the paper. Also see `help(gam)`, `help(smooth.terms)` and
`help(multinom)` in the `mgcv` package.

``` r
mod <- gam(
  formula = list(BMI_int ~ s(age, bs = "ps", k = 15), ~ s(age, bs = "ps", k = 15)),
  family = multinom(K = 2),
  data = bmi.data %>% mutate(BMI_int = BMI_cat %>% as.integer %>% "-"(1)),
  method = "REML")
```

### Predict age-specific prevalences
