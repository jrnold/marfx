---
title: "README.Rmd"
author: "Jeffrey Arnold"
date: "May 1, 2016"
output: html_document
---

Yet another package to calculate 

# Methods

- `postsimev` : Simulate the expected value of the response, $E(y)$.
- `postsimy`: Simulate the response, $y$.
- `postsim`: Simulate parameters of the model.
- `postsim_partialfx`: Simulate `\delta y/ \delta x`.
- `partialfx`: Calculate partial effects (discrete marginal effects)
- `avg_partialfx`: Calculate the average partial effect
- `marfx`: Calculate the marginal effects
- `avg_margx`: Calculate the average marginal effect

# Supported classes

- `lm`
- `glm`

This is still a work in progress and even for these objects, not all variations
may be supported.

# Similar packages

[margins](https://github.com/leeper/margins) is the most similar, but implements
the derivatives of the model matrix with respect to the data variables via symbolic differentiation,
and calculates the variance-covariance of the marginal effects via the Delta rule.
**marfx** uses numerical derivatives to calculate the derivatives of the model
matrix with respect to the variables, which is slower but will handle all formulae.
**marfx** uses simulation rather than the delta method to calculate the variance
covariance of the marginal effects, which is slower, but more flexible.

**marfx** adopts some inspiration from [simcf](http://faculty.washington.edu/cadolph/?page=60)
which was written by Chris Adolph while teaching the same course which prompted me
to write this package.

[Zelig](http://zeligproject.org/ allows for partial (finite) differences via specified counterfactuals,
but only works for models within its package.

Several other packages provide estimates of the marginal effects for some types
of models, including [car](https://cran.r-project.org/web/packages/car/index.html),
[alr](https://cran.r-project.org/web/packages/alr3/index.html), [mfx](https://cran.r-project.org/web/packages/mfx/index.html),
and [erer](https://cran.r-project.org/web/packages/erer/index.html).
But these packages do not account for interrelations between variables.
