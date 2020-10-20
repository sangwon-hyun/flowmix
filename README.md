flowmix
=============

Sparse mixture of regressions for flow cytometry data.

This package implements a sparse mixture of multivariate regressions model for
flow cytometry data. The main motivating application is continuous-time flow
cytometry data collected in the ocean, over space and time.

The paper preprint is here: [https://arxiv.org/abs/2008.11251](https://arxiv.org/abs/2008.11251).

## Installation

This R package can be installed using the following commands.

```{r}
library(devtools)
install_github("robohyun66/flowcy", subdir = "flowcy")
```

## Demo

For a short demo about using the package, see [demo.Rmd](demo.Rmd). A
longer vignette is in the works.
	
## Authors & Contributors

The contributors are (Sangwon Hyun)[http://sangwon-hyun.org/], Mattias Rolf Cape, Francois Ribalet, and Jacob Bien.
