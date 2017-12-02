# outference

This R package provides a tool for doing valid statistical inference corrected for outlier removal in linear regressions. See [arXiv:1711.10635](https://arxiv.org/abs/1711.10635) [stat.ME] for more details. 

This repo contains the development version of the package. To install it with the vignette, type 
```{r}
devtools::install_github("shuxiaoc/outference", build_vignettes = TRUE)
```
in R. If you do not want to include the vignette, simply set `build_vignettes = FALSE`.
