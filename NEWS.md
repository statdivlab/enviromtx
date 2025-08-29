# enviromtx 1.0.0

This is a major release that requires a new release of `raoBust` which performs many of the backend computation, changes the way data are input to `fit_mgx_model`. It also provides the option of centering covariates.

## Breaking changes

* `raoBust` version 1.1.1 is required, which updates the computation of robust score tests that are used in `fit_mgx_model()`
* the inputs to `fit_mgx_model()` are updated so that all vectors need to be included within the `enviro_df` data frame, and cannot be taken directly from the global environment.
* a `control` argument is added for `fit_mgx_model()` with the option `center = TRUE` by default, which will center all covariates in the model.

## Minor changes

* `fit_mgx_model` now includes a `control` argument. One possible inclusion in the `control` argument is `center`. This is by default `FALSE`, but when it is `TRUE` all covariates will be centered before fitting the model. This should have very small if any effects on estimation, but could affect stability. 
