# enviromtx 1.2.0

This is a minor release that provides a vignette and an example dataset. 

## Minor changes

* `example_data` is added and documented, "intro_enviromtx.Rmd" vignette is added. 

# enviromtx 1.1.0

This is a minor release that provides the option of centering covariates.

## Minor changes

* `fit_mgx_model` now includes a `control` argument. One possible inclusion in the `control` argument is `center`. This is by default `FALSE`. When `TRUE` all covariates will be centered before fitting the model. This should have negligible effects on estimation, but may improve stability in rare cases. 

# enviromtx 1.0.0

This is a major release that requires a new release of `raoBust` which performs many of the backend computation, changes the way data are input to `fit_mgx_model`. 

## Breaking changes

* `raoBust` version 1.1.1 is required, which updates the computation of robust score tests that are used in `fit_mgx_model()`
* the inputs to `fit_mgx_model()` are updated so that all vectors need to be included within the `enviro_df` data frame, and cannot be taken directly from the global environment.

