# API reference

All exported functions and types, grouped by topic. The package module:

```@docs
Occultations
```

```@contents
Pages = ["api.md"]
Depth = 2:2
```

## Time and angles

```@autodocs
Modules = [Occultations]
Pages = ["core/time.jl", "core/angles.jl", "core/constants.jl"]
Private = false
```

## Stars

```@autodocs
Modules = [Occultations]
Pages = ["astro/star.jl"]
Private = false
```

## Ephemerides and bodies

```@autodocs
Modules = [Occultations]
Pages = ["astro/ephemeris.jl"]
Private = false
```

## Observers

```@autodocs
Modules = [Occultations]
Pages = ["astro/observer.jl"]
Private = false
```

```@docs
AbstractObserver
```

## Geometry, shadow paths and predictions

```@autodocs
Modules = [Occultations]
Pages = ["astro/geometry.jl", "astro/prediction.jl"]
Private = false
```

## Orientation and shapes

```@autodocs
Modules = [Occultations]
Pages = ["body/rotation.jl", "body/shape.jl", "body/bilobed.jl"]
Private = false
```

```@docs
AbstractShape
```

## Outlines

```@autodocs
Modules = [Occultations]
Pages = ["body/outline.jl", "body/composite.jl"]
Private = false
```

## Light curves and photometry

```@autodocs
Modules = [Occultations]
Pages = ["photometry/lightcurve.jl", "photometry/diffraction.jl", "photometry/timing.jl"]
Private = false
```

## Atmospheres

```@autodocs
Modules = [Occultations]
Pages = ["photometry/atmosphere.jl"]
Private = false
```

## Events, chords and tracks

```@autodocs
Modules = [Occultations]
Pages = ["event/event.jl", "event/chord.jl", "event/track.jl", "event/cora_report.jl"]
Private = false
```

## Limb models and limb fits

```@autodocs
Modules = [Occultations]
Pages = ["limb/models.jl", "limb/binary.jl", "limb/rings.jl", "limb/problem.jl", "limb/results.jl"]
Private = false
```

## Inference: parameters, fits and posteriors

```@autodocs
Modules = [Occultations]
Pages = ["fit/params.jl", "fit/solve.jl", "fit/bayes.jl"]
Private = false
```

## Plots and maps

Plots.jl recipes (`plot(fit)`, `plot(posterior)`, `plot(lightcurve)`, …), maps with GeoMakie
and 3-D views with Makie.

```@autodocs
Modules = [Occultations]
Pages = ["plots/recipes.jl", "plots/maps.jl"]
Private = false
```

```@docs
geomap
meshview
occview
occanimate
```
