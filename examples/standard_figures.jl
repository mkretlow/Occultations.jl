# The standard figures written by every example (included by the example scripts)
#
#   <name>_limb_fit.png           the limb fit, default display: 1σ shape band and centre ellipse
#   <name>_limb_fit_region.png    the same with the 68.3 % (1σ) region of whole limbs
#   <name>_limb_fit_paper.png     the limb in the style of the paper's figure (band type, level)
#   <name>_posterior.png          the limbs of the posterior (if a posterior is given)
#   <name>_geomap_prediction.png  the predicted shadow path (ephemeris only) with its uncertainty
#   <name>_geomap_postfit.png     the shadow path reconstructed from the fit
#   <name>_geomap_comparison.png  both paths together
#   <name>_globe.png              the post-fit path on the globe: Earth limb, night side, terminator

using Occultations, Printf, Plots
import CairoMakie

"""
    standard_figures(name, fit; out, posterior=nothing, event=1, title=<event name>,
                     paper=(fit=fit, uncertainty=:pointwise, nsigma=1, label="pointwise 1σ region"),
                     error_mas=20.0, error_note="assumed", projection, limits, window=600, marks=60)

Write the standard figures of an example (see the top of this file).
"""
function standard_figures(name::AbstractString, fit; out::AbstractString, posterior=nothing, event::Integer=1,
                          title=nothing, paper=(;), error_mas::Real=20.0, error_note::AbstractString="assumed",
                          projection::AbstractString, limits, window::Real=600, marks::Real=60)
    f(x) = joinpath(out, "$(name)_$(x).png")
    title = something(title, fit.problem.setup.events[event].name)
    pf = get(paper, :fit, fit)
    pu = get(paper, :uncertainty, :pointwise)
    pn = get(paper, :nsigma, 1)
    pl = get(paper, :label, pu === :region ? "68.3 % region" : "pointwise $(pn)σ region")
    # limb fits
    savefig(plot(fit; event=event, title=title, size=(1000, 750)), f("limb_fit"))
    savefig(plot(fit; event=event, uncertainty=:region, title="$title (68.3 % region of whole limbs)", size=(1000, 750)),
            f("limb_fit_region"))
    savefig(plot(pf; event=event, uncertainty=pu, nsigma=pn, title="$title (as the paper: $pl)", size=(1000, 750)),
            f("limb_fit_paper"))
    posterior === nothing ||
        savefig(plot(posterior; event=event, title="$title: posterior limbs", size=(1000, 750)), f("posterior"))
    # maps: prediction (ephemeris only), post-fit, both; and the globe
    ev = fit.problem.setup.events[event]
    pred = occ_map(ev; window=window, marks=marks, error_mas=error_mas)
    post = occ_map(fit; event=event, window=window, marks=marks)
    kw = (projection=projection, limits=limits)
    CairoMakie.save(f("geomap_prediction"),
                    geomap(pred; kw..., title=@sprintf("%s: prediction (ephemeris, %s σ = %g mas)", title, error_note, error_mas)))
    CairoMakie.save(f("geomap_postfit"), geomap(post; kw..., title="$title: post-fit shadow path"))
    CairoMakie.save(f("geomap_comparison"), geomap(pred, post; kw..., title="$title: prediction and post-fit path"))
    CairoMakie.save(f("globe"), geomap(post; size=(1000, 1000), title="$title: post-fit path, night side shaded"))
    return nothing
end
