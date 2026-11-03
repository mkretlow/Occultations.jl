module OccultationsSPICEExt

using Occultations
using Occultations: EphemKernel, Epoch, et
using SPICE

function _load!(k::EphemKernel)
    k.loaded && return
    for f in k.kernels
        isfile(f) || error("SPICE kernel not found: $f")
        furnsh(f)
    end
    k.loaded = true
end

# Geocentric position, ICRF ("J2000" in SPICE), km; abcorr "CN" = astrometric (converged light time, as JPL Horizons)
function Occultations._spice_position(k::EphemKernel, t::Epoch)
    _load!(k)
    pos, _ = spkpos(k.target, et(t), "J2000", k.abcorr, "399")
    return collect(Float64, pos)
end

end # module
