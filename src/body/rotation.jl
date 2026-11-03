# ---------------------------------------------------------------------------
# Body-fixed (planetocentric) frames and body orientation
# ---------------------------------------------------------------------------

"""
    BodyFrame(; pole_ra, pole_dec, W0=0, rotation_rate=0, epoch=JD_J2000,
                pole_ra_rate=0, pole_dec_rate=0)

Body-fixed frame in the IAU convention (Archinal et al. 2018): the pole is at
`α = pole_ra + pole_ra_rate T`, `δ = pole_dec + pole_dec_rate T` and the prime
meridian at `W = W0 + rotation_rate d`, with `d` days and `T` Julian centuries
(TDB) since `epoch` (a TDB Julian date). Angles in degrees, `rotation_rate` in
deg/day (negative for retrograde rotators), pole rates in deg/century.

The frame is right-handed and longitudes are planetocentric, positive towards east.
"""
Base.@kwdef struct BodyFrame
    pole_ra::Float64
    pole_dec::Float64
    W0::Float64 = 0.0
    rotation_rate::Float64 = 0.0
    epoch::Float64 = JD_J2000
    pole_ra_rate::Float64 = 0.0
    pole_dec_rate::Float64 = 0.0
end

"""
    BodyFrame(pole_ra, pole_dec, rotation_period; W0=0, epoch=JD_J2000)

Frame from a fixed pole [deg] and a sidereal `rotation_period` [hours]
(negative for retrograde rotation).
"""
BodyFrame(pole_ra::Real, pole_dec::Real, rotation_period::Real; W0::Real=0.0, epoch::Real=JD_J2000) =
    BodyFrame(; pole_ra=pole_ra, pole_dec=pole_dec, W0=W0, epoch=epoch,
              rotation_rate=isfinite(rotation_period) ? 360 * 24 / rotation_period : 0.0)

Base.show(io::IO, f::BodyFrame) = @printf(io, "BodyFrame(pole = (%.4f°, %+.4f°), W = %.4f° + %.6f°/day)",
                                          f.pole_ra, f.pole_dec, f.W0, f.rotation_rate)

"""
    orientation_at(frame, t) -> (pole_ra, pole_dec, W)

Pole [deg] and prime-meridian angle [deg, in [0, 360)] at `t` (an [`Epoch`](@ref);
the TDB time is used).
"""
function orientation_at(fr::BodyFrame, t)
    t1, t2 = tdb_jd(_epoch(t))
    d = (t1 - fr.epoch) + t2
    T = d / 36525
    return fr.pole_ra + fr.pole_ra_rate * T, fr.pole_dec + fr.pole_dec_rate * T,
           mod(fr.W0 + fr.rotation_rate * d, 360.0)
end

# passive rotations
_R1(θ) = (c = cosd(θ); s = sind(θ); [1.0 0 0; 0 c s; 0 -s c])
_R3(θ) = (c = cosd(θ); s = sind(θ); [c s 0.0; -s c 0; 0 0 1])

"""
    body_matrix(frame, t) -> 3×3 matrix

Rotation from ICRS to the body-fixed frame at `t`: `x_body = M * x_icrs`.
"""
function body_matrix(fr::BodyFrame, t)
    α, δ, W = orientation_at(fr, t)
    return _R3(W) * _R1(90 - δ) * _R3(90 + α)
end

function _frame(b)
    b.frame === nothing && throw(ArgumentError("body $(b.name): no body frame (give `pole` and `rotation_period`)"))
    return b.frame
end

"""
    orientation(body, t; observer=GEOCENTER) -> NamedTuple

Orientation of `body` as seen by `observer` at `t` (body frame evaluated at the
light-time corrected epoch):

* `sub_observer_lon`, `sub_observer_lat` – planetocentric sub-observer point [deg]
* `sub_solar_lon`, `sub_solar_lat`       – sub-solar point [deg]
* `pole_position_angle` – position angle of the north pole, from north through east [deg]
* `pole_aperture_angle` – angle between the pole and the sky plane, positive when the
  north pole points towards the observer [deg]
* `matrix` – ICRS → body rotation, `light_time` [s]
"""
function orientation(b, t; observer::AbstractObserver=GEOCENTER)
    t = _epoch(t)
    fr = _frame(b)
    eph = _ephem(b)
    Bt = _topocentric_vector(eph, observer, t)
    D = norm(Bt)
    lt = D / C_KMS
    M = body_matrix(fr, t - lt)
    so = M * (-Bt ./ D)
    sun = -earth_heliocentric(t) .* AU_KM .- geocentric_vector(eph, t)
    ss = M * (sun ./ norm(sun))
    α, δ, _ = orientation_at(fr, t - lt)
    ra, dec, _ = vec2radec(Bt)
    return (sub_observer_lon=mod(atand(so[2], so[1]), 360), sub_observer_lat=asind(clamp(so[3], -1, 1)),
            sub_solar_lon=mod(atand(ss[2], ss[1]), 360), sub_solar_lat=asind(clamp(ss[3], -1, 1)),
            pole_position_angle=position_angle(ra, dec, α, δ), pole_aperture_angle=asind(clamp(so[3], -1, 1)),
            matrix=M, light_time=lt)
end

"""
    pole_position_angle(body, t; observer=GEOCENTER) -> (position_angle, aperture_angle) [deg]

Position angle (from north through east) and aperture angle of the pole of
`body` (its `frame`) as seen by `observer`.
"""
function pole_position_angle(b, t; observer::AbstractObserver=GEOCENTER)
    o = orientation(b, t; observer=observer)
    return o.pole_position_angle, o.pole_aperture_angle
end
