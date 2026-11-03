# ---------------------------------------------------------------------------
# A bilobed body: two triaxial ellipsoids on the body x-axis (a contact binary)
# ---------------------------------------------------------------------------

"""
    BilobedEllipsoid(lobe1::Ellipsoid, lobe2::Ellipsoid, d)

Contact binary of two triaxial ellipsoids with parallel axes, centred on the body x-axis at
the distance `d` [km] from each other (lobe 2 towards +x). The origin of the body frame is
the centre of mass for equal densities, i.e. lobe 1 lies at `x₁ = -d V₂/(V₁+V₂)` and lobe 2
at `x₂ = d V₁/(V₁+V₂)` (`Vᵢ` the volumes of the ellipsoids). With `d < a₁ + a₂` the lobes
overlap (interpenetrate). The projection ([`project`](@ref)) is exact: the union of the two
projected ellipses, a [`CompositeOutline`](@ref).

```jldoctest
julia> b = BilobedEllipsoid(Ellipsoid(20.0, 18.0, 10.0), Ellipsoid(15.0, 14.0, 10.0), 33.0);

julia> round(volume_equivalent_radius(b); digits=2)       # with the overlap counted once
17.85
```
"""
struct BilobedEllipsoid{T<:Real} <: AbstractShape
    lobe1::Ellipsoid{T}
    lobe2::Ellipsoid{T}
    d::T
end
function BilobedEllipsoid(e1::Ellipsoid, e2::Ellipsoid, d::Real)
    T = promote_type(typeof(e1.a), typeof(e2.a), typeof(d))
    return BilobedEllipsoid{T}(Ellipsoid(T(e1.a), T(e1.b), T(e1.c)), Ellipsoid(T(e2.a), T(e2.b), T(e2.c)), T(d))
end

Base.show(io::IO, b::BilobedEllipsoid) = print(io, "BilobedEllipsoid(", b.lobe1, ", ", b.lobe2,
                                               @sprintf(", d = %.3f km)", b.d))

# x positions of the lobe centres (centre of mass at the origin)
function _lobe_x(b::BilobedEllipsoid)
    V1, V2 = volume(b.lobe1), volume(b.lobe2)
    return -b.d * V2 / (V1 + V2), b.d * V1 / (V1 + V2)
end

function project(b::BilobedEllipsoid, M::AbstractMatrix, u::AbstractVector; n::Integer=720)
    S = _sky_matrix(u)
    x1, x2 = _lobe_x(b)
    p = S * (M' * [1, 0, 0])                    # the body x-axis in sky coordinates
    return CompositeOutline([project(b.lobe1, M, u), project(b.lobe2, M, u)],
                            [(x1 * p[1], x1 * p[2]), (x2 * p[1], x2 * p[2])])
end

# volume of the intersection of the two lobes: slices perpendicular to x are two concentric,
# aligned ellipses; their common area is ½∫ min(r₁², r₂²) dφ (Gauss–Legendre in x and φ)
function _overlap_volume(b::BilobedEllipsoid; n::Integer=48)
    x1, x2 = _lobe_x(b)
    e1, e2 = b.lobe1, b.lobe2
    lo, hi = max(x1 - e1.a, x2 - e2.a), min(x1 + e1.a, x2 + e2.a)
    lo < hi || return zero(b.d)
    xs, ws = _gauss_legendre(n)
    φs, wφ = _gauss_legendre(n)
    V = zero(b.d)
    for (ξ, w) in zip(xs, ws)
        x = (lo + hi) / 2 + (hi - lo) / 2 * ξ
        s1 = max(1 - ((x - x1) / e1.a)^2, 0.0)
        s2 = max(1 - ((x - x2) / e2.a)^2, 0.0)
        A = zero(V)
        for (η, v) in zip(φs, wφ)
            φ = π / 4 * (η + 1)                 # first quadrant, ×4 by symmetry
            c2, sn2 = cos(φ)^2, sin(φ)^2
            r1 = s1 / (c2 / e1.b^2 + sn2 / e1.c^2)
            r2 = s2 / (c2 / e2.b^2 + sn2 / e2.c^2)
            A += v * min(r1, r2)
        end
        V += w * 4 * (π / 4) * A / 2
    end
    return V * (hi - lo) / 2
end

volume(b::BilobedEllipsoid) = volume(b.lobe1) + volume(b.lobe2) - _overlap_volume(b)

"""
    ellipsoid_mesh(b::BilobedEllipsoid; subdivisions=4, name) -> TriMesh

Triangle mesh of a [`BilobedEllipsoid`](@ref): the meshes of the two lobes, placed on the
x-axis (they interpenetrate where the lobes overlap; fine for viewing and for OBJ export).
"""
function ellipsoid_mesh(b::BilobedEllipsoid; subdivisions::Integer=4,
                        name::AbstractString="bilobed ellipsoid")
    x1, x2 = _lobe_x(b)
    m1 = ellipsoid_mesh(b.lobe1; subdivisions=subdivisions)
    m2 = ellipsoid_mesh(b.lobe2; subdivisions=subdivisions)
    V = hcat(m1.vertices .+ [x1, 0, 0], m2.vertices .+ [x2, 0, 0])
    F = hcat(m1.faces, m2.faces .+ size(m1.vertices, 2))
    return TriMesh(V, F; name=name)
end

# Gauss–Legendre nodes and weights on [-1, 1] (Golub–Welsch)
const _GL_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
function _gauss_legendre(n::Integer)
    get!(_GL_CACHE, n) do
        β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
        E = eigen(SymTridiagonal(zeros(n), β))
        (E.values, 2 .* E.vectors[1, :] .^ 2)
    end
end
