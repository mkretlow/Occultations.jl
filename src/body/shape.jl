# ---------------------------------------------------------------------------
# 3-D shapes (triaxial ellipsoid, triangle meshes) and their projection on the sky
# ---------------------------------------------------------------------------


"""
    Ellipsoid(a, b=a, c=b)

Triaxial ellipsoid with semi-axes `a`, `b`, `c` [km] along the body x, y, z axes.
"""
struct Ellipsoid{T<:Real} <: AbstractShape
    a::T
    b::T
    c::T
end
Ellipsoid(a::Real, b::Real=a, c::Real=b) = Ellipsoid(promote(a, b, c)...)
Base.show(io::IO, e::Ellipsoid) = @printf(io, "Ellipsoid(%.3f, %.3f, %.3f km)", e.a, e.b, e.c)

volume(e::Ellipsoid) = 4π / 3 * e.a * e.b * e.c

"""
    TriMesh(vertices, faces; name="")

Triangle mesh: `vertices` is `3 × nv` [km, body frame], `faces` `3 × nf` (1-based,
counter-clockwise seen from outside). See [`read_obj`](@ref), [`write_obj`](@ref),
[`ellipsoid_mesh`](@ref).
"""
struct TriMesh <: AbstractShape
    name::String
    vertices::Matrix{Float64}
    faces::Matrix{Int}
    edges::Vector{Tuple{Int,Int}}
end

function TriMesh(vertices::AbstractMatrix, faces::AbstractMatrix; name::AbstractString="")
    size(vertices, 1) == 3 && size(faces, 1) == 3 || throw(DimensionMismatch("vertices and faces must be 3 × n"))
    edges = Set{Tuple{Int,Int}}()
    for j in axes(faces, 2), (a, b) in ((1, 2), (2, 3), (3, 1))
        i1, i2 = faces[a, j], faces[b, j]
        push!(edges, i1 < i2 ? (i1, i2) : (i2, i1))
    end
    return TriMesh(String(name), Matrix{Float64}(vertices), Matrix{Int}(faces), sort!(collect(edges)))
end

Base.show(io::IO, m::TriMesh) = print(io, "TriMesh(", isempty(m.name) ? "" : m.name * ", ",
                                      size(m.vertices, 2), " vertices, ", size(m.faces, 2), " faces)")

"""
    read_obj(file; scale=1, name=basename(file)) -> TriMesh

Read a Wavefront OBJ file (`v` and `f` records; polygons are fan-triangulated; texture
and normal indices are ignored), multiplying the coordinates by `scale` (e.g. to km).
"""
function read_obj(file::AbstractString; scale::Real=1.0, name::AbstractString=basename(file))
    V = Float64[]; F = Int[]
    nv = 0
    for line in eachline(file)
        w = split(line)
        isempty(w) && continue
        if w[1] == "v"
            append!(V, parse.(Float64, w[2:4])); nv += 1
        elseif w[1] == "f"
            idx = [(k = parse(Int, first(split(x, '/'))); k < 0 ? nv + 1 + k : k) for x in w[2:end]]
            for j in 2:length(idx)-1
                append!(F, (idx[1], idx[j], idx[j+1]))
            end
        end
    end
    nv == 0 && throw(ArgumentError("no vertices found in $file"))
    return TriMesh(scale .* reshape(V, 3, :), reshape(F, 3, :); name=name)
end

"""
    write_obj(file, mesh)

Write a mesh as a Wavefront OBJ file.
"""
function write_obj(file::AbstractString, m::TriMesh)
    open(file, "w") do io
        println(io, "# ", isempty(m.name) ? "mesh" : m.name, " written by Occultations.jl (km)")
        for j in axes(m.vertices, 2)
            @printf(io, "v %.6f %.6f %.6f\n", m.vertices[:, j]...)
        end
        for j in axes(m.faces, 2)
            @printf(io, "f %d %d %d\n", m.faces[:, j]...)
        end
    end
    return file
end

"""
    ellipsoid_mesh(a, b=a, c=b; subdivisions=4, name) -> TriMesh
    ellipsoid_mesh(e::Ellipsoid; subdivisions=4, name)

Triangle mesh of a triaxial ellipsoid: an icosahedron subdivided `subdivisions` times
(20·4ⁿ faces), projected on the unit sphere and scaled by the semi-axes.
"""
function ellipsoid_mesh(a::Real, b::Real=a, c::Real=b; subdivisions::Integer=4,
                        name::AbstractString=@sprintf("ellipsoid(%.2f, %.2f, %.2f)", a, b, c))
    t = (1 + sqrt(5)) / 2
    V = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
         [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
    V = [v ./ norm(v) for v in V]
    F = [(1, 12, 6), (1, 6, 2), (1, 2, 8), (1, 8, 11), (1, 11, 12), (2, 6, 10), (6, 12, 5), (12, 11, 3),
         (11, 8, 7), (8, 2, 9), (4, 10, 5), (4, 5, 3), (4, 3, 7), (4, 7, 9), (4, 9, 10), (5, 10, 6),
         (3, 5, 12), (7, 3, 11), (9, 7, 8), (10, 9, 2)]
    for _ in 1:subdivisions
        cache = Dict{Tuple{Int,Int},Int}()
        mid(i, j) = get!(cache, minmax(i, j)) do
            m = (V[i] .+ V[j]) ./ 2
            push!(V, m ./ norm(m))
            length(V)
        end
        F = reduce(vcat, [[(i, mid(i, j), mid(k, i)), (j, mid(j, k), mid(i, j)), (k, mid(k, i), mid(j, k)),
                           (mid(i, j), mid(j, k), mid(k, i))] for (i, j, k) in F])
    end
    X = reduce(hcat, V) .* [a, b, c]
    return TriMesh(X, reduce(hcat, [collect(f) for f in F]); name=name)
end
ellipsoid_mesh(e::Ellipsoid; kw...) = ellipsoid_mesh(e.a, e.b, e.c; kw...)

"Volume of a closed mesh [km³] (divergence theorem)."
function volume(m::TriMesh)
    v = 0.0
    for j in axes(m.faces, 2)
        p, q, r = (view(m.vertices, :, m.faces[k, j]) for k in 1:3)
        v += dot(p, cross(q, r))
    end
    return abs(v) / 6
end

"Radius of the sphere with the same volume."
volume_equivalent_radius(s::AbstractShape) = cbrt(3volume(s) / 4π)

"Mesh scaled by `s` (a number or a 3-vector of axis factors)."
scale(m::TriMesh, s) = TriMesh(m.vertices .* s, m.faces; name=m.name)

# ---------------------------------------------------------------------------
# projection on the sky
# ---------------------------------------------------------------------------

# matrix whose rows are the sky basis (e, n) and the line of sight u (observer → body)
function _sky_matrix(u::AbstractVector)
    e, n = sky_basis(u)
    return vcat(e', n', (u ./ norm(u))')
end

"""
    project(shape, M, u; n=720) -> Outline

Apparent outline of a 3-D shape with orientation `M` (ICRS → body rotation) seen along the
direction `u` (observer → body, ICRS). An [`Ellipsoid`](@ref) gives an exact
[`EllipseOutline`](@ref); a [`TriMesh`](@ref) a [`PolarOutline`](@ref) sampled in `n`
directions (the mesh must be star-shaped as seen from its origin).
"""
function project(e::Ellipsoid, M::AbstractMatrix, u::AbstractVector; n::Integer=720)
    P = _sky_matrix(u)'                                  # sky → ICRS
    A = Diagonal([1 / e.a^2, 1 / e.b^2, 1 / e.c^2])
    Q = P' * M' * A * M * P
    # silhouette: Schur complement of the line-of-sight coordinate
    S11 = Q[1, 1] - Q[1, 3]^2 / Q[3, 3]
    S12 = Q[1, 2] - Q[1, 3] * Q[2, 3] / Q[3, 3]
    S22 = Q[2, 2] - Q[2, 3]^2 / Q[3, 3]
    tr2, dt = (S11 + S22) / 2, S11 * S22 - S12^2
    d = sqrt(max(tr2^2 - dt, 0.0))
    λ1, λ2 = tr2 - d, tr2 + d                           # λ1 ≤ λ2: major, minor axis
    # minor-axis eigenvector (λ2)
    vx, vy = abs(S12) > 1e-300 ? (S12, λ2 - S11) : (S11 >= S22 ? (1.0, 0.0) : (0.0, 1.0))
    return EllipseOutline(1 / sqrt(λ1), 1 / sqrt(λ2), mod(atand(vx, vy), 180.0))
end

function project(s::TriMesh, M::AbstractMatrix, u::AbstractVector; n::Integer=720)
    X = _sky_matrix(u) * (M' * s.vertices)
    F, G = view(X, 1, :), view(X, 2, :)
    dθ = 2π / n
    θ = collect(range(-π, π; length=n + 1)[1:n])
    r = zeros(n)
    @inbounds for (i, j) in s.edges
        pf, pg = F[i], G[i]
        qf, qg = F[j] - pf, G[j] - pg
        θp = atan(pg, pf)
        Δ = mod(atan(G[j], F[j]) - θp + π, 2π) - π
        lo, hi = minmax(θp, θp + Δ)
        for k in ceil(Int, (lo + π) / dθ):floor(Int, (hi + π) / dθ)
            m = mod1(k + 1, n)
            c, s_ = cos(θ[m]), sin(θ[m])
            den = c * qg - s_ * qf
            abs(den) < 1e-12 && continue
            ρ = (pf * qg - pg * qf) / den
            ρ > r[m] && (r[m] = ρ)
        end
    end
    return PolarOutline(θ, r)
end

"""
    project(body, t; observer=GEOCENTER, n=720) -> Outline

Apparent outline of the body's shape model at `t` (orientation from its frame at the
light-time corrected epoch).
"""
function project(b::Body, t; observer::AbstractObserver=GEOCENTER, n::Integer=720)
    b.shape === nothing && throw(ArgumentError("body $(b.name) has no shape model"))
    o = orientation(b, t; observer=observer)
    return project(b.shape, o.matrix, _topocentric_vector(b.ephem, observer, _epoch(t)); n=n)
end

"""
    projected_axes_dunham(a, b, c, phase, aspect) -> (a′, b′)

Semi-axes of the apparent ellipse of a triaxial ellipsoid (Dunham et al. 2019, appendix):
`phase` is the rotational phase from light-curve **minimum** and `aspect` the angle between
the rotation axis and the line of sight [deg].

# Examples
Haumea (1161 × 852 × 513 km) at an aspect of 76.3°, at light-curve minimum; the formula of
[`projected_axes_magnusson`](@ref) gives the same with the phase shifted by 90°:
```jldoctest
julia> round.(projected_axes_dunham(1161.0, 852.0, 513.0, 0.0, 76.3); digits=1)
(852.0, 569.2)

julia> all(projected_axes_dunham(1161.0, 852.0, 513.0, 30.0, 76.3) .≈
           projected_axes_magnusson(1161.0, 852.0, 513.0, 120.0, 76.3))
true
```
"""
function projected_axes_dunham(a, b, c, phase, aspect)
    ψ, φ = phase * DEG2RAD, aspect * DEG2RAD
    cψ2, sψ2, tφ2 = cos(ψ)^2, sin(ψ)^2, tan(φ)^2
    K = cψ2 / a^2 + sψ2 / b^2
    Δ = 1 + c^2 * tφ2 * K
    P = sψ2 * (1 + c^2 / b^2 * tφ2)^2 / (Δ^2 * a^2) + cψ2 * (1 + c^2 / a^2 * tφ2)^2 / (Δ^2 * b^2) +
        tφ2 * sψ2 * cψ2 * (c^2 / b^2 - c^2 / a^2)^2 / (Δ^2 * c^2)
    R = K / (Δ^2 * cos(φ)^2) + c^2 * tφ2 * K^2 / (Δ^2 * cos(φ)^2)
    Q = 2 * sin(ψ) * cos(ψ) / (Δ^2 * c^2 * cos(φ)) * (c^2 / b^2 - c^2 / a^2) * (1 + c^2 * tφ2 * K)
    θ = atan(Q, R - P) / 2
    ia2 = P * cos(θ)^2 + R * sin(θ)^2 - Q * sin(θ) * cos(θ)
    ib2 = P * sin(θ)^2 + R * cos(θ)^2 + Q * sin(θ) * cos(θ)
    x, y = 1 / sqrt(ia2), 1 / sqrt(ib2)
    return max(x, y), min(x, y)
end

"""
    projected_axes_magnusson(a, b, c, phase, aspect) -> (a′, b′)

Semi-axes of the apparent ellipse of a triaxial ellipsoid (Magnusson 1986, eqs. 25–26):
`phase` from light-curve **maximum**, `aspect` as in [`projected_axes_dunham`](@ref) [deg].
"""
function projected_axes_magnusson(a, b, c, phase, aspect)
    φ, ψ = phase * DEG2RAD, aspect * DEG2RAD
    A = b^2 * c^2 * sin(ψ)^2 * sin(φ)^2 + a^2 * c^2 * sin(ψ)^2 * cos(φ)^2 + a^2 * b^2 * cos(ψ)^2
    B = -(a^2 * (cos(ψ)^2 * sin(φ)^2 + cos(φ)^2) + b^2 * (cos(ψ)^2 * cos(φ)^2 + sin(φ)^2) + c^2 * sin(ψ)^2)
    d = sqrt(B^2 - 4A)
    return sqrt(2A / (-B - d)), sqrt(2A / (-B + d))
end

"""
    rotational_amplitude(a, b, c, aspect) -> mag

Peak-to-peak amplitude of the rotational light curve of a triaxial ellipsoid with
geometric scattering at zero phase angle, for the `aspect` angle between the rotation
axis and the line of sight [deg] (e.g. Binzel et al. 1989):
`2.5 log₁₀(a/b) - 1.25 log₁₀[(a² cos²θ + c² sin²θ)/(b² cos²θ + c² sin²θ)]`.

# Examples
```jldoctest
julia> rotational_amplitude(300, 200, 150, 90) ≈ 2.5 * log10(300 / 200)   # equator-on
true

julia> round(rotational_amplitude(300, 200, 150, 60); digits=3)
0.233
```
"""
function rotational_amplitude(a, b, c, aspect)
    cθ2, sθ2 = cosd(aspect)^2, sind(aspect)^2
    return 2.5 * log10(a / b) - 1.25 * log10((a^2 * cθ2 + c^2 * sθ2) / (b^2 * cθ2 + c^2 * sθ2))
end

"""
    rotational_amplitude(shape, pole_ra, pole_dec, u; nphase=72) -> mag

Peak-to-peak amplitude of the rotational light curve of any shape (geometric scattering:
brightness ∝ projected area) with the pole at (`pole_ra`, `pole_dec`) [deg] seen along `u`
(observer → body), from the extreme projected areas over a rotation.
"""
function rotational_amplitude(s::AbstractShape, pole_ra, pole_dec, u::AbstractVector; nphase::Integer=72)
    areas = [area(project(s, _R3(W) * _R1(90 - pole_dec) * _R3(90 + pole_ra), u; n=360))
             for W in range(0, 360; length=nphase + 1)[1:nphase]]
    return 2.5 * log10(maximum(areas) / minimum(areas))
end

"Aspect angle [deg] between the pole (`pole_ra`, `pole_dec`) and the line of sight `u`."
aspect_angle(pole_ra, pole_dec, u::AbstractVector) = acosd(clamp(dot(radec2vec(pole_ra, pole_dec), u ./ norm(u)), -1, 1))
