module MeshViews

# A simple 3-D viewer for shape models (triangle meshes, ellipsoids) with Makie

using ..Occultations
using ..Occultations: TriMesh, Ellipsoid, AbstractShape, Body, AbstractObserver, GEOCENTER, _epoch, sky_basis,
                    orientation, topocentric_vector, project, outline_points, ellipsoid_mesh, isot
using LinearAlgebra
using Printf
using Makie
import ..Occultations: meshview

_mesh(m::TriMesh; kw...) = m
_mesh(e::Ellipsoid; subdivisions::Integer=4) = ellipsoid_mesh(e; subdivisions=subdivisions)

function _draw!(ax, V::AbstractMatrix, F::AbstractMatrix; color, wireframe, shading)
    mesh!(ax, permutedims(V), permutedims(F); color=color, shading=shading ? Makie.automatic : Makie.NoShading)
    if wireframe
        xs = Float64[]; ys = Float64[]; zs = Float64[]
        for j in axes(F, 2), (a, b) in ((1, 2), (2, 3), (3, 1))
            p, q = V[:, F[a, j]], V[:, F[b, j]]
            append!(xs, (p[1], q[1], NaN)); append!(ys, (p[2], q[2], NaN)); append!(zs, (p[3], q[3], NaN))
        end
        lines!(ax, xs, ys, zs; color=(:black, 0.25), linewidth=0.4)
    end
end

function _arrows!(ax, R, dirs, labels, colors)
    for (d, l, c) in zip(dirs, labels, colors)
        p = 1.35R .* d
        lines!(ax, [0, p[1]], [0, p[2]], [0, p[3]]; color=c, linewidth=3)
        text!(ax, p[1], p[2], p[3]; text=l, color=c, fontsize=16, font=:bold)
    end
end

"""
    meshview(shape; color=:gray75, wireframe=true, axes=true, subdivisions=4, title, size=(800, 750))

The shape (a [`TriMesh`](@ref) or an [`Ellipsoid`](@ref)) in its body frame: x (prime
meridian, red), y (green) and z (rotation axis, blue). With GLMakie the view can be rotated
and zoomed with the mouse.
"""
function meshview(s::AbstractShape; color=:gray75, wireframe::Bool=true, axes::Bool=true, subdivisions::Integer=4,
                  title::AbstractString=string(s), size=(800, 750), shading::Bool=true)
    m = _mesh(s; subdivisions=subdivisions)
    fig = Figure(; size=size)
    ax = Axis3(fig[1, 1]; aspect=:data, title=title, xlabel="x [km]", ylabel="y [km]", zlabel="z [km]")
    _draw!(ax, m.vertices, m.faces; color=color, wireframe=wireframe, shading=shading)
    R = maximum(sqrt.(sum(abs2, m.vertices; dims=1)))
    axes && _arrows!(ax, R, ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]), ("x", "y", "z (pole)"), (:red, :green, :blue))
    return fig
end

"""
    meshview(body, t; observer=GEOCENTER, color=:gray75, wireframe=false, size=(1300, 650))

The body's shape model at epoch `t` (orientation from its frame, light-time corrected):
left, in 3-D in sky coordinates (f towards east – plotted to the left as on the sky –, g
towards north, line of sight towards the observer), with the body axes; right, the projected limb as seen by the observer, with the
direction of the north pole.
"""
function meshview(b::Body, t; observer::AbstractObserver=GEOCENTER, color=:gray75, wireframe::Bool=false,
                  size=(1300, 650), subdivisions::Integer=4, shading::Bool=true)
    b.shape === nothing && throw(ArgumentError("body $(b.name) has no shape model"))
    t = _epoch(t)
    m = _mesh(b.shape; subdivisions=subdivisions)
    o = orientation(b, t; observer=observer)
    u = normalize(topocentric_vector(b, observer, t))
    e, n = sky_basis(u)
    S = vcat(e', n', -u')                         # body → sky: (f, g, towards the observer)
    V = S * (o.matrix' * m.vertices)
    fig = Figure(; size=size)
    ax = Axis3(fig[1, 1]; aspect=:data, title="$(b.name) – $(isot(t; digits=0)) UTC",
               xlabel="f (east) [km]", ylabel="g (north) [km]", zlabel="towards observer [km]",
               xreversed=true, azimuth=-0.5π, elevation=0.5π - 0.35)
    _draw!(ax, V, m.faces; color=color, wireframe=wireframe, shading=shading)
    R = maximum(sqrt.(sum(abs2, V; dims=1)))
    B = S * o.matrix'                             # body axes in sky coordinates
    _arrows!(ax, R, (B[:, 1], B[:, 2], B[:, 3]), ("x", "y", "pole"), (:red, :green, :blue))
    # projected limb as seen by the observer (east to the left, as on the sky)
    O = project(b.shape, o.matrix, u)
    ax2 = Axis(fig[1, 2]; aspect=DataAspect(), xreversed=true, xlabel="f (east) [km]", ylabel="g (north) [km]",
               title=@sprintf("seen from %s: sub-observer lat %.1f°, pole PA %.1f°", observer isa Occultations.Observer ?
                              observer.name : "observer", o.sub_observer_lat, o.pole_position_angle))
    x, y = outline_points(O)
    poly!(ax2, Point2f.(x, y); color=(:gray75, 1.0), strokecolor=:black, strokewidth=1)
    pp = B[:, 3]
    lines!(ax2, [0, 1.3R * pp[1]], [0, 1.3R * pp[2]]; color=:blue, linewidth=2)
    text!(ax2, 1.3R * pp[1], 1.3R * pp[2]; text=pp[3] > 0 ? "pole (towards us)" : "pole (away)", color=:blue)
    return fig
end

# rotation matrices about the world axes [deg]
_Rx(a) = (c = cosd(a); s = sind(a); [1.0 0 0; 0 c -s; 0 s c])
_Ry(a) = (c = cosd(a); s = sind(a); [c 0 s; 0 1.0 0; -s 0 c])
_Rz(a) = (c = cosd(a); s = sind(a); [c -s 0; s c 0; 0 0 1.0])


end # module
