# Build the manual:   julia --project=docs docs/make.jl
# (first time: julia --project=docs -e 'using Pkg; Pkg.instantiate()')
#
# Nothing slow is run here:
# * the example scripts are turned into pages with Literate, but not executed: the output shown
#   after each section of a script is the one stored in examples/output/<example>/output.txt
#   by the last run of the script;
# * every figure referenced as assets/figures/<file> in docs/src is copied from
#   examples/output/*/.

using Documenter, DocumenterCitations, Literate, Occultations

const ROOT = dirname(@__DIR__)
const SRC = joinpath(@__DIR__, "src")

const REPO = "https://github.com/mkretlow/Occultations.jl/blob/main/"

# the examples, in the order of the manual: synthetic, then the published events grouped by the
# kind of analysis (as in docs/src/examples.md and the paper)
const SYNTHETIC = ["tutorial", "sim_binary", "sim_contact_binary", "sim_contact_binary_3d", "sim_rings", "sim_ring_pole"]
const REPRODUCED = ["makemake_2011", "varda_2018", "huya_2019", "2003vs2_2019",      # one event
                    "2002ms4_2020",                                                  # a large campaign
                    "bienor_2019_2023", "ixion_2020_2023",                           # several events
                    "2003uy117_2020",                                                # three-dimensional shape
                    "haumea_2017", "chariklo_2017_2020", "quaoar_2011_2022",         # rings
                    "2007oc10_2022",                                                 # double star
                    "2002xv93_2024"]                                                 # atmosphere
const DIFFERENCES = ["2007uk126_2014", "2002gz32_2017", "2002tc302_2018", "chiron_2018_2019", "1995sm55_2024",  # few chords
                     "huya_2023", "huya_2023_satellite",                             # negative chords, satellite
                     "2003az84_2012_2014", "2002kx14_2020_2023"]                     # several events
const EXAMPLES = [SYNTHETIC; REPRODUCED; DIFFERENCES]

# the printed output of an example, by section (the titles printed by `section("…")`)
function example_output(name)
    text = read(joinpath(ROOT, "examples", "output", name, "output.txt"), String)
    parts = split(text, r"\n?={78}\n(.*)\n={78}\n")
    titles = [m.captures[1] for m in eachmatch(r"={78}\n(.*)\n={78}\n", text)]
    out = Dict{String,Vector{String}}()
    for (t, body) in zip(titles, parts[2:end])
        lines = filter(l -> !startswith(l, "figures written to"), split(rstrip(body), '\n'))
        out[t] = String.(lines)
    end
    return out
end

# Literate preprocessing of an example script: a note on where the output comes from after the
# introduction; each `section("…")` becomes a heading, and the stored output of the section is
# inserted at its end
function insert_output(src, name)
    out = example_output(name)
    block(t) = haskey(out, t) ? vcat("", "# ```text", ["# " * l for l in out[t]], "# ```", "") :
               error("examples/$name.jl: no output stored for section \"$t\"")
    note = ["# !!! note \"Output\"",
            "#     The code is not run when the manual is built. The output after each section is that of the",
            "#     last run of [`examples/$name.jl`]($(REPO)examples/$name.jl), stored in",
            "#     [`examples/output/$name/output.txt`]($(REPO)examples/output/$name/output.txt).", ""]
    res, cur, noted = String[], nothing, false
    for l in split(src, '\n')
        if !noted && !isempty(strip(l)) && !startswith(l, "#")
            append!(res, note); noted = true
        end
        m = match(r"^section\(\"(.*)\"\)\s*$", l)
        if m === nothing
            push!(res, l)
        else
            cur === nothing || append!(res, block(cur))
            cur = replace(m.captures[1], "\\\"" => "\"")
            append!(res, ["# ## " * cur, ""])
        end
    end
    cur === nothing || append!(res, block(cur))
    return join(res, '\n')
end

function literate_examples()
    dir = joinpath(SRC, "generated")
    rm(dir; force=true, recursive=true)
    for name in EXAMPLES
        Literate.markdown(joinpath(ROOT, "examples", name * ".jl"), dir; flavor=Literate.DocumenterFlavor(),
                          codefence="```julia" => "```", credit=false, preprocess=s -> insert_output(s, name))
    end
end
literate_examples()

# copy the figures used by the pages from the example outputs
function copy_figures()
    dest = mkpath(joinpath(SRC, "assets", "figures"))
    used = Set{String}()
    for (dir, _, files) in walkdir(SRC), f in files
        endswith(f, ".md") || continue
        for m in eachmatch(r"assets/figures/([^)\s]+)", read(joinpath(dir, f), String))
            push!(used, m.captures[1])
        end
    end
    outputs = joinpath(ROOT, "examples", "output")
    for name in sort!(collect(used))
        src = [joinpath(outputs, d, name) for d in readdir(outputs) if isfile(joinpath(outputs, d, name))]
        isempty(src) && error("figure $name not found in examples/output/*/")
        cp(first(src), joinpath(dest, name); force=true)
    end
    return length(used)
end
copy_figures()

const BIB = CitationBibliography(joinpath(SRC, "refs.bib"); style=:authoryear)

DocMeta.setdocmeta!(Occultations, :DocTestSetup, :(using Occultations); recursive=true)

makedocs(;
    sitename = "Occultations.jl",
    authors = "Mike Kretlow",
    modules = [Occultations],
    repo = Remotes.GitHub("mkretlow", "Occultations.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://mkretlow.github.io/Occultations.jl",
        assets = ["assets/citations.css"],
        edit_link = "main",
        size_threshold = 400 * 2^10,         # the comparison and API pages are long
        size_threshold_warn = 300 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "User guide" => [
            "guide/quickstart.md",
            "guide/lightcurves.md",
            "guide/limbfits.md",
            "guide/binaries.md",
            "guide/rings.md",
            "guide/atmospheres.md",
            "guide/plots.md",
            "guide/occview.md",
            "guide/design.md",
        ],
        "Methods" => [
            "theory/limbfit.md",
            "theory/residuals.md",
            "theory/sigma_limb.md",
            "theory/uncertainty_bands.md",
        ],
        "Examples" => [
            "Overview" => "examples.md",
            "Synthetic events" => ["generated/$name.md" for name in SYNTHETIC],
            "Published: reproduced" => ["generated/$name.md" for name in REPRODUCED],
            "Published: with differences" => ["generated/$name.md" for name in DIFFERENCES],
            "comparison/publications.md",
            "comparison/sora.md",
        ],
        "validation.md",
        "api.md",
        "references.md",
    ],
    plugins = [BIB],
    checkdocs = :exports,
    doctest = true,
    linkcheck = get(ENV, "LINKCHECK", "false") == "true",
    # the repository on GitHub is replaced only after this version is pushed
    linkcheck_ignore = [r"^https://github\.com/mkretlow/Occultations\.jl"],
)

deploydocs(;
    repo = "github.com/mkretlow/Occultations.jl.git",
    devbranch = "main",
    push_preview = false,
)
