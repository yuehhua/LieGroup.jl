#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
docs/make.jl

Render the `LieGroups.jl` documentation with optional arguments

Arguments
* `--exclude-tutorials` - exclude the tutorials from the menu of Documenter,
  this can be used if you do not have Quarto installed to still be able to render the docs
  locally on this machine. This option should not be set on CI.
* `--help`              - print this help and exit without rendering the documentation
* `--prettyurls`        – toggle the prettyurls part to true (which is otherwise only true on CI)
* `--quarto`            – run the Quarto notebooks from the `tutorials/` folder before generating the documentation
  this has to be run locally at least once for the `tutorials/*.md` files to exist that are included in
  the documentation (see `--exclude-tutorials`) for the alternative.
  If they are generated once they are cached accordingly.
  Then you can spare time in the rendering by not passing this argument.
  If quarto is not run, some tutorials are generated as empty files, since they
  are referenced from within the documentation. These are currently
  `Optimize.md` and `ImplementOwnManifold.md`.
""",
    )
    exit(0)
end

#
# (a) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

# (b) Did someone say render?
if "--quarto" ∈ ARGS
    using CondaPkg
    CondaPkg.withenv() do
        @info "Rendering Quarto"
        tutorials_folder = (@__DIR__) * "/../tutorials"
        # instantiate the tutorials environment if necessary
        Pkg.activate(tutorials_folder)
        # For a breaking release -> also set the tutorials folder to the most recent version
        Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))
        Pkg.resolve()
        Pkg.instantiate()
        Pkg.build("IJulia") # build `IJulia` to the right version.
        Pkg.activate(@__DIR__) # but return to the docs one before
        run(`quarto render $(tutorials_folder)`)
        return nothing
    end
else # fallback to at least create empty files for Optimize and Implement
    #    touch(joinpath(@__DIR__, "src/tutorials/Optimize.md"))
end

tutorials_in_menu = true
if "--exclude-tutorials" ∈ ARGS
    @warn """
    You are excluding the tutorials from the Menu,
    which might be done if you can not render them locally.

    Remember that this should never be done on CI for the full documentation.
    """
    tutorials_in_menu = false
end

# (c) load necessary packages for the docs
using Documenter
using DocumenterCitations, DocumenterInterLinks
using LieGroups

# (d) add contributing.md to docs
generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/JuliaManifolds/LieGroups.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)
for (md_file, doc_file) in [("CONTRIBUTING.md", "contributing.md"), ("NEWS.md", "news.md")]
    open(joinpath(generated_path, doc_file), "w") do io
        # Point to source license file
        println(
            io,
            """
            ```@meta
            EditURL = "$(base_url)$(md_file)"
            ```
            """,
        )
        # Write the contents out below the meta block
        for line in eachline(joinpath(dirname(@__DIR__), md_file))
            println(io, line)
        end
    end
end

## Build tutorials menu
tutorials_menu =
    "How to..." => [
        "🚀 Get Started with LieGroups.jl" => "tutorials/getstarted.md",
        "Transition from `GroupManifolds`" => "tutorials/transition.md",
    ]
# (e) finally make docs
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style=:alpha)
links = InterLinks(
    "ManifoldsBase" => ("https://juliamanifolds.github.io/ManifoldsBase.jl/stable/"),
    "Manifolds" => ("https://juliamanifolds.github.io/Manifolds.jl/stable/"),
)
makedocs(;
    format=Documenter.HTML(;
        prettyurls=(get(ENV, "CI", nothing) == "true") || ("--prettyurls" ∈ ARGS),
        assets=["assets/favicon.ico", "assets/citations.css"],
    ),
    modules=[LieGroups],
    authors="Seth Axen, Mateusz Baran, Ronny Bergmann, Olivier Verdier, and contributors",
    sitename="LieGroups.jl",
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        (tutorials_in_menu ? [tutorials_menu] : [])...,
        "Lie groups" => [
            "List of Lie groups" => "groups/index.md",
            "General Linear" => "groups/general_linear.md",
            "Heisenberg" => "groups/heisenberg_group.md",
            "Power group" => "groups/power_group.md",
            "Product group" => "groups/product_group.md",
            "Semidirect product group" => "groups/semidirect_product_group.md",
            "Special Euclidean group" => "groups/special_euclidean.md",
            "Special orthogonal group" => "groups/special_orthogonal.md",
            "Special unitary group" => "groups/special_unitary.md",
            "Translation group" => "groups/translation.md",
            "unitary group" => "groups/unitary.md",
        ],
        "Interfaces" => [
            "Lie group" => "interface/group.md",
            "Lie algebra" => "interface/algebra.md",
            "Group operation" => "interface/operations.md",
            "Group action" => "interface/actions.md",
        ],
        "Contributing to LieGroups.jl" => "contributing.md",
        "Notation" => "notation.md",
        "Changelog" => "news.md",
        "References" => "references.md",
    ],
    plugins=[bib, links],
)
deploydocs(; repo="github.com/JuliaManifolds/LieGroups.jl", push_preview=true)
#back to main env
Pkg.activate()
