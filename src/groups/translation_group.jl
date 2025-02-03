"""
    TranslationGroup{𝔽,T}

The translation group ``$(_math(:T))(n)`` is Lie group consisting of
the [`AdditionGroupOperation`](@ref) on some [`Euclidean`](@extref `Manifolds.Euclidean`) space.

# Constructor
    TranslationGroup(n₁,...,nᵢ; kwargs...)

Generate the translation group on ``𝔽^{n₁,…,nᵢ}`` = `Euclidean(n₁,...,nᵢ; field=𝔽)`,
which is isomorphic to the group itself. All keyword arguments in `kwargs...`
are passed on to [`Euclidean`](@extref `Manifolds.Euclidean`) as well

We denote the Lie Algebra of ``$(_math(:T))(n)`` by ``$(_math(:t))(n)``.
"""
const TranslationGroup{𝔽,T} = LieGroup{𝔽,AdditionGroupOperation,Manifolds.Euclidean{T,𝔽}}

function TranslationGroup(n::Int...; kwargs...)
    Rn = Manifolds.Euclidean(n...; kwargs...)
    return TranslationGroup{typeof(Rn).parameters[[2, 1]]...}(Rn, AdditionGroupOperation())
end

function Base.show(io::IO, G::TranslationGroup{𝔽}) where {𝔽}
    size = Manifolds.get_parameter(G.manifold.size)
    return print(io, "TranslationGroup($(join(size, ", ")); field=$(𝔽))")
end
