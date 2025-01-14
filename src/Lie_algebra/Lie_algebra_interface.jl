
"""
    LieAlgebra{𝔽, G} <: AbstractManifold{𝔽}

Represent the Lie algebra ``$(_math(:𝔤))``, that is a ``𝔽`` vector space with an associated
[`lie_bracket`](@ref) ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The Lie algebras considered here are those related to a [`LieGroup`](@ref) ``$(_math(:G))``,
namely the tangent space ``T_{$(_math(:e))}$(_math(:G))`` at the [`Identity`](@ref),
this is internally just a `const` of the corresponding $(_link(:TangentSpace)).



!!! note "Convention for representing tangent vectors in the Lie algebra"
    A vector field ``$(_tex(:Cal,"X")): $(_math(:G)) → T$(_math(:G))``, ``X(g) ∈ T_g$(_math(:G))``
    is called a left-invariant vector field if it satisfies

    ```math
    $(_tex(:Cal,"X"))(λ_g(h)) = Dλ_g(h)[$(_tex(:Cal,"X"))(h)], $(_tex(:quad))$(_tex(:text, "for all"))$(_tex(:quad)) g, h ∈ $(_math(:G)),
    ```

    where ``λ_g: $(_math(:G)) → $(_math(:G))`` is the left multiplication by ``g``.
    Hence ``$(_tex(:Cal,"X"))`` is determined already when ``X ∈ $(_math(:𝔤))`` is given,
    since ``$(_tex(:Cal,"X"))(g) = Dλ_g(e)[X]``, cf [HilgertNeeb:2012; Definition 9.1.7](@cite).

    Throughout `LieGroups.jl`, we use this left-invariant convention to store tangent vectors at points on a Lie group as elements of the corresponding Lie algebra.

# Constructor

    LieAlgebra(G::LieGroup)

Return the Lie Algebra belonging to the [`LieGroup`](@ref) `G`.
"""
const LieAlgebra{𝔽,O<:AbstractGroupOperation,G<:LieGroup{𝔽,O}} = ManifoldsBase.Fiber{
    𝔽,ManifoldsBase.TangentSpaceType,G,Identity{O}
}

function LieAlgebra(G::LieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation}
    return LieAlgebra{𝔽,O,typeof(G)}(G, Identity(G), ManifoldsBase.TangentSpaceType())
end

function ManifoldsBase.get_coordinates(𝔤::LieAlgebra, X, B::ManifoldsBase.AbstractBasis)
    G = 𝔤.manifold
    return get_coordinates(base_manifold(G), identity_element(G), X, B)
end
function ManifoldsBase.get_coordinates!(𝔤::LieAlgebra, c, X, B::ManifoldsBase.AbstractBasis)
    G = 𝔤.manifold
    get_coordinates!(base_manifold(G), c, identity_element(G), X, B)
    return c
end

function ManifoldsBase.get_vector(𝔤::LieAlgebra, c, B::ManifoldsBase.AbstractBasis)
    G = 𝔤.manifold
    return get_vector(base_manifold(G), identity_element(G), c, B)
end
function ManifoldsBase.get_vector!(𝔤::LieAlgebra, X, c, B::ManifoldsBase.AbstractBasis)
    G = 𝔤.manifold
    get_vector!(base_manifold(G), X, identity_element(G), c, B)
    return X
end

"""
    is_point(𝔤::LieAlgebra, X; kwargs...)

Check whether `X` is a valid point on the Lie Algebra `𝔤`.
This falls back to checking whether `X` is a valid point on the tangent space
at the [`identity_element`](@ref)`(G)` on `G.manifold` on the [`LieGroup`](@ref)
of `G`
"""
function ManifoldsBase.is_point(𝔤::LieAlgebra, X; kwargs...)
    # the manifold stored in the Fiber / Lie algebra is the Lie group G
    G = 𝔤.manifold
    e = identity_element(G)
    return ManifoldsBase.is_vector(G.manifold, e, X; kwargs...)
end

_doc_lie_bracket = """
    lie_bracket!(𝔤::LieAlgebra, X, Y)
    lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The computation can be done in-place of `Z`.
"""
function lie_bracket end
@doc "$(_doc_lie_bracket)"
function lie_bracket(𝔤::LieAlgebra, X, Y)
    Z = ManifoldsBase.allocate_result(𝔤, lie_bracket, X, Y)
    return lie_bracket!(𝔤, Z, X, Y)
end

function lie_bracket! end
@doc "$(_doc_lie_bracket)"
lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

function LinearAlgebra.norm(
    G::LieGroup{𝔽,O}, ::Identity{O}, X
) where {𝔽,O<:AbstractGroupOperation}
    return LinearAlgebra.norm(G, identity_element(G), X)
end

_doc_rand_algebra = """
    rand(::LieGroup; vector_at=nothing, σ=1.0, kwargs...)
    rand(::LieAlgebra; σ=1.0, kwargs...)
    rand!(::LieGroup, gX; vector_at=nothing, kwargs...)
    rand!(::LieAlgebra, X; σ=1.0, kwargs...)

Compute a random point or tangent vector on a Lie group.

For points this just means to generate a random point on the
underlying manifold itself.

For tangent vectors, an element in the Lie Algebra is generated,
see also [`rand(::LieAlgebra; kwargs...)`](@ref)
"""

@doc "$(_doc_rand_algebra)"
Random.rand(::LieAlgebra; kwargs...)

@doc "$(_doc_rand_algebra)"
Random.rand!(::LieAlgebra, X; kwargs...)

function Base.show(io::IO, 𝔤::LieAlgebra)
    return print(io, "LieAlgebra( $(𝔤.manifold) )")
end

function ManifoldsBase.zero_vector(𝔤::LieAlgebra)
    return ManifoldsBase.zero_vector(𝔤.manifold.manifold, identity_element(𝔤.manifold))
end

function ManifoldsBase.zero_vector!(𝔤::LieAlgebra, X)
    return ManifoldsBase.zero_vector!(𝔤.manifold.manifold, X, identity_element(𝔤.manifold))
end
