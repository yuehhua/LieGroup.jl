#
#
# Generic types for the interface of a Lie group
@doc """
    AbstractGroupOperation

Represent a type of group operation for a [`LieGroup`](@ref) ``$(_math(:G))``, that is a
smooth binary operation ``$(_math(:∘)) : $(_math(:G)) × $(_math(:G)) → $(_math(:G))``
on elements of a Lie group ``$(_math(:G))``.
"""
abstract type AbstractGroupOperation end

"""
    DefaultLieAlgebraOrthogonalBasis{𝔽} <: ManifoldsBase.AbstractOrthogonalBasis{𝔽,ManifoldsBase.TangentSpaceType}

Specify an orthogonal basis for a Lie algebra.
This is used as the default within [`hat`](@ref) and [`vee`](@ref).

If not specifically overwritten/implemented for a Lie group, the [`DefaultOrthogonalBasis`](@extref `ManifoldsBase.DefaultOrthogonalBasis`)
at the [`identity_element`](@ref) on the [`base_manifold](@ref base_manifold(::LieGroup)) acts as a fallback.

!!! note
    In order to implement the corresponding [`get_coordinates`](@ref) and [`get_vector`](@ref) functions,
    define `get_coordiinates_lie(::LieGroup, p, X, N)` and `get_vector_lie(::LieGroup, p, X, N)`, resp.
"""
struct DefaultLieAlgebraOrthogonalBasis{𝔽} <:
       ManifoldsBase.AbstractOrthogonalBasis{𝔽,ManifoldsBase.TangentSpaceType} end
function DefaultLieAlgebraOrthogonalBasis(𝔽::ManifoldsBase.AbstractNumbers=ℝ)
    return DefaultLieAlgebraOrthogonalBasis{𝔽}()
end

"""
    LieGroup{𝔽, O<:AbstractGroupOperation, M<:AbstractManifold{𝔽}} <: AbstractManifold{𝔽}

Represent a Lie Group ``$(_math(:G))``.

A *Lie Group* ``$(_math(:G))`` is a group endowed with the structure of a manifold such that the
group operations ``$(_math(:∘)): $(_math(:G))×$(_math(:G)) → $(_math(:G))``, see [`compose`](@ref)
and the inverse operation ``⋅^{-1}: $(_math(:G)) → $(_math(:G))``, see [`inv`](@ref) are smooth,
see for example [HilgertNeeb:2012; Definition 9.1.1](@cite).

Lie groups are named after the Norwegian mathematician [Marius Sophus Lie](https://en.wikipedia.org/wiki/Sophus_Lie) (1842–1899).

# Fields

* `manifold`: an $(_link(:AbstractManifold)) ``$(_math(:M))``
* `op`: an [`AbstractGroupOperation`](@ref) ``$(_math(:∘))`` on that manifold

# Constructor

    LieGroup(M::AbstractManifold, op::AbstractGroupOperation)

Generate a Lie group based on a manifold `M` and a group operation `op`, where vectors by default are stored in the Lie Algebra.
"""
struct LieGroup{𝔽,O<:AbstractGroupOperation,M<:ManifoldsBase.AbstractManifold{𝔽}} <:
       ManifoldsBase.AbstractManifold{𝔽}
    manifold::M
    op::O
end

@doc """
    Identity{O<:AbstractGroupOperation}

Represent the group identity element ``e ∈ $(_math(:G))`` on a [`LieGroup`](@ref) ``$(_math(:G))``
with [`AbstractGroupOperation`](@ref) of type `O`.

Similar to the philosophy that points are agnostic of their group at hand, the identity
does not store the group ``$(_math(:G))`` it belongs to. However it depends on the type of the [`AbstractGroupOperation`](@ref) used.

See also [`identity_element`](@ref) on how to obtain the corresponding [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`) or array representation.

# Constructors

    Identity(::LieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation}
    Identity(o::AbstractGroupOperation)
    Identity(::Type{AbstractGroupOperation})

create the identity of the corresponding subtype `O<:`[`AbstractGroupOperation`](@ref)
"""
struct Identity{O<:AbstractGroupOperation} end

Identity(::LieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation} = Identity{O}()
Identity(::O) where {O<:AbstractGroupOperation} = Identity(O)
Identity(::Type{O}) where {O<:AbstractGroupOperation} = Identity{O}()

#
#
# --- Functions ---

# Internal pass through for coordinates and vectors

@inline function ManifoldsBase._get_coordinates(
    G::LieGroup, p, X, B::DefaultLieAlgebraOrthogonalBasis
)
    return get_coordinates_lie(G, p, X, number_system(B))
end
@inline function ManifoldsBase._get_coordinates!(
    G::LieGroup, Y, p, X, B::DefaultLieAlgebraOrthogonalBasis
)
    return get_coordinates_lie!(G, Y, p, X, number_system(B))
end
@inline function ManifoldsBase._get_vector(
    G::LieGroup, p, c, B::DefaultLieAlgebraOrthogonalBasis
)
    return get_vector_lie(G, p, c, number_system(B))
end
@inline function ManifoldsBase._get_vector!(
    G::LieGroup, Y, p, c, B::DefaultLieAlgebraOrthogonalBasis
)
    return get_vector_lie!(G, Y, p, c, number_system(B))
end

_doc_adjoint = """
    adjoint(G::LieGroup, g X)
    adjoint!(G::LieGroup, Y, g, X)

Compute the adjoint ``$(_math(:Ad))(g): $(_math(:𝔤)) → $(_math(:𝔤))``, which is defined as
the differential [`diff_conjugate`](@ref) of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``
evaluated at the [`Identity`](@ref) ``h=$(_math(:e))``.
The operation can be performed in-place of `Y`.

```math
  $(_math(:Ad))(g)[X] = D c_g($(_math(:e))) [X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```

see [HilgertNeeb:2012; Section 9.2.3](@cite).

On matrix Lie groups the adjoint reads ``$(_math(:Ad))(g)[X] = g$(_math(:∘))X$(_math(:∘))g^{-1}``.
"""

@doc "$(_doc_adjoint)"
function Base.adjoint(G::LieGroup, g, X)
    Y = ManifoldsBase.allocate_result(G, adjoint, g, X)
    return adjoint!(G, Y, g, X)
end

function adjoint! end
@doc "$(_doc_adjoint)"
function adjoint!(G::LieGroup, Y, g, X)
    diff_conjugate!(G, Y, g, Identity(G), X)
    return Y
end

@doc """
    base_manifold(G::LieGroup)

Return the manifold stored within the [`LieGroup`](@ref) `G`.
"""
Manifolds.base_manifold(G::LieGroup) = G.manifold

function ManifoldsBase.check_point(
    G::LieGroup{𝔽,O}, g; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.check_point(G.manifold, g; kwargs...)
end
function ManifoldsBase.check_point(
    G::LieGroup{𝔽,O}, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return nothing
end
function ManifoldsBase.check_point(
    G::LieGroup{𝔽,O}, e::Identity{O2}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,O2<:AbstractGroupOperation}
    return DomainError(
        e,
        """
        The provided point $e is not the Identity on $G.
        Expected an Identity corresponding to $(G.op).
        """,
    )
end

ManifoldsBase.check_size(G::LieGroup, ::Identity) = nothing

function ManifoldsBase.check_vector(G::LieGroup, X; kwargs...)
    return ManifoldsBase.check_vector(G.manifold, identity_element(G), X; kwargs...)
end

# compose g ∘ h
_doc_compose = """
    compose(G::LieGroup, g, h)
    compose!(G::LieGroup, k, g, h)

Perform the group oepration ``g $(_math(:∘)) h`` for two ``g, h ∈ $(_math(:G))``
on the [`LieGroup`](@ref) `G`. This can also be done in-place of `h`.

!!! info
    This function also handles the case where `g` or/and `h` are the [`Identity`](@ref)`(G)`.
    Since this would lead to ambiguities when implementing a new group operations,
    this function calls `_compose` and `_compose!`, respectively, which is meant for the actual computation of
    group operations on (non-[`Identity`](@ref)` but maybe its numerical representation) elements.
"""
@doc "$(_doc_compose)"
compose(G::LieGroup, g, h) = _compose(G, g, h)
compose(::LieGroup{𝔽,O}, g::Identity{O}, h) where {𝔽,O<:AbstractGroupOperation} = h
compose(::LieGroup{𝔽,O}, g, h::Identity{O}) where {𝔽,O<:AbstractGroupOperation} = g
function compose(
    ::LieGroup{𝔽,O}, g::Identity{O}, h::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return g
end

function _compose(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose, g, h)
    return _compose!(G, k, g, h)
end

function compose! end

@doc "$(_doc_compose)"
compose!(G::LieGroup, k, g, h) = _compose!(G, k, g, h)
function compose!(G::LieGroup{𝔽,O}, k, ::Identity{O}, h) where {𝔽,O<:AbstractGroupOperation}
    return copyto!(G, k, h)
end
function compose!(G::LieGroup{𝔽,O}, k, g, ::Identity{O}) where {𝔽,O<:AbstractGroupOperation}
    return copyto!(G, k, g)
end
function compose!(
    G::LieGroup{𝔽,O}, k, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return identity_element!(G, k)
end
function compose!(
    ::LieGroup{𝔽,O}, k::Identity{O}, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return k
end

function _compose! end

_doc_conjugate = """
    conjugate(G::LieGroup, g, h)
    conjugate!(G::LieGroup, k, g, h)

Compute the conjugation map ``c_g: $(_math(:G)) → $(_math(:G))`` given by ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_conjugate)"
function conjugate(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, conjugate, h, g)
    return conjugate!(G, k, g, h)
end

function conjugate! end
@doc "$(_doc_conjugate)"
function conjugate!(G::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, k) # `h∘k` in-place of k
    compose!(G, k, g, k) # `g∘k` in-place of k
    return k
end

ManifoldsBase.copyto!(G::LieGroup, h, g) = copyto!(G.manifold, h, g)
function ManifoldsBase.copyto!(
    G::LieGroup{𝔽,O}, h, g::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.copyto!(G.manifold, h, identity_element(G))
end
function ManifoldsBase.copyto!(
    ::LieGroup{𝔽,O}, h::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return h
end
function ManifoldsBase.copyto!(
    G::LieGroup{𝔽,O}, h::Identity{O}, g
) where {𝔽,O<:AbstractGroupOperation}
    (is_identity(G, g)) && return h
    throw(
        DomainError(
            g,
            "copyto! into the identity element of $G ($h) is not defined for a non-identity element g ($g)",
        ),
    )
end
_doc_diff_conjugate = """
    diff_conjugate(G::LieGroup, g, h, X)
    diff_conjugate!(G::LieGroup, Y, g, h, X)

Compute the differential of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``,
which can be performed in-place of `Y`.

```math
  D(c_g(h))[X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```
"""
@doc "$(_doc_diff_conjugate)"
function diff_conjugate(G::LieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_conjugate, g, h, X)
    return diff_conjugate!(G, Y, g, h, X)
end

function diff_conjugate! end
@doc "$(_doc_diff_conjugate)"
diff_conjugate!(::LieGroup, Y, g, h, X)

_doc_diff_inv = """
    diff_inv(G::LieGroup, g, X)
    diff_inv!(G::LieGroup, Y, g, X)

Compute the differential of the function ``ι_{$(_math(:G))}(g) = g^{-1}``, where
``Dι_{$(_math(:G))}(g): $(_math(:𝔤)) → $(_math(:𝔤))``.
This can be done in-place of `Y`.
"""

@doc "$_doc_diff_inv"
function diff_inv(G::LieGroup, g, X)
    Y = allocate_result(G, diff_inv, g, X)
    return diff_inv!(G, Y, g, X)
end

function diff_inv! end
@doc "$_doc_diff_inv"
diff_inv!(G::LieGroup, Y, g, X)

_doc_diff_left_compose = """
    diff_left_compose(G::LieGroup, g, h, X)
    diff_left_compose!(G::LieGroup, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dλ_g(h)[X]``, ``X ∈ 𝔤``.
This can be done in-place of `Y`.
"""
@doc "$(_doc_diff_left_compose)"
function diff_left_compose(G::LieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_left_compose, g, h, X)
    return diff_left_compose!(G, Y, g, h, X)
end

function diff_left_compose! end
@doc "$(_doc_diff_left_compose)"
diff_left_compose!(::LieGroup, Y, g, h, X)

_doc_diff_right_compose = """
    diff_right_compose(G::LieGroup, h, g, X)
    diff_right_compose!(G::LieGroup, Y, h, g, X)

Compute the differential of the right group multiplication ``ρ_g(h) = h$(_math(:∘))g``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dρ_g(h)[X]``, ``X ∈ 𝔤``
This can be done in-place of `Y`.
"""
@doc "$(_doc_diff_right_compose)"
function diff_right_compose(G::LieGroup, h, g, X)
    Y = ManifoldsBase.allocate_result(G, diff_right_compose, h, g, X)
    return diff_right_compose!(G, Y, h, g, X)
end

function diff_right_compose! end
@doc "$(_doc_diff_right_compose)"
diff_right_compose!(::LieGroup, Y, g, h, X)

_doc_exp = """
    exp(G::LieGroup, g, X, t::Number=1)
    exp!(G::LieGroup, h, g, X, t::Number=1)

Compute the Lie group exponential map for ``g∈$(_math(:G))`` and ``X∈$(_math(:𝔤))``,
where ``$(_math(:𝔤))`` denotes the [`LieAlgebra`](@ref) of ``$(_math(:G))``.
It is given by

```math
$(_tex(:exp))_g X = g$(_math(:∘))$(_tex(:exp))_{$(_math(:G))}(X)
```

where `X` can be scaled by `t`, the computation can be performed in-place of `h`,
and ``$(_tex(:exp))_{$(_math(:G))}`` denotes the  [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)).

If `g` is the [`Identity`](@ref) the [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)) ``$(_tex(:exp))_{$(_math(:G))}`` is computed directly.
Implementing the Lie group exponential function introduces a default implementation with the formula above.

!!! note
    The Lie group exponential map is usually different from the exponential map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) exponential map, use `exp(`[`base_manifold`](@ref)`(G), g, X)`.
"""

@doc "$_doc_exp"
function ManifoldsBase.exp(G::LieGroup, g, X, t::Number=1)
    h = allocate_result(G, exp, g)
    exp!(G, h, g, X, t)
    return h
end

@doc "$_doc_exp"
function ManifoldsBase.exp!(G::LieGroup, h, g, X, t::Number=1)
    exp!(G, h, Identity(G), X, t)
    compose!(G, h, g, h)
    return h
end

_doc_exp_id = """
    exp(G::LieGroup, e::Identity, X)
    exp!(G::LieGroup, h, e::Identity, X)

Compute the (Lie group) exponential function

```math
$(_tex(:exp))_{$(_math(:G))}: $(_math(:𝔤)) → $(_math(:G)),$(_tex(:qquad)) $(_tex(:exp))_{$(_math(:G))}(X) = γ_X(1),
```

where ``γ_X`` is the unique solution of the initial value problem

```math
γ(0) = $(_math(:e)), $(_tex(:quad)) γ'(t) = γ(t)$(_math(:act))X.
```

The computation can be performed in-place of `h`.

See also [HilgertNeeb:2012; Definition 9.2.2](@cite).
"""

@doc "$(_doc_exp_id)"
function ManifoldsBase.exp(
    G::LieGroup{𝔽,O}, e::Identity{O}, X, t::Number=1
) where {𝔽,O<:AbstractGroupOperation}
    h = identity_element(G)
    exp!(G, h, e, X, t)
    return h
end

# Fallback to a MethodError to avoid stack overflow
@doc "$(_doc_exp_id)"
function ManifoldsBase.exp!(
    G::LieGroup{𝔽,Op}, h, e::Identity{Op}, X, t::Number=1
) where {𝔽,Op<:AbstractGroupOperation}
    throw(MethodError(exp!, (typeof(G), typeof(h), typeof(e), typeof(X), typeof(t))))
end

_doc_get_coordinates = """
    get_coordinates(G::LieGroup, g, X, B::AbstractBasis)
    get_coordinates(𝔤::LieAlgebra, X, B::AbstractBasis)
    get_coordinates!(G::LieGroup, c, g, X, B::AbstractBasis)
    get_coordinates!(𝔤::LieAlgebra, c, X, B::AbstractBasis)

Return the vector of coordinates to the decomposition of `X` with respect to an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`)
of the [`LieAlgebra`](@ref) `𝔤`.
Since all tangent vectors are assumed to be represented in the Lie algebra,
both signatures are equivalent.
The operation can be performed in-place of `c`.

By default this function requires [`identity_element`](@ref)`(G)` and calls
the corresponding [`get_coordinates`](@extref ManifoldsBase :jl:function:`ManifoldsBase.get_coordinates`) function
of the Riemannian manifold the Lie group is build on.

The inverse operation is [`get_vector`](@ref).

See also [`vee`](@ref).
"""

@doc "$(_doc_get_coordinates)"
ManifoldsBase.get_coordinates(G::LieGroup, g, X, B::ManifoldsBase.AbstractBasis)

@doc "$(_doc_exp_id)"
ManifoldsBase.get_coordinates!(G::LieGroup, c, g, X, B::ManifoldsBase.AbstractBasis)

function get_coordinates_lie(G::LieGroup, g, X, N)
    c = allocate_result(
        G, get_coordinates, identity_element(G), X, DefaultLieAlgebraOrthogonalBasis(N)
    )

    return get_coordinates_lie!(G, c, g, X, N)
end
function get_coordinates_lie!(G::LieGroup, c, g, X, N)
    return get_coordinates!(
        base_manifold(G), c, identity_element(G), X, ManifoldsBase.DefaultOrthogonalBasis(N)
    )
end

_doc_get_vector = """
    get_vector(G::LieGroup, g, c, B::AbstractBasis)
    get_vector(𝔤::LieAlgebra, c, B::AbstractBasis)
    get_vector!(G::LieGroup, X, g, c, B::AbstractBasis)
    get_vector!(𝔤::LieAlgebra, X, c, B::AbstractBasis)

Return the vector corresponding to a set of coefficients in an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`)
of the [`LieAlgebra`](@ref) `𝔤`.
Since all tangent vectors are assumed to be represented in the Lie algebra,
both signatures are equivalend.
The operation can be performed in-place of a tangent vector `X`.

By default this function requires [`identity_element`](@ref)`(G)` and calls
the corresponding [`get_vector`](@extref ManifoldsBase :jl:function:`ManifoldsBase.get_vectors`) function
of the Riemannian manifold the Lie group is build on.

The inverse operation is [`get_coordinates`](@ref).

See also [`hat`](@ref)
"""

@doc "$(_doc_get_vector)"
ManifoldsBase.get_vector(G::LieGroup, g, c, B::ManifoldsBase.AbstractBasis)

@doc "$(_doc_exp_id)"
ManifoldsBase.get_vector!(G::LieGroup, X, g, c, B::ManifoldsBase.AbstractBasis)

@inline function get_vector_lie(G::LieGroup, g, c, N)
    X = zero_vector(G, Identity(G))
    return get_vector_lie!(G, X, g, c, N)
end
@inline function get_vector_lie!(G::LieGroup, Y, g, c, N)
    return get_vector!(
        base_manifold(G), Y, identity_element(G), c, ManifoldsBase.DefaultOrthogonalBasis(N)
    )
end

_doc_hat = """
    hat(G::LieGroup, c)
    hat!(G::LieGroup, X, c)

Compute the hat map ``(⋅)^̂ `` that maps a vector of coordinates ``c_i``
with respect to a certain basis to a tangent vector in the Lie algebra

```math
X = $(_tex(:sum))_{i∈$(_tex(:Cal,"I"))} c_iB_i,
```

where ``$(_tex(:Set, "B_i"))_{i∈$(_tex(:Cal,"I"))}`` is a basis of the Lie algebra
and ``$(_tex(:Cal,"I"))`` a corresponding index set, which is usually ``$(_tex(:Cal,"I"))=$(_tex(:Set,raw"1,\ldots,n"))``.

The computation can be performed in-place of `X`.
The inverse of `hat` is [`vee`](@ref).

Technically, `hat` is a specific case of [`get_vector`](@ref) and is implemented using the
[`DefaultLieAlgebraOrthogonalBasis`](@ref)
"""

# function hat end
@doc "$(_doc_hat)"
function hat(G::LieGroup{𝔽}, c) where {𝔽}
    return get_vector_lie(G, Identity(G), c, 𝔽)
end

# function hat! end
@doc "$(_doc_hat)"
function hat!(G::LieGroup{𝔽}, X, c) where {𝔽}
    get_vector_lie!(G, X, Identity(G), c, 𝔽)
    return X
end

_doc_identity_element = """
    identity_element(G::LieGroup)
    identity_element!(G::LieGroup, e)

Return a point representation of the [`Identity`](@ref) on the [`LieGroup`](@ref) `G`.
By default this representation is the default array or number representation.
It should return the corresponding default representation of ``e`` as a point on `G` if
points are not represented by arrays.
This can be performed in-place of `e`.
"""
# `function identity_element end`
@doc "$(_doc_identity_element)"
function identity_element(G::LieGroup)
    e = ManifoldsBase.allocate_result(G, identity_element)
    return identity_element!(G, e)
end

function identity_element! end
@doc "$(_doc_identity_element)"
identity_element!(G::LieGroup, e)

_doc_inv = """
    inv(G::LieGroup, g)
    inv!(G::LieGroup, h, g)

Compute the inverse group element ``g^{-1}`` with respect to the [`AbstractGroupOperation`](@ref) ``$(_math(:∘))``
on the [`LieGroup`](@ref) ``$(_math(:G))``,
that is, return the unique element ``h=g^{-1}`` such that ``h$(_math(:∘))g=$(_math(:e))``, where ``$(_math(:e))`` denotes the [`Identity`](@ref).

This can be done in-place of `h`, without side effects, that is you can do `inv!(G, g, g)`.
"""

@doc "$_doc_inv"
function Base.inv(G::LieGroup, g)
    h = allocate_result(G, inv, g)
    return inv!(G, h, g)
end

function inv! end
@doc "$_doc_inv"
inv!(G::LieGroup, h, g)

function Base.inv(::LieGroup{𝔽,O}, e::Identity{O}) where {𝔽,O<:AbstractGroupOperation}
    return e
end

function inv!(G::LieGroup{𝔽,O}, g, ::Identity{O}) where {𝔽,O<:AbstractGroupOperation}
    return identity_element!(G, g)
end

_doc_inv_left_compose = """
    inv_left_compose(G::LieGroup, g, h)
    inv_left_compose!(G::LieGroup, k, g, h)

Compute the inverse of the left group operation ``λ_g(h) = g$(_math(:∘))h``,
on the [`LieGroup`](@ref) `G`, that is, compute ``λ_g^{-1}(h) = g^{-1}$(_math(:∘))h``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_left_compose)"
function inv_left_compose(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, inv_left_compose, g, h)
    return inv_left_compose!(G, k, g, h)
end

function inv_left_compose! end
@doc "$(_doc_compose)"
function inv_left_compose!(G::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, k, h) # compose `k∘h` in-place of k
    return k
end

_doc_inv_right_compose = """
    inv_right_compose(G::LieGroup, h, g)
    inv_right_compose!(G::LieGroup, k, h, g)

Compute the inverse of the right group operation ``ρ_g(h) = h$(_math(:∘))g``,
on the [`LieGroup`](@ref) `G`, that is compute ``ρ_g^{-1}(h) = h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_right_compose)"
function inv_right_compose(G::LieGroup, h, g)
    k = ManifoldsBase.allocate_result(G, inv_right_compose, h, g)
    return inv_right_compose!(G, k, h, g)
end

function inv_right_compose! end
@doc "$(_doc_inv_right_compose)"
function inv_right_compose!(G::LieGroup, k, h, g)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, k) # compose `h∘k` in-place of k
    return k
end

function is_identity end
@doc """
    is_identity(G::LieGroup, q; kwargs)

Check whether `q` is the identity on the [`LieGroup`](@ref) ``$(_math(:G))``.
This means it is either the [`Identity`](@ref)`{O}` with the respect to the corresponding
[`AbstractGroupOperation`](@ref) `O`, or (approximately) the correct point representation.

# See also

[`identity_element`](@ref), [`identity_element!`](@ref)
"""
is_identity(G::LieGroup, q)

function is_identity(G::LieGroup{𝔽,O}, h; kwargs...) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.isapprox(G, Identity{O}(), h; kwargs...)
end
function is_identity(
    ::LieGroup{𝔽,O}, ::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return true
end
# any other identity than the fitting one
function is_identity(
    G::LieGroup{𝔽,<:AbstractGroupOperation},
    h::Identity{<:AbstractGroupOperation};
    kwargs...,
) where {𝔽}
    return false
end

"""
    is_point(G::LieGroup, g; kwargs...)

Check whether `g` is a valid point on the Lie Group `G`.
This falls back to checking whether `g` is a valid point on `G.manifold`,
unless `g` is an [`Identity`](@ref). Then, it is checked whether it is the
identity element corresponding to `G`.
"""
ManifoldsBase.is_point(G::LieGroup, g; kwargs...)

_doc_is_vector = """
    is_vector(G::LieGroup, X; kwargs...)
    is_vector(G::LieGroup{𝔽,O}, e::Indentity{O}, X; kwargs...)

Check whether `X` is a tangent vector, that is an element of the [`LieAlgebra`](@ref)
of `G`.
The first variant calls [`is_point`](@extref `ManifoldsBase.is_point-Tuple{AbstractManifold, Any, Bool}`) on the [`LieAlgebra`](@ref) `𝔤` of `G`.
The second variant calls [`is_vector`](@extref ManifoldsBase.is_vector) on the $(_link(:AbstractManifold)) at the [`identity_element`](@ref).

All keyword arguments are passed on to the corresponding call
"""

@doc "$(_doc_is_vector)"
function ManifoldsBase.is_vector(G::LieGroup, X; kwargs...)
    return is_point(LieAlgebra(G), X; kwargs...)
end

@doc "$(_doc_is_vector)"
function ManifoldsBase.is_vector(
    G::LieGroup{𝔽,O}, e::Identity{O}, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return is_vector(G.manifold, identity_element(G), X; kwargs...)
end

"""
    isapprox(M::LieGroup, g, h; kwargs...)

Check if points `g` and `h` from [`LieGroup`](@ref) are approximately equal.
this function calls the corresponding $(_link(:isapprox)) on the $(_link(:AbstractManifold))
after handling the cases where one or more
of the points are the [`Identity`](@ref).
All keyword argments are passed to this function as well.
"""
ManifoldsBase.isapprox(G::LieGroup, g, h; kwargs...) = isapprox(G.manifold, g, h; kwargs...)
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g::Identity{O}, h; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.isapprox(G.manifold, identity_element(G), h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.isapprox(G.manifold, g, identity_element(G); kwargs...)
end
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g::Identity{O}, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return true
end
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g::Identity{O}, h::Identity{O2}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,O2<:AbstractGroupOperation}
    return false
end

_doc_jacobian_conjugate = """
    jacobian_conjugate(G::LieGroup, g, h, B::ManifoldsBase.AbstractBasis=DefaultLieAlgebraOrthogonalBasis())
    jacobian_conjugate!(G::LieGroup, J, g, h, B::ManifoldsBase.AbstractBasis=DefaultLieAlgebraOrthogonalBasis())

Compute the Jacobian of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``,
with respect to an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`).

This can be seen as a matrix representation of the [`diff_conjugate`](@ref) ``D(c_g(h))[X]``
with respect to the given basis.

!!! note
    For the case that `h` is the [`Identity`](@ref) and the relation of ``D(c_g(h))[X]``
    to the [`adjoint`](@ref) ``$(_math(:Ad))(g)``, the Jacobian then sometimes called “adjoint matrix”,
    e.g. in [SolaDerayAtchuthan:2021](@cite), when choosing as a basis the
    [`DefaultLieAlgebraOrthogonalBasis`](@ref)`()` that is used for [`hat`](@ref) and [`vee`](@ref).
"""
@doc "$(_doc_jacobian_conjugate)"
function jacobian_conjugate(
    G::LieGroup, g, h, B::ManifoldsBase.AbstractBasis=DefaultLieAlgebraOrthogonalBasis()
)
    J = ManifoldsBase.allocate_result(G, jacobian_conjugate, g, h, B)
    return jacobian_conjugate!(G, J, g, h, B)
end

function jacobian_conjugate! end
@doc "$(_doc_jacobian_conjugate)"
jacobian_conjugate!(
    ::LieGroup, J, g, h; B::ManifoldsBase.AbstractBasis=DefaultLieAlgebraOrthogonalBasis()
)

_doc_log = """
    log(G::LieGroup, g, h)
    log!(G::LieGroup, X, g, h)

Compute the Lie group logarithmic map ``$(_tex(:log))_g: $(_math(:G)) → $(_math(:𝔤))``,
where ``$(_math(:𝔤))`` denotes the [`LieAlgebra`](@ref) of ``$(_math(:G))``.
It is given by

```math
$(_tex(:log))_g h = $(_tex(:log))_{$(_math(:G))}(g^{-1}$(_math(:∘))h)
```

where ``$(_tex(:log))_{$(_math(:G))}`` denotes the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any))
The computation can be performed in-place of `X`.

If `g` is the [`Identity`](@ref) the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any)) ``$(_tex(:log))_{$(_math(:G))}`` is computed directly.
Implementing the Lie group logarithmic function introduces a default implementation for this function with the formula above.

!!! note
    The Lie group logarithmic map is usually different from the logarithmic map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) logarithmic map, use `log(`[`base_manifold`](@ref)`G, g, h)`.
"""

@doc "$_doc_log"
function ManifoldsBase.log(G::LieGroup, g, h)
    X = allocate_result(G, log, g, h)
    log!(G, X, g, h)
    return X
end
function ManifoldsBase.log(
    G::LieGroup{𝔽,O}, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return zero_vector(LieAlgebra(G))
end

@doc "$_doc_log"
function ManifoldsBase.log!(G::LieGroup, X, g, h)
    log!(G, X, Identity(G), compose(G, inv(G, g), h))
    return h
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,O}, X, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    zero_vector!(LieAlgebra(G), X)
    return X
end

_doc_log_id = """
    log(G::LieGroup, e::Identity, g)
    log!(G::LieGroup, X, e::Identity, g)

Compute the (Lie group) logarithmic function ``$(_tex(:log))_{$(_math(:G))}: $(_math(:G)) → $(_math(:𝔤))``,
which is the inverse of the [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any))
The computation can be performed in-place of `X`.
"""

@doc "$(_doc_log_id)"
function ManifoldsBase.log(G::LieGroup, e::Identity, g)
    X = allocate_result(G, log, g)
    log!(G, X, e, g)
    return X
end

ManifoldsBase.manifold_dimension(G::LieGroup) = manifold_dimension(G.manifold)

ManifoldsBase.project!(G::LieGroup, h, g) = project!(G.manifold, h, g)
# Since tangent vectors are always in the Lie algebra – do project always on TeG
function ManifoldsBase.project!(G::LieGroup, Y, g, X)
    return project!(G.manifold, Y, identity_element(G), X)
end

_doc_rand = """
    rand(::LieGroup; vector_at=nothing, σ::Real=1.0, kwargs...)
    rand(::LieAlgebra; σ::Real=1.0, kwargs...)
    rand!(::LieGroup, gX; vector_at=nothing, σ::Real=1.0, kwargs...)
    rand!(::LieAlgebra, X; σ::Real=1.0, kwargs...)

Compute a random point or tangent vector on a Lie group.

For points this just means to generate a random point on the
underlying manifold itself.

For tangent vectors, an element in the Lie Algebra is generated,
see also [`rand(::LieAlgebra; kwargs...)`](@ref)
"""

@doc "$(_doc_rand)"
Random.rand(::LieGroup; kwargs...)

@doc "$(_doc_rand)"
function Random.rand!(G::LieGroup, pX; kwargs...)
    return rand!(Random.default_rng(), G, pX; kwargs...)
end

function Random.rand!(rng::AbstractRNG, G::LieGroup, pX; vector_at=nothing, kwargs...)
    M = base_manifold(G)
    if vector_at === nothing # for points -> pass to manifold
        rand!(rng, M, pX, kwargs...)
    else # for tangent vectors -> materialize identity, pass to tangent space there.
        rand!(rng, M, pX; vector_at=identity_element(G), kwargs...)
    end
end

function ManifoldsBase.representation_size(G::LieGroup)
    return representation_size(G.manifold)
end

function Base.show(io::IO, G::LieGroup)
    return print(io, "LieGroup($(G.manifold), $(G.op))")
end

_doc_vee = """
    vee(G::LieGroup, X)
    vee!(G::LieGroup, c, X)

Compute the vee map ``(⋅)^∨`` that maps a tangent vector `X` from the [`LieAlgebra`](@ref)
to its coordinates with respect to the [`DefaultLieAlgebraOrthogonalBasis`](@ref) basis in the Lie algebra

```math
X = $(_tex(:sum))_{i∈$(_tex(:Cal,"I"))} c_iB_i,
```

where ``$(_tex(:Set, "B_i"))_{i∈$(_tex(:Cal,"I"))}`` is a basis of the Lie algebra
and ``$(_tex(:Cal,"I"))`` a corresponding index set, which is usually ``$(_tex(:Cal,"I"))=$(_tex(:Set,raw"1,\ldots,n"))``.

The computation can be performed in-place of `c`.
The inverse of `vee` is [`hat`](@ref).

Technically, `vee` is a specific case of [`get_coordinates`](@ref) and is implemented using
the [`DefaultLieAlgebraOrthogonalBasis`](@ref).
"""

# function vee end
@doc "$(_doc_vee)"
function vee(G::LieGroup{𝔽}, X) where {𝔽}
    return get_coordinates_lie(G, Identity(G), X, 𝔽)
end

# function vee! end
@doc "$(_doc_vee)"
function vee!(G::LieGroup{𝔽}, c, X) where {𝔽}
    get_coordinates_lie!(G, c, Identity(G), X, 𝔽)
    return c
end

function ManifoldsBase.zero_vector(
    G::LieGroup{𝔽,<:O}, ::Identity{<:O}
) where {𝔽,O<:AbstractGroupOperation}
    return zero_vector(G.manifold, identity_element(G))
end
function ManifoldsBase.zero_vector!(
    G::LieGroup{𝔽,<:O}, X, ::Identity{<:O}
) where {𝔽,O<:AbstractGroupOperation}
    return zero_vector!(G.manifold, X, identity_element(G))
end

#
# Allocation hints - mainly pass-through, especially for power manifolds
function ManifoldsBase.allocate_result(
    G::LieGroup,
    f::Union{typeof(compose),typeof(inv),typeof(conjugate),typeof(exp)},
    args...,
)
    return ManifoldsBase.allocate_result(G.manifold, ManifoldsBase.exp, args...)
end
function ManifoldsBase.allocate_result(G::LieGroup, f::typeof(log), args...)
    return ManifoldsBase.allocate_result(G.manifold, f, args...)
end
function ManifoldsBase.allocate_result(G::LieGroup, f::typeof(zero_vector), g)
    return ManifoldsBase.allocate_result(G.manifold, f, g)
end
function ManifoldsBase.allocate_result(
    G::LieGroup, f::Union{typeof(rand),typeof(identity_element)}
)
    # both get a type allocated like rand
    return ManifoldsBase.allocate_result(G.manifold, rand)
end
