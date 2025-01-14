# for (g, t)
const LeftSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:LeftSemidirectProductGroupOperation{
        <:MatrixMultiplicationGroupOperation,
        <:AdditionGroupOperation,
        LeftGroupOperationAction,
    },
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Rotations{T},<:Manifolds.Euclidean{T,ℝ}}
    },
}

# for (t, g)
const RightSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:RightSemidirectProductGroupOperation{
        <:AdditionGroupOperation,
        <:MatrixMultiplicationGroupOperation,
        LeftGroupOperationAction,
    },
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Euclidean{T,ℝ},<:Manifolds.Rotations{T}}
    },
}

"""
    SpecialEuclideanGroup{T}

The special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``. is the Lie group consisting of the
[`LeftSemidirectProductGroupOperation`](@ref) of the [`SpecialOrthogonalGroup`](@ref) and the
[`TranslationGroup`](@ref) together with the [`GroupOperationAction`](@ref)`{`[`LeftGroupOperationAction`](@ref)`}`.

To be precise, the group operation is defined on ``$(_math(:SO))(n) ⋉ $(_math(:T))(n)`` as follows:

```math
(r_1, t_1) ⋅ (r_2, t_2) = (r_1$(_math(:∘))r_2, t_1 + r_1$(_math(:⋅))t_2),
```

where ``r_1,r_2 ∈ $(_math(:SO))(n)`` and ``t_1,t_2 ∈ $(_math(:T))(n)``

Analogously you can write this on elements of ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)`` as

```math
(t_1, r_1) ⋅ (t_2, r_2) = (t_1 + r_1$(_math(:⋅))s_2, r_1$(_math(:∘))r_2)
```

Both these cases can be represented in a single matrix in [affine form](https://en.wikipedia.org/wiki/Affine_group#Matrix_representation)

```math
g = $(_tex(:pmatrix, "r & t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")),
$(_tex(:qquad)) r ∈ $(_math(:SO))(n), t ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

We refer also in general to elements on ``$(_math(:SE))(n)`` as ``g``
and their rotation and translation components as ``r`` and ``t``, respectively.

# Constructor
    SpecialEuclideanGroup(n; variant=:left, kwargs...)
    SpecialOrthogonalGroup(n; kwargs...) ⋉ TranslationGroup(n; kwargs...)
    TranslationGroup(n; kwargs...) ⋊ SpecialOrthogonalGroup(n; kwargs...)

Generate special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``, where the first
constructor is equivalent to the second.

All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.

The default representation for ``$(_math(:SE))(n)`` is the affine form. Alternatively
you can use the `ArrayPartition` from [`RecursiveArrayTools.jl`](https://docs.sciml.ai/RecursiveArrayTools/stable/) to work on ``(r,t)``.
or for ``$(_math(:T))(n) ⋊ $(_math(:SO))(n)`` using the `ArrayPartition`s ``(t,r)``;
which corresponds to setting `variant=:right` in the first constructor.
"""
const SpecialEuclideanGroup{T} = Union{
    <:LeftSpecialEuclideanGroup{T},<:RightSpecialEuclideanGroup{T}
}

const SpecialEuclideanOperation = Union{
    <:LeftSemidirectProductGroupOperation{
        <:MatrixMultiplicationGroupOperation,
        <:AdditionGroupOperation,
        LeftGroupOperationAction,
    },
    <:RightSemidirectProductGroupOperation{
        <:AdditionGroupOperation,
        <:MatrixMultiplicationGroupOperation,
        LeftGroupOperationAction,
    },
}

# This union we can also use for the matrix case where we do not care

function SpecialEuclideanGroup(n; variant=:left, kwargs...)
    SOn = SpecialOrthogonalGroup(n; kwargs...)
    Tn = TranslationGroup(n; kwargs...)
    variant ∉ [:left, :right] && error(
        "SE(n) requires a  variant ∉ [:left, :right] but you provided `variant=:$variant`",
    )
    return variant === :left ? SOn ⋉ Tn : Tn ⋊ SOn
end

"""
    default_left_action(G::SpecialOrthogonalGroup, ::TranslationGroup)

Return the default left action for the special Euclidean group ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)``,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
default_left_action(::SpecialOrthogonalGroup, ::TranslationGroup) =
    LeftGroupOperationAction()

"""
    default_right_action(::TranslationGroup, G::SpecialOrthogonalGroup)

Return the default right action for the special Euclidean group,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
function default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)
    return LeftGroupOperationAction()
end

function _check_matrix_affine(p, n; v=1, kwargs...)
    if !isapprox(p[end, :], [zeros(size(p, 2) - 1)..., v]; kwargs...)
        return nothing
        DomainError(p[end, :], "The last row of $p is not of form [0,..,0,$v].")
    end
    return nothing
end
# Order in a unified way
function ManifoldsBase.check_point(G::LeftSpecialEuclideanGroup, p; kwargs...)
    return _check_point(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], p; kwargs...)
end
function ManifoldsBase.check_point(G::RightSpecialEuclideanGroup, p; kwargs...)
    return _check_point(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], p; kwargs...)
end

function _check_point(
    G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, p; kwargs...
) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(p, n; v=1, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SOn
    errS = ManifoldsBase.check_point(
        Rotn, ManifoldsBase.submanifold_component(G, p, :Rotation); kwargs...
    )
    !isnothing(errS) && push!(errs, errS)
    # translate part
    errT = ManifoldsBase.check_point(
        Rn, ManifoldsBase.submanifold_component(G, p, :Translation); kwargs...
    )
    !isnothing(errT) && push!(errs, errT)
    if length(errs) > 1
        return ManifoldsBase.CompositeManifoldError(errs)
    end
    return length(errs) == 0 ? nothing : first(errs)
end

# Order in a unified way
function ManifoldsBase.check_vector(G::LeftSpecialEuclideanGroup, X; kwargs...)
    return _check_vector(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], X; kwargs...)
end
function ManifoldsBase.check_vector(G::RightSpecialEuclideanGroup, X; kwargs...)
    return _check_vector(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], X; kwargs...)
end

function _check_vector(
    G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, X; kwargs...
) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(p, n; v=0, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SO(n)  part
    SOn = LieGroup(Rotn, op1)
    errS = ManifoldsBase.check_vector(
        SOn, ManifoldsBase.submanifold_component(G, X, :Rotation); kwargs...
    )
    !isnothing(errS) && push!(errs, errS)
    # T(n) part
    Tn = LieGroup(Rn, G.op2)
    errT = ManifoldsBase.check_vector(
        Tn, ManifoldsBase.submanifold_component(G, X, :Translation); kwargs...
    )
    !isnothing(errT) && push!(errs, errT)
    (length(errs) > 1) && (return ManifoldsBase.CompositeManifoldError(errs))
    return length(errs) == 0 ? nothing : first(errs)
end

function _compose!(
    ::SpecialEuclideanGroup, k::AbstractMatrix, g::AbstractMatrix, h::AbstractMatrix
)
    mul!(k, g, h)
    return k
end

#
#
# We overwrite the `submanifold_component` a bit different, since they are not ordered
# as (1) and (2) but more semantically as :Rotation and :Translation, we also use that here

@inline function ManifoldsBase.submanifold_component(G::SpecialEuclideanGroup, g, s::Symbol)
    return ManifoldsBase.submanifold_component(G, g, Val(s))
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{:Rotation}
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p, 1:n, n + 1)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{:Translation}
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p, 1:n, 1:n)
end

_doc_inv_SEn = """
    inv(G::SpecialEuclideanGroup, g)
    inv(G::SpecialEuclideanGroup, h, g)

Compute the inverse element of a ``g ∈ $(_math(:SE))(n)`` from the [`SpecialEuclideanGroup`](@ref)`(n)`.

In affine form, ``g = $(_tex(:pmatrix, "r & t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1"))``,
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

The inverse reads

```math
g^{-1} = $(_tex(:pmatrix, "r^{$(_tex(:transp))} & -r^{$(_tex(:transp))}t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")).
```

This computation can be done in-place of `h`.
"""

@doc "$(_doc_inv_SEn)"
Base.inv(G::SpecialEuclideanGroup, g)

@doc "$(_doc_inv_SEn)"
function inv!(G::SpecialEuclideanGroup, h, g)
    rg = submanifold_component(G, g, :Rotation)
    tg = submanifold_component(G, g, :Translation)
    rh = submanifold_component(G, h, :Rotation)
    th = submanifold_component(G, h, :Translation)
    copyto!(rh, transpose(rg))
    copyto!(th, -rh * tg)
    return h
end
function inv!(G::SpecialEuclideanGroup, q, ::Identity{<:SpecialEuclideanOperation})
    return identity_element!(G, q)
end

function ManifoldsBase.isapprox(
    G::SpecialEuclideanGroup, g::AbstractMatrix, h::AbstractMatrix; kwargs...
)
    return isapprox(g, h; kwargs...)
end

function Random.rand!(
    rng::AbstractRNG,
    G::SpecialEuclideanGroup,
    pX::AbstractMatrix;
    vector_at=nothing,
    kwargs...,
)
    SOn, Tn = _SOn_and_Tn(G)
    if vector_at === nothing # for points -> pass to manifold
        rand!(
            rng,
            SOn,
            submanifold_component(G, pX, :Rotation);
            vector_at=vector_at,
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, pX, :Translation);
            vector_at=vector_at,
            kwargs...,
        )
    else # for tangent vectors -> subset the vector_at as well.
        rand!(
            rng,
            SOn,
            submanifold_component(G, pX, :Rotation);
            vector_at=submanifold_component(G, vector_at, :Rotation),
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, pX, :Translation);
            vector_at=submanifold_component(G, vector_at, :Translation),
            kwargs...,
        )
    end
    return pX
end

function _SOn_and_Tn(G::LeftSpecialEuclideanGroup)
    return map(LieGroup, G.manifold.manifolds, G.op.operations)
end
function _SOn_and_Tn(G::RightSpecialEuclideanGroup)
    A = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return A[2], A[1]
end

function Base.show(io::IO, G::SpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[2].size)[1]
    return print(io, "SpecialEuclideanGroup($(size))")
end
function Base.show(io::IO, G::RightSpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[2].size)[1]
    return print(io, "SpecialEuclideanGroup($(size); variant=:right)")
end
