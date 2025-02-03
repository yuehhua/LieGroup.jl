#
#
# Together with a product group operation, that is a tuple of operations,
# One for each factor, we define the product operation as acting element wise.

"""
    AbstractProductGroupOperation <: AbstractGroupOperation

An abstract type to model group operations on a product manifold
"""
abstract type AbstractProductGroupOperation <: AbstractGroupOperation end

"""
    ProductGroupOperation{O<:<:NTuple{N,AbstractGroupOperation} where N} <: AbstractProductGroupOperation

A struct do model a tuple of group operations, one for each factor of a product group,
that together forms a new group operation.

Access to the single operations can be done by `pgo[i]`.

# Constructor

    ProductGroupOperation(o::AbstractGroupOperation...)
    ×(o::AbstractGroupOperation...) = ProductGroupOperation(o...)
"""
struct ProductGroupOperation{OTM<:NTuple{N,AbstractGroupOperation} where {N}} <:
       AbstractProductGroupOperation
    operations::OTM
end
function ProductGroupOperation(operations::AbstractGroupOperation...)
    return ProductGroupOperation(operations)
end

@inline Base.getindex(pgo::ProductGroupOperation, i::Integer) = pgo.ooperations[i]

@doc raw"""
    cross(O1::AbstractGroupOperation, O2::AbstractGroupOperation)
    O1 × O2
    O1 × O2 × O3 × ...

Return the [`ProductGroupOperation`](@ref) For two [AbstractGroupOperation`](@ref) `O1` and `O2`,
where for the case that one of them is a [`ProductGroupOperation`](@ref) itself,
the other is either prepended (if `O1` is a product) or appended (if `O2` is).
If both are product operations, they are combined into one, keeping the order of operations.

For the case that more than two are concatenated with `×` this is iterated.
"""
cross(::AbstractGroupOperation...)
function LinearAlgebra.cross(O1::AbstractGroupOperation, O2::AbstractGroupOperation)
    return ProductGroupOperation(O1, O2)
end
function LinearAlgebra.cross(P::ProductGroupOperation, O::AbstractGroupOperation)
    return ProductGroupOperation(P.operations..., O)
end
function LinearAlgebra.cross(O::AbstractGroupOperation, P::ProductGroupOperation)
    return ProductGroupOperation(O, P.operations...)
end
function LinearAlgebra.cross(P1::ProductGroupOperation, P2::ProductGroupOperation)
    return ProductGroupOperation(P1.operations..., P2.operations...)
end

"""
    ProductLieGroup(G, H, ...)

Return the [`LieGroup`](@ref) of the product of Lie groups `G` and `H`.

Alternatively, the short hand `G × H` can be used.
"""
function ProductLieGroup(G::LieGroup, H::LieGroup)
    return LieGroup(G.manifold × H.manifold, G.op × H.op)
end

function ManifoldsBase.submanifold_components(
    ::LieGroup{𝔽,Op,M}, op::ProductGroupOperation
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    return op.operations
end

function _compose!(
    PrG::LieGroup{𝔽,Op,M}, k, g, h
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    map(
        compose!,
        map(LieGroup, PrG.manifold.manifolds, PrG.op.operations),
        submanifold_components(PrG.manifold, k),
        submanifold_components(PrG.manifold, g),
        submanifold_components(PrG.manifold, h),
    )
    return k
end

function ManifoldsBase.check_size(
    PrG::LieGroup{𝔽,Op,M}, g
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, g)
end
function ManifoldsBase.check_size(
    ::LieGroup{𝔽,Op,M}, ::Identity
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    return nothing
end
function ManifoldsBase.check_size(
    PrG::LieGroup{𝔽,Op,M}, g, X
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, g, X)
end
function ManifoldsBase.check_size(
    PrG::LieGroup{𝔽,Op,M}, ::Identity, X
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, identity_element(PrG, typeof(X)), X)
end

function conjugate!(
    PrG::LieGroup{𝔽,Op,M}, k, g, h
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        conjugate!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, k),
        submanifold_components(PrM, g),
        submanifold_components(PrM, h),
    )
    return k
end

@doc raw"""
    cross(G::LieGroup, H::LieGroup)
    G × H
    G1 × G2 × G3 × ...

Return the [`ProductLieGroup`](@ref) For two [`LieGroups`](@ref) `G` and `H`,
where for the case that one of them is a [`ProductLieGroup`](@ref) itself,
the other is either prepended (if `H` is a product) or appended (if `G` is).
If both are product Lie groups, they are combined into one, keeping the order of operations.

For the case that more than two are concatenated with `×` this is iterated.
"""
LinearAlgebra.cross(::LieGroup...)
function LinearAlgebra.cross(G::LieGroup, H::LieGroup)
    return ProductLieGroup(G, H)
end

function diff_conjugate!(
    PrG::LieGroup{𝔽,Op,M}, Y, g, h, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        diff_conjugate!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_inv!(
    PrG::LieGroup{𝔽,Op,M}, Y, g, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        diff_inv!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_left_compose!(
    PrG::LieGroup{𝔽,Op,M}, Y, g, h, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        diff_left_compose!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_right_compose!(
    PrG::LieGroup{𝔽,Op,M}, Y, g, h, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        diff_right_compose!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function ManifoldsBase.exp!(
    PrG::LieGroup{𝔽,Op,M}, h, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        (M, h, X) -> exp!(M, h, X), # introduce a function with “hard coded” t
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, X),
    )
    return h
end

function ManifoldsBase.exp!(
    PrG::LieGroup{𝔽,Op,M}, h, g, X
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        (M, h, g, X) -> exp!(M, h, g, X), # introduce a function with “hard coded” t
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, g),
        submanifold_components(PrM, X),
    )
    return h
end

function get_vector_lie!(
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}}, X, c, N
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, submanifold_components(PrM, X), Prc)
    map(ts) do t
        return get_vector_lie!(t..., N)
    end
    return X
end

function identity_element!(
    PrG::LieGroup{𝔽,Op,M}, e
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        identity_element!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, e),
    )
    return e
end

function inv!(
    PrG::LieGroup{𝔽,Op,M}, h, g
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, g),
    )
    return h
end
function inv!(
    PrG::LieGroup{𝔽,Op,M}, h, ::Identity{Op}
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        map(Identity, PrG.op.operations),
    )
    return h
end

function lie_bracket!(
    PrA::LieAlgebra{𝔽,Op,<:LieGroup{𝔽,Op,M}}, Z, X, Y
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrA.manifold.manifold
    map(
        lie_bracket!,
        map(LieAlgebra, LieGroup.(PrM.manifolds, PrA.manifold.op.operations)),
        submanifold_components(PrM, Z),
        submanifold_components(PrM, X),
        submanifold_components(PrM, Y),
    )
    return Z
end

function ManifoldsBase.log!(
    PrG::LieGroup{𝔽,Op,M}, X, ::Identity{Op}
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    zero_vector!(LieAlgebra(PrG), X)
    return X
end

function ManifoldsBase.log!(
    PrG::LieGroup{𝔽,Op,M}, X, g, h
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        log!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        submanifold_components(PrM, g),
        submanifold_components(PrM, h),
    )
    return X
end
function ManifoldsBase.log!(
    PrG::LieGroup{𝔽,Op,M}, X, e::Identity{Op}, h
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        log!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        submanifold_components(PrM, e),
        submanifold_components(PrM, h),
    )
    return X
end
function ManifoldsBase.log!(
    PrG::LieGroup{𝔽,Op,M}, X, h
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = PrG.manifold
    map(
        log!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        submanifold_components(PrM, h),
    )
    return X
end

function Base.show(
    io::IO, G::LieGroup{𝔽,Op,M}
) where {𝔽,Op<:ProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = G.manifold.manifolds
    ops = G.op.operations
    return print(io, "ProductLieGroup($(join(PrM, " × ")), $(join(ops, " × ")))")
end

function get_coordinates_lie!(
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}}, c, X, N
) where {𝔽,Op<:AbstractProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc, submanifold_components(PrM, X))
    map(ts) do t
        return get_coordinates_lie!(t..., N)
    end
    return c
end
