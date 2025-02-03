module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
# Implement SE(n) also on an Array Partition

Base.convert(::Type{<:ArrayPartition}, g::SpecialEuclideanProductPoint) = g.value
function Base.convert(::Type{SpecialEuclideanProductPoint}, p::ArrayPartition)
    return SpecialEuclideanProductPoint(p)
end
Base.convert(::Type{<:ArrayPartition}, X::SpecialEuclideanProductTangentVector) = X.value
function Base.convert(::Type{SpecialEuclideanProductTangentVector}, X::ArrayPartition)
    return SpecialEuclideanProductTangentVector(X)
end
# convert between both representations
function Base.convert(
    ::Type{AbstractMatrix}, g::SpecialEuclideanProductPoint{<:ArrayPartition}
)
    n = size(g.value.x[1])[1]
    A = zeros(n + 1, n + 1)
    A[n + 1, n + 1] = 1.0
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(g.value.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= g.value.x[1]
        A[1:n, n + 1] .= g.value.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= g.value.x[2]
        A[1:n, n + 1] .= g.value.x[1]
    end
end
function Base.convert(
    ::Type{AbstractMatrix}, X::SpecialEuclideanProductTangentVector{<:ArrayPartition}
)
    n = size(X.value.x[1])[1]
    A = zeros(n + 1, n + 1)
    A[n + 1, n + 1] = 0.0
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(X.value.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= X.value.x[1]
        A[1:n, n + 1] .= X.value.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= X.value.x[2]
        A[1:n, n + 1] .= X.value.x[1]
    end
end
# The reverse is hence not unique, but we can still cast to the special matrix point
function Base.convert(::Type{SpecialEuclideanMatrixPoint}, g::ArrayPartition)
    n = size(g.x[1])[1]
    A = zeros(n + 1, n + 1)
    A[n + 1, n + 1] = 1.0
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(g.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= g.x[1]
        A[1:n, n + 1] .= g.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= g.x[2]
        A[1:n, n + 1] .= g.x[1]
    end
    return SpecialEuclideanMatrixPoint(A)
end
function Base.convert(::Type{SpecialEuclideanMatrixTangentVector}, g::ArrayPartition)
    n = size(g.x[1])[1]
    A = zeros(n + 1, n + 1)
    A[n + 1, n + 1] = 0.0
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(g.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= g.x[1]
        A[1:n, n + 1] .= g.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= g.x[2]
        A[1:n, n + 1] .= g.x[1]
    end
    return SpecialEuclideanMatrixPoint(A)
end
# convert between both representation explicitly
# the inverse is again not unique since both (R,t) and (t,R) are possible results.
function Base.convert(
    ::Type{SpecialEuclideanMatrixPoint}, g::SpecialEuclideanProductPoint{<:ArrayPartition}
)
    return SpecialEuclideanMatrixPoint(convert(AbstractMatrix, g))
end
function Base.convert(
    ::Type{SpecialEuclideanMatrixTangentVector},
    g::SpecialEuclideanProductTangentVector{<:ArrayPartition},
)
    return SpecialEuclideanMatrixTangentVector(convert(AbstractMatrix, g))
end

function LieGroups.identity_element(
    G::LieGroups.LeftSpecialEuclideanGroup, ::Type{<:ArrayPartition}
)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(SOn), identity_element(Tn))
end
function LieGroups.identity_element(
    G::LieGroups.RightSpecialEuclideanGroup, ::Type{<:ArrayPartition}
)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(Tn), identity_element(SOn))
end
# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; v=1) = nothing

function ManifoldsBase.submanifold_component(
    G::SE,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{I},
) where {I<:Int,SE<:SpecialEuclideanGroup}
    # pass down to manifold by default
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), I
    )
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Rotation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 1
    )
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 2
    )
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 1
    )
end

function ManifoldsBase.zero_vector(
    G::LieAlgebra{
        𝔽,LieGroups.SpecialEuclideanGroupOperation,LieGroups.SpecialEuclideanGroup
    },
    ::ArrayPartition,
) where {𝔽}
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(n), zeros(n, n))
end
function ManifoldsBase.zero_vector(
    G::LieAlgebra{
        𝔽,LieGroups.SpecialEuclideanGroupOperation,LieGroups.SpecialEuclideanGroup
    },
    ::SpecialEuclideanProductTangentVector,
) where {𝔽}
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return SpecialEuclideanProductTangentVector(ArrayPartition(zeros(n), zeros(n, n)))
end
end
