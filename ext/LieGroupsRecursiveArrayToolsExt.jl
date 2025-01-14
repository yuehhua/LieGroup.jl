module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
# Implement SE(n) also on an Array Partition

# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; v=1) = nothing

function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 1)
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Translation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Translation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 1)
end

# TODO: Implement?
# The following three should also work due to the
# AbstractProductGroup s implementations
# * `check_point`
# * `check_size`
# * `check_vector`
end
