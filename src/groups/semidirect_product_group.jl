
"""
    SemiDirectProductGroupOperation{
        O1<:AbstractGroupOperation,
        O2<:AbstractGroupOperation,
        A<:AbstractGroupActionType
    } <: AbstractProductGroupOperation

An abstract type for all semdirect product group operations.
"""
abstract type SemiDirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: AbstractProductGroupOperation end

"""
    LeftSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a semidirect Lie group product.

Let ``($(_tex(:Cal, "N")), ⋄)`` and ``($(_tex(:Cal, "H")), ⋆)`` be two Lie groups
with group operations ``⋄`` and ``⋆``, respectively, as well as a group action
``σ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``, cf [`AbstractLeftGroupActionType`](@ref).

We use here as well use the notation ``σ_h: $(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``
as a family of maps on ``$(_tex(:Cal, "N"))``

Then we define a group operation ``∘`` on the product manifold ``$(_tex(:Cal, "N"))×$(_tex(:Cal, "H"))`` by

```math
    (h_1,n_1) ∘ (h_2,n_2) := (h_1 ⋆ h_2, σ_{h_2}(n_1) ⋄ n_2).
```

See [HilgertNeeb:2012; Definition 9.2.22](@cite), second definition for more details.

# Constructor

    LeftSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `op2::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "N"))``
* `action::AbstractGroupActionType` The group action ``σ`` of ``$(_tex(:Cal, "H"))`` on ``$(_tex(:Cal, "N"))``

"""
struct LeftSemidirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: SemiDirectProductGroupOperation{O1,O2,A}
    operations::Tuple{O1,O2}
    action_type::A
    function LeftSemidirectProductGroupOperation(
        op1::O1, op2::O2, action::A
    ) where {
        O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
    }
        return new{O1,O2,A}((op1, op2), action)
    end
end
@inline Base.getindex(spgo::SemiDirectProductGroupOperation, i::Integer) =
    spgo.operations[i]

"""
    RightSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a right semidirect Lie group product.

Let ``($(_tex(:Cal, "N")), ⋄)`` and ``($(_tex(:Cal, "H")), ⋆)`` be two Lie groups
with group operations ``⋄`` and ``⋆``, respectively, as well as a group action
``σ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``, cf [`AbstractGroupActionType`](#ref).

We use here as well use the notation ``σ_h: $(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``
as a family of maps on ``$(_tex(:Cal, "N"))``

Then we define a group operation ``∘`` on the product manifold ``$(_tex(:Cal, "N"))×$(_tex(:Cal, "H"))`` by

```math
    (n_1,h_1) ∘ (n_2,h_2) := (n_1 ⋄ σ_{h_1}(n_2), h_1 ⋆ h_2)
```

See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition for more details.

# Constructor

    RightSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "N"))``
* `op2::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `action::AbstractGroupActionType`: The group action ``σ`` of ``$(_tex(:Cal, "H"))`` on ``$(_tex(:Cal, "N"))``

"""
struct RightSemidirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: SemiDirectProductGroupOperation{O1,O2,A}
    operations::Tuple{O1,O2}
    action_type::A
    function RightSemidirectProductGroupOperation(
        op1::O1, op2::O2, action::A
    ) where {
        O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
    }
        return new{O1,O2,A}((op1, op2), action)
    end
end

"""
    LeftSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_left_action(N, H)
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) = N ⋉ H`` for an [`AbstractLeftGroupActionType`](@ref)
using the [`LeftSemidirectProductGroupOperation`](@ref) for the group operation definition
as well as [HilgertNeeb:2012; Definition 9.2.22](@cite), second definition, for more details.

The short form `N `[`⋉`](@ref ⋉(L1::LieGroup, L2::LieGroup))` H` can be used if the
corresponding [`default_left_action`](@ref)`(N,H)` is the one you want to use.
"""
function LeftSemidirectProductLieGroup(
    N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_left_action(N, H)
)
    return LieGroup(
        N.manifold × H.manifold, LeftSemidirectProductGroupOperation(N.op, H.op, action)
    )
end

"""
    RightSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_right_action(N,H)
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) = N ⋊ H`` for an [`AbstractLeftGroupActionType`](@ref)
using the [`RightSemidirectProductGroupOperation`](@ref) for the group operation definition
as well as [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition, for more details.

The short form `N `[`⋊`](@ref ⋊(L1::LieGroup, L2::LieGroup))` H` can be used if the
corresponding [`default_right_action`](@ref)`(N,H)` is the one you want to use.
"""
function RightSemidirectProductLieGroup(
    N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_right_action(N, H)
)
    return LieGroup(
        N.manifold × H.manifold, RightSemidirectProductGroupOperation(N.op, H.op, action)
    )
end

"""
    L1 ⋉ L2
    ⋉(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`LeftSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_left_action`](@ref)`(L1, L2)` is used.
"""
function ⋉(L1::LieGroup, L2::LieGroup)
    return LeftSemidirectProductLieGroup(L1, L2, default_left_action(L1, L2))
end

"""
    L1 ⋊ L2
    ⋊(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`RightSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_right_action`](@ref)`(L1, L2)` is used.
"""
function ⋊(L1::LieGroup, L2::LieGroup)
    return RightSemidirectProductLieGroup(L1, L2, default_right_action(L1, L2))
end

#
#
# Functions

"""
    compose(L::LieGroup{𝔽,LeftSemidirectProductGroupOperation}, g, h)

Compute the group operation ``$(_math(:∘))`` on the semidirect product Lie group ``L = G ⋉ H``,
that is for `g` `` = (g_1,h_1)``, `h` ``= (g_2,h_2)`` with ``g_1,g_2 ∈ G``, ``h_1,h_2 ∈ H``
this computes

```math
    (g_1,h_1) ∘ (g_2,h_2) := (g_1 ⋄ g_2, h_1 ⋆ σ_{g_1}(h_2)).
```
where ``∘`` denotes the group operation on ``L``, ``⋄`` and ``⋆`` those on ``G`` and ``H``,
respectively, and ``σ`` is the group action specified by the [`AbstractGroupActionType`](#ref) within the [`LeftSemidirectProductLieGroup`](@ref)  ``L``.
"""
compose(
    SDPG::LieGroup{𝔽,LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
    ::Any,
    ::Any,
) where {𝔽}

function _compose!(
    SDPG::LieGroup{𝔽,<:LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
    k,
    g,
    h,
) where {𝔽}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    # We have to perform 3 steps applying the group action
    # 1) x = σ_g[1](h[2]) (a point on H)
    # 2) compose g[1] and h[1] (a point on G, neither input needed afterwards)
    # 3) compose g[2] and x (a point on H, neither input needed afterwards)
    # to avoid to overwrite elements in case k=g or k=h: allocate for result of (1)
    # especially after (1) we still need g[2] (in case k=g)
    x = copy(H, submanifold_component(SDPG, k, 2))
    # (1)
    apply!(A, x, submanifold_component(SDPG, g, 1), submanifold_component(SDPG, h, 2))
    # (2)
    _compose!(
        G,
        submanifold_component(SDPG, k, 1),
        submanifold_component(SDPG, g, 1),
        submanifold_component(SDPG, h, 1),
    )
    _compose!(H, submanifold_component(SDPG, k, 2), submanifold_component(SDPG, g, 2), x)
    return k
end
"""
    compose(L::LieGroup{𝔽,RightSemidirectProductGroupOperation}, g, h)

Compute the group operation ``$(_math(:∘))`` on the semidirect product Lie group ``L = G ⋊ H``,
that is for `g` `` = (g_1,h_1)``, `h` ``= (g_2,h_2)`` with ``g_1,g_2 ∈ G``, ``h_1,h_2 ∈ H``
this computes

```math
    (g_1,h_1) ∘ (g_2,h_2) := (g_1 ⋄ σ_{h_1}(g_2), h_1 ⋆ h_2),
```

where ``∘`` denotes the group operation on ``L``, ``⋄`` and ``⋆`` those on ``G`` and ``H``,
respectively, and ``σ`` is the group action specified by the [`AbstractGroupActionType`](#ref) within the [`RightSemidirectProductLieGroup`](@ref) ``L``.
"""
compose(
    SDPG::LieGroup{𝔽,RightSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
    ::Any,
    ::Any,
) where {𝔽}
function _compose!(
    SDPG::LieGroup{
        𝔽,<:RightSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold
    },
    k,
    g,
    h,
) where {𝔽}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, H, G)
    # We have to perform 3 steps applying the group action
    # 1) x = σ_g[2](h[1]) (a point on G)
    # 2) compose g[1] and x (a point on G)
    # 3) compose g[2] and h[2] (a point on H)
    # to avoid to overwrite elements in case k=g or k=h: allocate for result of (1)
    # especially after (1) we still need g[1] (in case k=g)
    x = copy(G, submanifold_component(SDPG, k, 1))
    # (1)
    apply!(A, x, submanifold_component(SDPG, g, 2), submanifold_component(SDPG, h, 1))
    # (2)
    _compose!(G, submanifold_component(SDPG, k, 1), submanifold_component(SDPG, g, 1), x)
    # (3)
    _compose!(
        H,
        submanifold_component(SDPG, k, 2),
        submanifold_component(SDPG, g, 2),
        submanifold_component(SDPG, h, 2),
    )
    return k
end

function get_vector_lie!(
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}}, X, c, N
) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, submanifold_components(PrG, X), Prc)
    map(ts) do t
        return get_vector_lie!(t..., N)
    end
    return X
end

"""
    inv(SDPG::LieGroup{𝔽,Op,M}, g) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ProductManifold}

Compute the inverse element of an element ``g = (g_1,g_2)`` given by

```math
g^{-1} = (g_1^{-1}, σ_{g_1^{-1}}g_2^{-1}).
```

for the left variant and

```math
g^{-1} = (σ_{g_2^{-1}} g_1^{-1}, g_2^{-1})
```

for the right variant, respectively. See also [HilgertNeeb:2012; Proof of Lemma 2.2.3](@cite).
"""
Base.inv(
    SDPG::LieGroup{𝔽,Op,M}, g
) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}

function inv!(
    SDPG::LieGroup{𝔽,O,M}, k, g
) where {𝔽,O<:LeftSemidirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    inv!(G, submanifold_component(SDPG, k, 1), submanifold_component(PM, g, 1))
    inv!(H, submanifold_component(SDPG, k, 2), submanifold_component(PM, g, 2))
    apply!( # Apply the group action with g1^-1 to g2^-1
        A,
        submanifold_component(SDPG, k, 2),
        submanifold_component(SDPG, k, 1),
        submanifold_component(SDPG, k, 2),
    )
    return k
end
function inv!(
    SDPG::LieGroup{𝔽,O,M}, k, g
) where {𝔽,O<:RightSemidirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    inv!(G, submanifold_component(SDPG, k, 1), submanifold_component(PM, g, 1))
    inv!(H, submanifold_component(SDPG, k, 2), submanifold_component(PM, g, 2))
    apply!( # Apply the group action with g2^-1 to g1^-1
        A,
        submanifold_component(SDPG, k, 1),
        submanifold_component(SDPG, k, 2),
        submanifold_component(SDPG, k, 1),
    )
    return k
end
function inv!(
    SDPG::LieGroup{𝔽,O,M}, k, ::Identity{O}
) where {𝔽,O<:LeftSemidirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = SDPG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, SDPG.op.operations),
        submanifold_components(PrM, k),
        map(Identity, SDPG.op.operations),
    )
    return k
end
function inv!(
    SDPG::LieGroup{𝔽,O,M}, k, ::Identity{O}
) where {𝔽,O<:RightSemidirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrM = SDPG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, SDPG.op.operations),
        submanifold_components(PrM, k),
        map(Identity, SDPG.op.operations),
    )
    return k
end
function identity_element!(
    SDPG::LieGroup{𝔽,Op,M}, e
) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    GH = map(LieGroup, SDPG.manifold.manifolds, SDPG.op.operations)
    identity_element!.(GH, submanifold_components(SDPG.manifold, e))
    return e
end

function Base.show(
    io::IO,
    SDPG::LieGroup{𝔽,<:LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = SDPG.op.action_type
    return print(io, "LeftSemidirectProductLieGroup($G, $H, $at)")
end
function Base.show(
    io::IO,
    SDPG::LieGroup{
        𝔽,<:RightSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold
    },
) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = SDPG.op.action_type
    return print(io, "RightSemidirectProductLieGroup($G, $H, $at)")
end

function get_coordinates_lie!(
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}}, c, X, N
) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc, submanifold_components(PrG, X))
    map(ts) do t
        return get_coordinates_lie!(t..., N)
    end
    return c
end
