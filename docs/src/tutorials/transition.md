# Transition from `GroupManifolds` in `Manifolds.jl`

One predecessor of `LieGroups.jl` are the [`GroupManifold`](@extref `Manifolds.GroupManifold`)s in `Manifolds.jl`.
While this package provides the same features, one reason for a new package is,
that a “restart” offers the opportunity to put the main focus for the functions in this package
really on Lie groups.

This tutorial provides an overview of the necessary changes to your code if you based it on the predecessor.

## Table of function names and its successors

The following table lists all functions related to `GroupManifolds` and their new names
or replacements here in `LieGroups.jl`. In this code `G` always refers to the `GroupManifold`
in the first column and the [`LieGroup`](@ref) in the second.
Lie group elements (points) are always `g,h`,
Lie algebra elements (vectors) always `X, Y`.

New functions and types in this package are only mentioned, if they are worth a comment and if something changed.

The list is alphabetical, but first lists types, then functions

| `Manifolds.jl` | `LieGroups.jl` | Comment |
|:---------- |:---------- |:-------------- |
| `AdditionOperation` | [`AdditionGroupOperation`](@ref) | |
| `LeftForwardAction` | [`LeftGroupOperation`](@ref)
| `RightBackwardAction` | [`RightGroupOperation`](@ref) | |
| `LeftBackwardAction` | [`InverseRightGroupOperation`](@ref) | note that this is now also aa [`AbstractLeftGroupActionType`](@ref) |
| | [`LieAlgebra`](@ref)`(G)` | new alias to emphasize its manifold- and vector structure as well as for a few dispatch methods. |
| `GroupManifold(M, op)` | [`LieGroup`](@ref)`(M, op)` | |
| `RightForwardAction` | [`InverseLeftGroupOperation`](@ref) | note that this is an [`AbstractRightGroupActionType`](@ref) |
| `adjoint` | [`adjoint`](@ref) | now implemented with a default, when you provide [`diff_conjugate!`](@ref).
| `apply_diff` | [`diff_apply`](@ref) | modifiers (diff) come first, consistent with [`ManifoldsDiff.jl`](https://juliamanifolds.github.io/ManifoldDiff.jl/stable/) |
| `apply_diff_group` | [`diff_group_apply`](@ref) | modifiers (diff/group) come first, consistent with [`ManifoldsDiff.jl`](https://juliamanifolds.github.io/ManifoldDiff.jl/stable/) |
| | [`conjugate`](@ref), [`diff_conjugate`](@ref) | a new function to model ``c_g: \mathcal G → \mathcal G`` given by ``c_g(h) = g∘h∘g^{-1}`` |
| `exp(G, g, X)` | `exp(`[`base_manifold`](@ref base_manifold(G::LieGroup))`(G), g, X)` | the previous defaults whenever not agreeing with the invariant one can now be accessed on the internal manifold |
| `exp_inv(G, g, X)` | [`exp`](@ref exp(G::LieGroup, g, X, t::Number))`(G, g, X)`  | the exponential map invariant to the group operation is the default on Lie groups here |
| `exp_lie(G, X)` | [`exp`](@ref exp(G::LieGroup, e::Identity, X, t::Number))`(G, `[`Identity`](@ref)`(G), X)` | the (matrix) exponential is now the one at the [`Identity`](@ref)`(G)`, since there it agrees with the invariant one |
| `inverse_translate(G, g, h, c)` | [`inv_left_compose`](@ref)`(G, g, h)`, [`inv_right_compose`](@ref)`(G, g, h)` | compute ``g^{-1}∘h`` and ``g∘h^{-1}``, resp. |
| `inverse_tranlsate_diff(G, g, h, X, LeftForwardAction())` | - | discontinued, use `diff_left_compose(G, inv(G,g), h)` |
| `inverse_tranlsate_diff(G, g, h, X, RightBackwardAction())` | - | discontinued, use `diff_left_compose(G, h, inv(G,g))` |
| `log(G, g, h)` | `log(`[`base_manifold`](@ref base_manifold(G::LieGroup))`(G), g, h)` | you can now access the previous defaults on the internal manifold whenever they do not agree with the invariant one |
| `log_inv(G, g, h)` | [`log`](@ref log(G::LieGroup, g, h))`(G, g, h)` | the logarithmic map invariant to the group operation is the default on Lie groups here |
| `log_lie(G, g)` | [`log`](@ref log(G::LieGroup, e::Identity, g))`(G, `[`Identity`](@ref)`(G), g)` | the (matrix) logarithm is now the one at the identity, since there it agrees with the invariant one |
| `switch_direction(A)` | [`inv`](@ref inv(::AbstractGroupAction))`(A)` | switches from an action to its inverse action (formerly the direction forward/backward, sometimes even left/right, do not confuse with the side left/right). |
| `switch_side(A)` | [`switch`](@ref switch(::AbstractGroupAction))`(A)` | switches from a left action to its corresponding right action. |
| `translate(G, g, h)` | [`compose`](@ref)`(G, g, h)` | unified to `compose` |
| `translate_diff(G, g, X, c)` | [`diff_left_compose`](@ref)`(G, g, h, X)`, [`diff_right_compose`](@ref)`(G, g, h, X)` | for compose ``∘(g, h)`` the functions now specify whether the derivative is taken w.r.t. to the left (`g`) or right (`h`) argument |
|`VeeOrthogonalBasis` | [`LieAlgebraOrthogonalBasis`](@ref) | |

# Notable changes

* The [`GeneralLinearGroup`](@ref) (formerly `GeneralLinear`) switched to using its Lie algebra to represent tangent vectors.