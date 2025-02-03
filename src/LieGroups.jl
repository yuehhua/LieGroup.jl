@doc raw"""
LieGroups.jl: Lie groups and Lie algebras in Julia.

The package is named after the Norwegian mathematician [Marius Sophus Lie](https://en.wikipedia.org/wiki/Sophus_Lie) (1842–1899).

* 📚 Documentation: [juliamanifolds.github.io/LieGroups.jl/](https://juliamanifolds.github.io/LieGroups.jl/)
* 📦 Repository: [github.com/JuliaManifolds/LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl)
* 💬 Discussions: [github.com/JuliaManifolds/LieGroups.jl/discussions](https://github.com/JuliaManifolds/LieGroups.jl/discussions)
* 🎯 Issues: [github.com/JuliaManifolds/LieGroups.jl/issues](https://github.com/JuliaManifolds/LieGroups.jl/issues)
"""
module LieGroups

using LinearAlgebra, ManifoldsBase, Manifolds, StaticArrays, Random

import LinearAlgebra: adjoint, adjoint!

using ManifoldsBase: RealNumbers, ComplexNumbers, ℝ, ℂ

#
#
# = Compatibility (and a bit of type piracy for now)
# The following imports are necessary to use Manifolds.jl 0.10 with Lie groups
# The line is removed when the Groups are removed from possibly 0.11
import Manifolds:
    apply, apply!, compose, identity_element, is_identity, hat, hat!, vee, vee!
# Both define the following structs, so these for now lead to asking for explicit prefixes
# Manifolds: Identity, TranslationGroup
include("documentation_glossary.jl")
include("utils.jl")
include("interface.jl")
include("Lie_algebra/Lie_algebra_interface.jl")
# Generic Operations
include("group_operations/addition_operation.jl")
include("group_operations/multiplication_operation.jl")

# Actions
include("group_actions/group_action_interface.jl")
include("group_actions/group_operation_action.jl")

# Meta Lie groups
include("groups/power_group.jl")
include("groups/product_group.jl")
include("groups/semidirect_product_group.jl")

# Lie groups

include("groups/translation_group.jl")
include("groups/general_linear_group.jl")
include("groups/heisenberg_group.jl")

# includes generic implementations for O(n), U(n), SO(n), SO(n), so we load this first
include("groups/unitary_group.jl")
include("groups/orthogonal_group.jl")
include("groups/special_unitary_group.jl")
include("groups/special_orthogonal_group.jl")

# Products of Groups
include("groups/special_euclidean_group.jl")

export LieGroup, LieAlgebra
export PowerLieGroup, ProductLieGroup
export LeftSemidirectProductLieGroup, RightSemidirectProductLieGroup
export DefaultLieAlgebraOrthogonalBasis
export ×, ^, ⋉, ⋊
#
#
# Group Operations
export AbstractGroupOperation, Identity
export AdditionGroupOperation
export AbstractMultiplicationGroupOperation
export MatrixMultiplicationGroupOperation
export ProductGroupOperation
export LeftSemidirectProductGroupOperation, RightSemidirectProductGroupOperation

#
#
# Group Actions
export AbstractGroupActionType
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export LeftGroupOperationAction, RightGroupOperationAction
export GroupAction, GroupOperationAction
export InverseLeftGroupOperationAction, InverseRightGroupOperationAction

#
#
# Specific groups
export GeneralLinearGroup
export HeisenbergGroup
export OrthogonalGroup
export SpecialEuclideanGroup, SpecialOrthogonalGroup, SpecialUnitaryGroup
export TranslationGroup
export UnitaryGroup

# Points and Tangent representations
export AbstractLieGroupPoint, AbstractLieAlgebraTangentVector
export SpecialEuclideanMatrixPoint, SpecialEuclideanMatrixTangentVector
export SpecialEuclideanProductPoint, SpecialEuclideanProductTangentVector
# Functions
export adjoint, adjoint!, apply, apply!
export base_lie_group, base_manifold
export compose, compose!
export default_left_action,
    default_right_action,
    diff_apply,
    diff_apply!,
    diff_group_apply,
    diff_group_apply!,
    diff_left_compose,
    diff_left_compose!,
    diff_right_compose,
    diff_right_compose!
export get_coordinates, get_coordinates!, get_vector, get_vector!
export hat, hat!
export inv, inv!, inv_left_compose, inv_left_compose!, inv_right_compose, inv_right_compose!
export isapprox, is_point, is_vector
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export exp, exp!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export jacobian_conjugate, jacobian_conjugate!
export lie_bracket, lie_bracket!, log, log!
export manifold_dimension
export norm
export injectivity_radius
export rand, rand!
export switch
export vee, vee!
export zero_vector, zero_vector!
end # module LieGroups
