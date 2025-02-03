# An overview of Lie groups

# Alphabetical list of Lie groups

| Group | Manifold | ``∘`` | Comment |
|:------|:---------|:---------:|:------|
| [`GeneralLinearGroup`](@ref) | [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`HeisenbergGroup`](@ref) | [`HeisenbergMatrices`](@extref `Manifolds.HeisenbergMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`OrthogonalGroup`](@ref) | [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | This can be interpreted as all rotations and reflections. |
| [`SpecialEuclideanGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`)[`⋉`](@ref LeftSemidirectProductGroupOperation)[`Euclidean`](@extref `Manifolds.Rotations`) | [`∘`](@ref LeftSemidirectProductGroupOperation) | Analogously you can also use a [`⋊`](@ref RightSemidirectProductGroupOperation) if you prefer tuples `(t,R)` having the rotation matrix in the second component |
| [`SpecialOrthogonalGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SpecialUnitaryGroup`](@ref) | [`GeneralUnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`TranslationGroup`](@ref) | [`Euclidean`](@extref `Manifolds.Euclidean`) | [`+`](@ref AdditionGroupOperation) | |
| [`UnitaryGroup`](@ref) | [`UnitaryMatrices`](@extref `Manifolds.UnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
