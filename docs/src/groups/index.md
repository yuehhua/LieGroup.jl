# An overview of Lie groups

# Alphabetical list of Lie groups

| Group | Manifold | ``∘`` | Comment |
|:------|:---------|:---------:|:------|
| [`GeneralLinearGroup`](@ref) | [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`OrthogonalGroup`](@ref) | [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | This can be interpreted as all rotations and reflections. |
| [`SpecialEuclideanGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`)`×`[`Euclidean`](@extref `Manifolds.Rotations`) | [`∘`](@ref LeftSemidirectProductGroupOperation) | |
| [`SpecialOrthogonalGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SpecialUnitaryGroup`](@ref) | [`GeneralUnitary`](@extref `Manifolds.GeneralUnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`TranslationGroup`](@ref) | [`Euclidean`](@extref `Manifolds.Euclidean`) | [`+`](@ref AdditionGroupOperation) | |
| [`UnitaryGroup`](@ref) | [`UnitaryMatrices`](@extref `Manifolds.UnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
