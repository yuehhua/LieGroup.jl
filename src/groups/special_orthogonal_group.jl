"""
    SpecialOrthogonalGroup{T}

The special orthogonal group ``$(_math(:SO))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`Rotations`](@extref `Manifolds.Rotations`).

# Constructor
    SpecialOrthogonalGroup(n; kwargs...)

Generate  special orthogonal group ``$(_math(:SO))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialOrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ,MatrixMultiplicationGroupOperation,Manifolds.Rotations{T}
}

function SpecialOrthogonalGroup(n; kwargs...)
    R = Manifolds.Rotations(n; kwargs...)
    return SpecialOrthogonalGroup{typeof(R).parameters[1]}(
        R, MatrixMultiplicationGroupOperation()
    )
end

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X)
@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g, X)

@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X)
@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g, X)

@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X)
@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g, X)

@doc "$(_doc_get_coordinates_On)"
ManifoldsBase.get_coordinates(
    G::SpecialOrthogonalGroup,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_coordinates_On)"
ManifoldsBase.get_coordinates!(
    G::SpecialOrthogonalGroup,
    c,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_vector_On)"
ManifoldsBase.get_vector(
    G::SpecialOrthogonalGroup,
    ::Identity{MatrixMultiplicationGroupOperation},
    c,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_vector_On)"
ManifoldsBase.get_vector!(
    G::SpecialOrthogonalGroup,
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    c,
    ::DefaultLieAlgebraOrthogonalBasis,
)

inv!(G::SpecialOrthogonalGroup, k, g) = copyto!(G, k, transpose(g))
function inv!(
    G::SpecialOrthogonalGroup, q, ::Identity{O}
) where {O<:AbstractMultiplicationGroupOperation}
    return identity_element!(G, q)
end

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g)

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X, g)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, g)

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g)

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X, g)

function Base.show(io::IO, G::SpecialOrthogonalGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "SpecialOrthogonalGroup($(size))")
end
