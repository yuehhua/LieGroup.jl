"""
    AbstractMultiplicationGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `*` being overloaded, for example
`_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, a, b) = a * b`
"""
abstract type AbstractMultiplicationGroupOperation <: AbstractGroupOperation end

"""
    AbstractMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a matrix multiplication.
"""
struct MatrixMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

Base.:*(::Identity{MatrixMultiplicationGroupOperation}, p::Union{AbstractMatrix,Number}) = p
function Base.:*(
    p::Union{AbstractMatrix,Number}, ::Identity{MatrixMultiplicationGroupOperation}
)
    return p
end
function Base.:*(
    e::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return e
end
function Base.:*(
    e::Identity{<:AbstractMultiplicationGroupOperation}, ::Identity{AdditionGroupOperation}
)
    return e
end
function Base.:*(
    ::Identity{AdditionGroupOperation}, e::Identity{<:AbstractMultiplicationGroupOperation}
)
    return e
end

Base.:/(p, ::Identity{<:AbstractMultiplicationGroupOperation}) = p
Base.:/(::Identity{<:AbstractMultiplicationGroupOperation}, p) = inv(p)
function Base.:/(
    e::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return e
end

Base.:\(p, ::Identity{<:AbstractMultiplicationGroupOperation}) = inv(p)
Base.:\(::Identity{<:AbstractMultiplicationGroupOperation}, p) = p
function Base.:\(
    e::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return e
end

_doc_compose_mult = """
    compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h)

Compute the group operation composition of `g` and `h` with respect to
an [`AbstractMultiplicationGroupOperation`](@ref) on `G`, which falls back to calling
`g*h`, where `*` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_mult)"
compose(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_mult)"
compose!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}

function _compose!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}
    # perform the multiplication “safe”, that is, even when providing
    # one of the inputs `g,h`` and as output `k`
    (k === g || k === h) ? copyto!(k, g * h) : mul!(k, g, h)
    return k
end

Base.inv(e::Identity{<:AbstractMultiplicationGroupOperation}) = e

LinearAlgebra.det(::Identity{<:AbstractMultiplicationGroupOperation}) = true

_doc_diff_conjugate_add = """
    diff_conjugate(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the conjutage ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = ghg^{-1}``,
which simplifies for an [`AbstractMultiplicationGroupOperation`](@ref) to ``D(c_g(h))[X] = gXg^{-1}``.
"""

@doc "$(_doc_diff_conjugate_add)"
diff_conjugate(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_conjugate_add)"
function diff_conjugate!(
    G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    inv_right_compose!(G, Y, X, g) # Y = Xg^{-1}
    compose!(G, Y, g, Y) # Y = gY
    return Y
end

_doc_diff_inv_mult = """
    diff_inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, X)
    diff_inv!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, X)

Compute the value of differential ``Dι_{$(_math(:G))}(g)[X]`` of matrix inversion ``ι_{$(_math(:G))}(g) := g^{-1}`` at ``X ∈ 𝔤``
in the [`LieAlgebra`](@ref) ``𝔤`` of the [`LieGroup`](@ref) `G`.

The formula is given by

```math
Dι_{$(_math(:G))}(g)[X] = -g^{$(_tex(:transp))}Xg^{-1},
```

which stems from using the differential of the inverse from [Giles:2008](@cite) given by
``D(g^{-1})[X] = -g^{-1}Xg^{-1}`` composed with the push forward of the left composition
``Dλ_$(_math(:e))(g)[X] = gX`` mapping from the Liea algebra into the tangent space at ``g``,
and its adjoint ``D^*λ_$(_math(:e))(g)[X] = g^{$(_tex(:transp))}X``.
Then we get ``g^{$(_tex(:transp))}(g^{-1}(gX)g^{-1})`` which simplifies to ``-g^{$(_tex(:transp))}Xg^{-1}`` from above.
"""

@doc "$(_doc_diff_inv_mult)"
diff_inv(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, X) where {𝔽}

function diff_inv(
    ::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation},
    p::AbstractArray{<:Number,0},
    X::AbstractArray{<:Number,0},
) where {𝔽}
    p_inv = inv(p[])
    return -(p[] * X * p_inv)
end

@doc "$(_doc_diff_inv_mult)"
function diff_inv!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, p, X) where {𝔽}
    p_inv = inv(p)
    Z = X * p_inv
    mul!(Y, p', Z)
    Y .*= -1
    return Y
end

_doc_diff_left_compose_mult = """
    diff_left_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_left_compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h``,
which simplifies for an [`AbstractMultiplicationGroupOperation`](@ref) to ``Dλ_g(h)[X] = gX``.
"""

@doc "$(_doc_diff_left_compose_mult)"
diff_left_compose(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_left_compose_mult)"
function diff_left_compose!(
    G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, g * X)
end

_doc_diff_right_compose_mult = """
    diff_right_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g, X)
    diff_right_compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, h, g, X)

Compute the differential of the right group multiplication ``ρ_g(h) = h$(_math(:∘))g``,
which simplifies for an [`AbstractMultiplicationGroupOperation`](@ref) to ``Dρ_g(h)[X] = Xg``.
"""

@doc "$(_doc_diff_right_compose_mult)"
diff_right_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g, X) where {𝔽}

@doc "$(_doc_diff_right_compose_mult)"
function diff_right_compose!(
    G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X * g)
end

_doc_exp_mult = """
    exp(G::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, e::Identity{MatrixMultiplicationGroupOperation}, X, t::Number=1)
    exp!(G::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, g, e::Identity{MatrixMultiplicationGroupOperation}, X, t::Number=1)

Compute the Lie group exponential on a [`LieGroup`](@ref) with a [`MatrixMultiplicationGroupOperation`](@ref),
which simplifies to the [matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential).

This can be computed in-place of `g`.
"""

@doc "$(_doc_exp_mult)"
Base.exp(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    t::Number=1,
) where {𝔽} = exp(t * X)

@doc "$(_doc_exp_mult)"
function ManifoldsBase.exp!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    g,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    t::Number=1,
) where {𝔽}
    copyto!(g, exp(t .* X))
    return g
end

_doc_identity_element_mult = """
    identity_element(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`AbstractMultiplicationGroupOperation`](@ref) is the one-element or identity array.
"""

@doc "$(_doc_identity_element_mult)"
identity_element(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}) where {𝔽}

@doc "$(_doc_identity_element_mult)"
identity_element!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(
    ::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e::AbstractMatrix
) where {𝔽}
    return copyto!(e, LinearAlgebra.I)
end

_doc_inv_mult = """
    inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperationroupOperation}, g)
    inv!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for an [`AbstractMultiplicationGroupOperation`](@ref)
simplifies to the multiplicative inverse ``g^{-1}``. This can be done in-place of `h`.
"""

@doc "$(_doc_inv_mult)"
Base.inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g) where {𝔽}

@doc "$(_doc_inv_mult)"
function inv!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g) where {𝔽}
    copyto!(h, inv(g))
    return h
end
function inv!(
    G::LieGroup{𝔽,O}, q, ::Identity{O}
) where {𝔽,O<:AbstractMultiplicationGroupOperation}
    return identity_element!(G, q)
end

# Compute g^{-1}h more efficient than inverting g
function inv_left_compose!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, g \ h)
end
# Compute g∘h^{-1} more efficient than inverting h
function inv_right_compose!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, g / h)
end

_doc_lie_bracket_mult = """
    lie_bracket(::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, X, Y)
    lie_bracket!(::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`MatrixMultiplicationGroupOperation`](@ref) yields the
commutator bracket

```math
[X, Y] = XY-YX
```

The computation can be done in-place of `Z`.
"""

@doc "$(_doc_lie_bracket_mult)"
lie_bracket(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, X, Y) where {𝔽}

@doc "$(_doc_lie_bracket_mult)"
function lie_bracket!(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, Z, X, Y) where {𝔽}
    mul!(Z, X, Y)
    mul!(Z, Y, X, -1, true)
    return Z
end

_doc_log_mult = """
    log(G::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, e::Identity{MatrixMultiplicationGroupOperation}, g)
    log!(G::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, X, e::Identity{MatrixMultiplicationGroupOperation}, g)

Compute the Lie group logarithm on a [`LieGroup`](@ref) with a [`MatrixMultiplicationGroupOperation`](@ref),
which simplifies to the [matrix logarithm](https://en.wikipedia.org/wiki/Logarithm_of_a_matrix).

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_mult)"
Base.log(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
) where {𝔽} = log(g)
function Base.log(
    G::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    e::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
) where {𝔽}
    return zero_vector(G, e)
end

@doc "$(_doc_log_mult)"
function ManifoldsBase.log!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
) where {𝔽}
    copyto!(X, log(g))
    return X
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,MatrixMultiplicationGroupOperation},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
) where {𝔽}
    zero_vector!(LieAlgebra(G), X)
    return X
end

LinearAlgebra.mul!(q, ::Identity{<:AbstractMultiplicationGroupOperation}, p) = copyto!(q, p)
function LinearAlgebra.mul!(
    q::AbstractMatrix, p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}
)
    return copyto!(q, p)
end
function LinearAlgebra.mul!(
    q::Union{AbstractMatrix},
    ::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return copyto!(q, I)
end
function LinearAlgebra.mul!(
    q::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return q
end
Base.one(e::Identity{<:AbstractMultiplicationGroupOperation}) = e
