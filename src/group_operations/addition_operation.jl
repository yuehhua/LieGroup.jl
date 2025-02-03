"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`_compose(G::LieGroup{𝔽,AdditionGroupOperation}, a, b) = a + b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

#
#
# Handle interactions of `+` and `-` with the identity element, though they are
# also already handled on the `compose()` level
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(g, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(g, ::Identity{AdditionGroupOperation}) = g

_doc_compose_add = """
    compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h)

Compute the group operation composition of `g` and `h` with respect to
the [`AdditionGroupOperation`](@ref) on `G`, which falls back to calling
`g+h`, where `+` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_add)"
compose(::LieGroup{𝔽,AdditionGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_add)"
compose!(::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}

function _compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}
    k .= g .+ h
    return k
end

_doc_diff_conjugate_add = """
    diff_conjugate(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X)

Compute the differential of the conjutage ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = g+h-g = h``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``D(c_g(h))[X] = X``.
"""

@doc "$(_doc_diff_conjugate_add)"
diff_conjugate(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_conjugate_add)"
function diff_conjugate!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_diff_inv_add = """
    diff_inv(G::LieGroup{𝔽,AdditionGroupOperation}, g, X)
    diff_inv!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, X)

Compute the differential of the inverse operation ``ι_{$(_math(:G))}(g) = g^-1 = -g``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``Dι_{$(_math(:G))}(g)[X] = -X``
"""

@doc "$(_doc_diff_inv_add)"
diff_inv(G::LieGroup{𝔽,AdditionGroupOperation}, g, X) where {𝔽}

@doc "$(_doc_diff_inv_add)"
function diff_inv!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, X) where {𝔽}
    Y .= (-1) .* X
    return Y
end

_doc_diff_left_compose_add = """
    diff_left_compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X)
    diff_left_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``Dλ_g(h)[X] = X``.
"""

@doc "$(_doc_diff_left_compose_add)"
diff_left_compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_left_compose_add)"
function diff_left_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_diff_right_compose_add = """
    diff_right_compose(G::LieGroup{𝔽,AdditionGroupOperation}, h, g, X)
    diff_right_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, h, g, X)

Compute the differential of the right group multiplication ``ρ_g(h) = h$(_math(:∘))g``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``Dρ_g(h)[X] = X``.
"""

@doc "$(_doc_diff_right_compose_add)"
diff_right_compose(::LieGroup{𝔽,AdditionGroupOperation}, ::Any, ::Any, ::Any) where {𝔽}

@doc "$(_doc_diff_right_compose_add)"
function diff_right_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_exp_add = """
    exp(G::LieGroup{𝔽,AdditionGroupOperation}, X)
    exp!(G::LieGroup{𝔽,AdditionGroupOperation}, g, X)

Compute the Lie group exponential on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `g`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``g=0+X=X``.
"""

@doc "$(_doc_exp_add)"
ManifoldsBase.exp(::LieGroup{𝔽,AdditionGroupOperation}, X) where {𝔽} = X

@doc "$(_doc_exp_add)"
function ManifoldsBase.exp!(::LieGroup{𝔽,AdditionGroupOperation}, g, X) where {𝔽}
    g .= X
    return g
end

_doc_identity_element_add = """
    identity_element(G::LieGroup{𝔽,AdditionGroupOperation})
    identity_element!(G::LieGroup{𝔽,AdditionGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for the [`AdditionGroupOperation`](@ref) is the zero element or array.
"""

@doc "$(_doc_identity_element_add)"
identity_element(::LieGroup{𝔽,AdditionGroupOperation}) where {𝔽}

@doc "$(_doc_identity_element_add)"
function identity_element!(::LieGroup{𝔽,AdditionGroupOperation}, e) where {𝔽}
    return fill!(e, 0)
end

_doc_inv_add = """
    inv(G::LieGroup{𝔽,AdditionGroupOperation}, g)
    inv!(G::LieGroup{𝔽,AdditionGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for the [`AdditionGroupOperation`](@ref)
simplifies to ``-g``. This can be done in-place of `h`.
"""

@doc "$(_doc_inv_add)"
Base.inv(G::LieGroup{𝔽,AdditionGroupOperation}, g) where {𝔽}

@doc "$(_doc_inv_add)"
function inv!(::LieGroup{𝔽,AdditionGroupOperation}, h, g) where {𝔽}
    h .= (-1) .* g
    return h
end
# Resolve ambiguity
function inv!(
    G::LieGroup{𝔽,AdditionGroupOperation}, q, ::Identity{AdditionGroupOperation}
) where {𝔽}
    return identity_element!(G, q)
end

_doc_lie_bracket_add = """
    lie_bracket!(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, X, Y)
    lie_bracket!(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`AdditionGroupOperation`](@ref) simplifies to the
corresponding $(_link(:zero_vector)).
The computation can be done in-place of `Z`.
"""

@doc "$(_doc_lie_bracket_add)"
lie_bracket(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, X, Y) where {𝔽}

@doc "$(_doc_lie_bracket_add)"
function lie_bracket!(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, Z, X, Y) where {𝔽}
    return zero_vector!(𝔤, Z)
end

_doc_log_add = """
    log(G::LieGroup{𝔽,AdditionGroupOperation}, g)
    log!(G::LieGroup{𝔽,AdditionGroupOperation}, X, g)

Compute the Lie group logarithm on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `X`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``X=g-0=g``.
"""

@doc "$(_doc_log_add)"
ManifoldsBase.log(::LieGroup{𝔽,AdditionGroupOperation}, ::Any) where {𝔽}
@doc "$(_doc_log_add)"
function ManifoldsBase.log!(G::LieGroup{𝔽,AdditionGroupOperation}, X, g) where {𝔽}
    return copyto!(G, X, g)
end
function ManifoldsBase.log!(
    ::LieGroup{𝔽,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}
) where {𝔽}
    return fill!(X, 0)
end
