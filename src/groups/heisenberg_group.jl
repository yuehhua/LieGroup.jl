
@doc """
    HeisenbergGroup{T}

The `HeisenbergGroup(n)` is the group of ``(n+2)×(n+2)`` matrices,
see also [BinzPods:2008](@cite) or [Heisenberg group](https://en.wikipedia.org/wiki/Heisenberg_group)
where `T` specifies the `eltype` of the matrix entries. `

```math
$(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
)),
```

where ``I_n`` is the ``n×n`` unit matrix, ``$(_tex(:vec, "a")), $(_tex(:vec, "b")) ∈ ℝ^n`` are vectors of length ``n``,
``$(_tex(:vec, 0))_n`` is the zero vector of length ``n``, and ``c ∈ ℝ`` is a real number.
The group operation is matrix multiplication.

The Lie Algebra consists of the elements
```math
$(_tex(:pmatrix,
    "0 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & Z_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 0"
)),
```
where additionally ``Z_n`` denotes the ``n×n`` zero matrix.
"""
const HeisenbergGroup{T} = LieGroup{
    ℝ,MatrixMultiplicationGroupOperation,Manifolds.HeisenbergMatrices{T}
}

function HeisenbergGroup(n::Int; parameter::Symbol=:type)
    Hm = Manifolds.HeisenbergMatrices(n; parameter=parameter)
    return HeisenbergGroup{typeof(Hm).parameters...}(
        Hm, MatrixMultiplicationGroupOperation()
    )
end

function _heisenberg_a_view(G::HeisenbergGroup, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return view(g, 1, 2:(n + 1))
end
function _heisenberg_b_view(G::HeisenbergGroup, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return view(g, 2:(n + 1), n + 2)
end

@doc """
    exp(G::HeisenbergGroup, X)
    exp!(G::HeisenbergGroup, g, X)

Compute the Lie group exponential for the [`HeisenbergGroup`](@ref) `G` of the vector `X`.

For ``X = $(_tex(:pmatrix,
    "0 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & Z_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 0"
))``
from the Lie algebra of the Heisenberg group,
where ``$(_tex(:vec, "a")), $(_tex(:vec, "b")) ∈ ℝ^n`` vectors of length ``n``,
``$(_tex(:vec, 0))_n`` is the zero vector of length ``n``, ``c ∈ ℝ``, and
``Z_n`` denotes the ``n×n`` zero matrix.

Then the

```math
$(_tex(:exp))_{$(_tex(:Cal, "G"))}(X) =
$(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c + $(_tex(:frac,"1","2"))$(_tex(:vec, "a"))^{$(_tex(:transp))}$(_tex(:vec, "b"))",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
)),
```
where ``I_n`` is the ``n×n`` unit matrix.

This can be computed in-place of the Lie group element `g`.
"""
ManifoldsBase.exp(G::HeisenbergGroup, X)

function ManifoldsBase.exp!(G::HeisenbergGroup, h, X)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    copyto!(h, I)
    a_view = _heisenberg_a_view(G, X)
    b_view = _heisenberg_b_view(G, X)
    h[1, 2:(n + 1)] .= a_view
    h[2:(n + 1), n + 2] .= b_view
    h[1, n + 2] = X[1, n + 2] + dot(a_view, b_view) / 2
    return h
end

@doc """
    exp(G::HeisenbergGroup, g, X)

Exponential map on the [`HeisenbergGroup`](@ref) `G` with the left-invariant metric.

We denote by `g` a point on the Heisenberg group and by ``X`` a vector from the Lie algebra.
These are of the form

```math
g = $(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
))
$(_tex(:qquad))
X = $(_tex(:pmatrix,
    "0 & $(_tex(:vec, "d"))^{$(_tex(:transp))} & f",
    "$(_tex(:vec, 0))_n & Z_n & $(_tex(:vec, "e"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 0"
)),
```
where ``I_n`` is the ``n×n`` unit matrix, ``Z_n`` is the ``n×n`` zero matrix,
``$(_tex(:vec, "a")), $(_tex(:vec, "b")), $(_tex(:vec, "d")), $(_tex(:vec, "e")) ∈ ℝ^n`` are vectors of length ``n``,
``$(_tex(:vec, 0))_n`` is the zero vector of length ``n``, and ``c,f ∈ ℝ`` are real numbers.

Then the formula reads
```math
$(_tex(:exp))_g(X) =
$(_tex(:pmatrix,
    "1 & ($(_tex(:vec, "a"))+$(_tex(:vec, "d")))^{$(_tex(:transp))} & c+f+$(_tex(:frac,"1","2"))$(_tex(:vec, "d"))^{$(_tex(:transp))}$(_tex(:vec, "e")) + $(_tex(:vec, "a"))^{$(_tex(:transp))}$(_tex(:vec, "e"))",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))+$(_tex(:vec, "e"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
)).
```
"""
function ManifoldsBase.exp(G::HeisenbergGroup, g, X)
    h = similar(X)
    exp!(G, h, g, X)
    return h
end

function ManifoldsBase.exp!(G::HeisenbergGroup, h, g, X)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    copyto!(h, I)
    a_g_view = _heisenberg_a_view(G, g)
    b_g_view = _heisenberg_b_view(G, g)
    a_X_view = _heisenberg_a_view(G, X)
    b_X_view = _heisenberg_b_view(G, X)
    h[1, 2:(n + 1)] .= a_g_view .+ a_X_view
    h[2:(n + 1), n + 2] .= b_g_view .+ b_X_view
    h[1, n + 2] =
        g[1, n + 2] + X[1, n + 2] + dot(a_X_view, b_X_view) / 2 + dot(a_g_view, b_X_view)
    return h
end

@doc raw"""
    injectivity_radius(G::HeisenbergGroup)

Return the injectivity radius on the [`HeisenbergGroup`](@ref) `G`, which is ``∞``.
"""
ManifoldsBase.injectivity_radius(::HeisenbergGroup) = Inf

@doc """
    log(G::HeisenbergGroup, g, h)

Compute the logarithmic map on the [`HeisenbergGroup`](@ref) group.

We denote two points ``g, h`` from the Heisenberg by

```math
g = $(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
))
$(_tex(:qquad))
h = $(_tex(:pmatrix,
    "1 & $(_tex(:vec, "d"))^{$(_tex(:transp))} & f",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "e"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
)),
```

where ``I_n`` is the ``n×n`` unit matrix, ``$(_tex(:vec, "a")), $(_tex(:vec, "b")), $(_tex(:vec, "d")), $(_tex(:vec, "e")) ∈ ℝ^n`` are vectors of length ``n``,
``$(_tex(:vec, 0))_n`` is the zero vector of length ``n``, and ``c,f ∈ ℝ`` are real numbers.

Then formula reads
```math
$(_tex(:log))_g(h) = $(_tex(:pmatrix,
    "0 & ($(_tex(:vec, "d"))-$(_tex(:vec, "q")))^{$(_tex(:transp))} & f - c + $(_tex(:vec, "a"))^{$(_tex(:transp))}$(_tex(:vec, "b")) - $(_tex(:vec, "d"))^{$(_tex(:transp))}$(_tex(:vec, "e")) - $(_tex(:frac,"1","2"))($(_tex(:vec, "d"))-$(_tex(:vec, "a")))^{$(_tex(:transp))}($(_tex(:vec, "e"))-$(_tex(:vec, "b")))",
    "$(_tex(:vec, 0))_n & Z_n & $(_tex(:vec, "e")) - $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 0"
)),
```
where additionally ``Z_n`` denotes the ``n×n`` zero matrix.
"""
ManifoldsBase.log(::HeisenbergGroup, g, h)

function ManifoldsBase.log!(G::HeisenbergGroup, X, g, h)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    fill!(X, 0)
    a_g_view = _heisenberg_a_view(G, g)
    b_g_view = _heisenberg_b_view(G, g)
    a_h_view = _heisenberg_a_view(G, h)
    b_h_view = _heisenberg_b_view(G, h)
    X[1, 2:(n + 1)] .= a_h_view .- a_g_view
    X[2:(n + 1), n + 2] .= b_h_view .- b_g_view
    pinvq_c = dot(a_g_view, b_g_view) - g[1, n + 2] + h[1, n + 2] - dot(a_g_view, b_h_view)
    X[1, n + 2] = pinvq_c - dot(a_h_view - a_g_view, b_h_view - b_g_view) / 2
    return X
end
function ManifoldsBase.log!(
    ::HeisenbergGroup,
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
)
    fill!(X, 0)
    return X
end

@doc """
    log(G::HeisenbergGroup, g)
    log!(G::HeisenbergGroup, X, g)

Compute the Lie group logarithm for the [`HeisenbergGroup`](@ref) `G`.

For ``g = $(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c",
    "$(_tex(:vec, 0))_n & I_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
))``
from the Lie algebra of the Heisenberg group,
where ``$(_tex(:vec, "a")), $(_tex(:vec, "b")) ∈ ℝ^n`` vectors of length ``n``,
``$(_tex(:vec, 0))_n`` is the zero vector of length ``n``, ``c ∈ ℝ``, and
``I_n`` is the ``n×n`` unit matrix.

Then the

```math
$(_tex(:log))_{$(_tex(:Cal, "G"))}(g) =
$(_tex(:pmatrix,
    "1 & $(_tex(:vec, "a"))^{$(_tex(:transp))} & c - $(_tex(:frac,"1","2"))$(_tex(:vec, "a"))^{$(_tex(:transp))}$(_tex(:vec, "b"))",
    "$(_tex(:vec, 0))_n & Z_n & $(_tex(:vec, "b"))",
    "0 & $(_tex(:vec, 0))_n^{$(_tex(:transp))} & 1"
)),
```
where ``Z_n`` denotes the ``n×n`` zero matrix.

This can be computed in-place of the Lie algebra vector `X`.
"""
ManifoldsBase.log(G::HeisenbergGroup, g)

function ManifoldsBase.log!(G::HeisenbergGroup, X, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    fill!(X, 0)
    # Obtain views to parts of x
    view_a_X = _heisenberg_a_view(G, X)
    view_b_X = _heisenberg_b_view(G, X)
    # Set then to views of g
    view_a_X .= _heisenberg_a_view(G, g)
    view_b_X .= _heisenberg_b_view(G, g)
    # Set first row last entry – since these contain g in X
    X[1, n + 2] = g[1, n + 2] - dot(view_a_X, view_b_X) / 2
    return X
end
function ManifoldsBase.log!(
    ::HeisenbergGroup, X, ::Identity{MatrixMultiplicationGroupOperation}
)
    fill!(X, 0)
    return X
end

function Base.show(
    io::IO, ::HeisenbergGroup{ManifoldsBase.TypeParameter{Tuple{n}}}
) where {n}
    return print(io, "HeisenbergGroup($(n))")
end
function Base.show(io::IO, G::HeisenbergGroup{Tuple{Int}})
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return print(io, "HeisenbergGroup($(n); parameter=:field)")
end
