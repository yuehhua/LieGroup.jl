"""
    LieGroupsTestSuite.jl

This module provides tools and dummy structures to test functionality provided
within `LieGroups.jl`.

For every test function, several interactions to other functions can be activated.
The following functions are expected to be available, since their defaults just pass through to the manifold
* `is_point` both on the Lie group `G` and the Lie algebra `𝔤`
* `isapprox(G, g, h)` and `isapprox(𝔤, X, Y)`
* `copy(G, g)`
* `norm(𝔤, X)`
"""
module LieGroupsTestSuite
using LieGroups
using Test, Random
using LinearAlgebra: I

include("so4_edge_cases.jl")

#
#
# === Dummy Types ===
struct DummyOperation <: AbstractGroupOperation end
struct DummySecondOperation <: AbstractGroupOperation end
struct DummyManifold <: LieGroups.AbstractManifold{LieGroups.ℝ} end
struct DummyActionType <: AbstractGroupActionType end
struct DummyLeftActionType <: AbstractLeftGroupActionType end
struct DummyRightActionType <: AbstractRightGroupActionType end
const DummyLieGroup = LieGroup{LieGroups.ℝ,DummyOperation,DummyManifold}
LieGroups.switch(a::DummyActionType) = a
LieGroups.switch(::DummyLeftActionType) = DummyRightActionType()
LieGroups.switch(::DummyRightActionType) = DummyLeftActionType()

"""
    rotate_matrix(n, k1, k2, α)

Generate a rotation matrix in ``ℝ^{n×n}`` with a rotation about ``α`` (in radians)
in the ``k_1-k_2``-plane.
"""
function rotation_matrix(n, k1, k2, α)
    R = Matrix{Float64}(I, n, n)
    R[k1, k1] = cos(α)
    R[k2, k2] = cos(α)
    R[k1, k2] = sin(α)
    R[k2, k1] = -sin(α)
    return R
end

# === Test single functions ===
#
#
# --- A
"""
    test_adjoint(G::LieGroup, g, X; kwargs...)

Test  `adjoint` function for a given Lie group element `g` and a Lie Algebra vector `X`

# Keyword arguments
* `expected=missing` provide the value expected. If none is provided, the
  default from `diff_conjugate` is used
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_adjoint(G::LieGroup, g, X; expected=missing, test_mutating::Bool=true)
    @testset "adjoint" begin
        v = if ismissing(expected)
            diff_conjugate(G, g, identity_element(G, typeof(X)), X)
        else
            expected
        end
        𝔤 = LieAlgebra(G)
        Y1 = adjoint(G, g, X)
        @test is_point(𝔤, Y1)
        if test_mutating
            Y2 = copy(𝔤, X)
            adjoint!(G, Y2, g, X)
            @test isapprox(𝔤, Y1, Y2)
        end
        @test isapprox(𝔤, Y1, v)
    end
    return nothing
end

"""
    test_apply(A::GroupAction, g, p; expected=missing)

Test  `apply`.

# Keyword arguments
* `expected=missing`: the result of the application of the group action.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_apply(A::GroupAction, g, p; expected=missing, test_mutating::Bool=true)
    @testset "apply" begin
        q1 = apply(A, g, p)
        M = base_manifold(A)
        @test is_point(M, q1)
        if test_mutating
            q2 = copy(M, p)
            apply!(A, q2, g, p)
            @test isapprox(M, q1, q2)
        end
        !ismissing(expected) && @test isapprox(M, q1, expected)
    end
end
#
#
# --- C
"""
    test_compose(G::LieGroup, g, h; kwargs...)

Test  `compose` for given Lie group elements `g`, `h`.

# Keyword arguments

* `atol::Real=0`: the absolute tolerance for the tests of zero-vectors
`
* `test_identity::Bool=true`: test that composing with the identity yields the identity (requires `identity_element`)
* `test_inverse::Bool=true`: test that `g^{-1}g` is the identity (requires `inv`, `inv!`, and `is_identity`)
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_compose(
    G::LieGroup,
    g,
    h;
    atol::Real=0,
    test_inverse::Bool=true,
    test_identity::Bool=true,
    test_mutating::Bool=true,
)
    @testset "compose" begin
        k1 = compose(G, g, h)
        k2 = copy(G, g)
        compose!(G, k2, g, h)
        @test isapprox(G, k1, k2)
        if test_inverse
            for g_ in [g, h]
                g_inv = inv(G, g_)
                k1 = compose(G, g_inv, g_)
                @test is_identity(G, k1; atol=atol)
                if test_mutating
                    compose!(G, k2, g_inv, g_)
                    @test isapprox(G, k1, k2; atol=atol)
                    @test is_identity(G, k2; atol=atol)
                end
            end
        end
        if test_identity && test_mutating
            for g_ in [g, h]
                for e in [Identity(G), identity_element(G, typeof(g))]
                    k1 = compose(G, g_, e)
                    compose!(G, k2, g_, e)
                    @test isapprox(G, k1, k2; atol=atol)
                    k1 = compose(G, e, g_)
                    compose!(G, k2, e, g_)
                    @test isapprox(G, k1, k2; atol=atol)
                end
            end
            e = Identity(G)
            k3 = copy(G, g)
            compose!(G, k3, e, e)
            @test is_identity(G, k3; atol=atol)
        end
    end
    return nothing
end

"""
    test_conjugate(G::LieGroup, g, h; expected=missing)

Test  `conjugate`.

# Keyword arguments
* `expected=missing`: the result of the conjugate with a concrete expected result.
* `test_default=false`: Since the conjugate can be computed using
  `compose` and `inv` – this test can check that this default provides the same result
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_conjugate(
    G::LieGroup, g, h; expected=missing, test_default::Bool=true, test_mutating::Bool=true
)
    @testset "conjugate" begin
        k1 = conjugate(G, g, h)
        @test is_point(G, k1, true)
        if test_mutating
            k2 = copy(G, g)
            conjugate!(G, k2, g, h)
            @test isapprox(G, k1, k2)
        end
        if !ismissing(expected)
            @test isapprox(G, k1, expected)
        end
        if test_default
            @test isapprox(G, k1, compose(G, g, compose(G, h, inv(G, g))))
        end
    end
    return nothing
end

#
#
# --- D
"""
    test_diff_apply(A::GroupAction, g, p, X; expected=missing)

Test  `diff_apply`.

# Keyword arguments
* `expected=missing`: the result of the application of the group action.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_diff_apply(
    A::GroupAction, g, p, X; expected=missing, test_mutating::Bool=true
)
    @testset "diff_apply" begin
        Y1 = diff_apply(A, g, p, X)
        M = base_manifold(A)
        q = apply(A, g, p)
        G = base_lie_group(A)
        is_vector(G, q, Y1)
        if test_mutating
            Y2 = copy(M, p, X)
            diff_apply!(A, Y2, g, p, X)
            @test isapprox(M, q, Y1, Y2)
        end
        !ismissing(expected) && @test isapprox(M, q, Y1, expected)
    end
end

"""
    test_diff_apply(A::GroupAction, g, p, X; expected=missing)

Test  `diff_group_apply`.

# Keyword arguments
* `expected=missing`: the result of the application of the group action.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_diff_group_apply(
    A::GroupAction, g, p, X; expected=missing, test_mutating::Bool=true
)
    @testset "diff_group_apply" begin
        Y1 = diff_group_apply(A, g, p, X)
        G = base_lie_group(A)
        @test is_vector(G, g, Y1; error=:error)
        𝔤 = LieAlgebra(G)
        if test_mutating
            Y2 = copy(𝔤, X)
            diff_group_apply!(A, Y2, g, p, X)
            @test isapprox(𝔤, Y1, Y2)
        end
        !ismissing(expected) && @test isapprox(𝔤, Y1, expected)
    end
end

"""
    test_diff_inv(G::LieGroup, g, X; expected=missing)

Test  `diff_inv`.

# Keyword arguments
* `expected=missing`: the result of the differential of the inverse, if not provided,
  only consistency between the allocating and the in-place variant is checked.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_diff_inv(G::LieGroup, g, X; expected=missing, test_mutating::Bool=true)
    @testset "diff_inv" begin
        𝔤 = LieAlgebra(G)
        # Check that at identity it is in the Lie algebra
        Ye = diff_inv(G, identity_element(G, typeof(X)), X)
        @test is_point(𝔤, Ye; error=:error)
        Y1 = diff_inv(G, g, X)
        if test_mutating
            Y2 = zero_vector(𝔤, typeof(X))
            Y2 = diff_inv!(G, Y2, g, X)
            @test isapprox(G, g, Y1, Y2)
        end
        if !ismissing(expected)
            @test isapprox(𝔤, Y1, expected)
        end
    end
end

"""
    test_diff_left_compose(G::LieGroup, g, h, X; expected=missing)

Test  `diff_left_compose`.

# Keyword arguments
* `expected=missing`: the result of the differential of the compose's left argument,
  if not provided, only consistency between the allocating and the in-place variant is checked.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_diff_left_compose(
    G::LieGroup, g, h, X; expected=missing, test_mutating::Bool=true
)
    @testset "diff_left_compose" begin
        𝔤 = LieAlgebra(G)
        Ye = diff_left_compose(G, identity_element(G, typeof(X)), h, X)
        @test is_point(𝔤, Ye; error=:error)
        # check that in-place and allocating agree
        Y1 = diff_left_compose(G, g, h, X)
        if test_mutating
            Y2 = zero_vector(𝔤, typeof(X))
            diff_left_compose!(G, Y2, g, h, X)
            @test isapprox(LieAlgebra(G), Y1, Y2)
        end
        if !ismissing(expected)
            @test isapprox(LieAlgebra(G), Y1, expected)
        end
    end
end

"""
    test_diff_right_compose(G::LieGroup, g, h, X; expected=missing)

Test  `diff_right_compose`.

# Keyword arguments
* `expected=missing`: the result of the differential of the compose's right argument,
  if not provided, only consistency between the allocating and the in-place variant is checked.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_diff_right_compose(
    G::LieGroup, g, h, X; expected=missing, test_mutating::Bool=true
)
    @testset "diff_right_compose" begin
        𝔤 = LieAlgebra(G)
        Ye = diff_right_compose(G, identity_element(G, typeof(X)), h, X)
        @test is_point(𝔤, Ye; error=:error)
        # check that in-place and allocating agree
        Y1 = diff_right_compose(G, g, h, X)
        if test_mutating
            Y2 = zero_vector(𝔤, typeof(X))
            diff_right_compose!(G, Y2, g, h, X)
            @test isapprox(𝔤, Y1, Y2)
        end
        if !ismissing(expected)
            @test isapprox(𝔤, Y1, expected)
        end
    end
end

"""
    test_copyto(G::LieGroup, g)

Test that `copyto!` works also when copying over an `Identity`.

The point `g` can be any point _but_ the `identity_element`.
The group has to be a mutating one, that is, not work on isbit types.
"""
function test_copyto(G::LieGroup, g)
    @testset "copyto!" begin
        k = copy(G, g)
        e = Identity(G)
        copyto!(G, k, e)
        # also check that this holds both ways
        @test isapprox(G, k, Identity(G))
        @test isapprox(G, Identity(G), k)
        # Test that copying back also works
        copyto!(G, k, g)
        @test isapprox(G, k, g)
        # copy into identity only works if provided the identity
        @test copyto!(G, e, Identity(G)) == e
        # but then also always returns e
        copyto!(G, k, e)
        @test copyto!(G, e, k) == e
        # and fails if the point to copy in is not (numerically) e
        @test_throws DomainError copyto!(G, e, g)
        return nothing
    end
end

#
#
# --- D
"""
    test_diff_conjugate(A::GroupAction, g, p, X; expected=missing)

Test  `diff_conjugate`
"""
function test_diff_conjugate(
    G::LieGroup, g, h, X; expected=missing, test_mutating::Bool=true
)
    𝔤 = LieAlgebra(G)
    @testset "diff_conjugate" begin
        Y1 = diff_conjugate(G, g, h, X)
        @test is_point(𝔤, Y1; error=:error)
        if test_mutating
            Y2 = zero_vector(𝔤, typeof(X))
            diff_conjugate!(G, Y2, g, h, X)
            @test isapprox(𝔤, Y1, Y2)
        end
        if !ismissing(expected)
            @test isapprox(𝔤, Y1, expected)
        end
        return nothing
    end
end
#
#
# --- E
"""
    test_exp_log(G::LieGroup, g, h, X)

Test  `exp` and `log` for given Lie group elements `g`, `h` and
a vector `X` from the Lie Algebra.

!!! note
    This function requires `is_point(G, g)` and `is_point(LieAlgebra(G), X)` to be implemented

# Keyword arguments

* `atol::Real=0`: the absolute tolerance for the tests of zero-vectors
* `test_exp::Bool=true`: test the exponential map yields a point on `G`
* `test_log::Bool=true`: test the logarithmic map.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_exp_log(
    G::LieGroup,
    g,
    h,
    X;
    atol::Real=0,
    test_exp::Bool=true,
    test_mutating::Bool=true,
    test_log::Bool=true,
)
    @testset "(Lie group) exp & log" begin
        𝔤 = LieAlgebra(G)
        e = Identity(G)
        if test_exp
            # Lie group exp
            k1 = exp(G, X)
            if test_mutating
                k2 = copy(G, g)
                exp!(G, k2, X)
                @test isapprox(G, k1, k2)
            end
            @test is_point(G, k1; error=:error)
            # exp
            k1 = exp(G, g, X)
            if test_mutating
                k2 = copy(G, g)
                exp!(G, k2, g, X)
                @test isapprox(G, k1, k2)
            end
            @test is_point(G, k1; error=:error)
        end
        if test_log
            # Lie group log
            Y1 = log(G, g)
            if test_mutating
                Y2 = zero_vector(𝔤, typeof(X))
                log!(G, Y2, g)
                @test isapprox(𝔤, Y1, Y2; atol=atol)
                log!(G, Y2, e)
                @test isapprox(𝔤, Y2, 0 * Y2; atol=atol)
            end
            @test is_point(𝔤, Y1; error=:error)
            @test norm(𝔤, log(G, g, g)) ≈ 0 atol = atol
            @test norm(𝔤, log(G, h, h)) ≈ 0 atol = atol
            # log
            Y1 = log(G, g, h)
            if test_mutating
                Y2 = copy(G, g)
                log!(G, Y2, g, h)
                @test isapprox(𝔤, Y1, Y2)
            end
            @test is_point(𝔤, Y1; error=:error)
            Y3 = zero_vector(𝔤, typeof(X))
            @test isapprox(𝔤, Y3, log(G, e, typeof(X)); atol=atol)
            log!(G, Y3, e)
            @test isapprox(G, e, Y3, log(G, e, typeof(Y3)); atol=atol)
            @test isapprox(𝔤, log(G, g, g), Y3; atol=atol)
            @test isapprox(𝔤, log(G, h, h), Y3; atol=atol)
        end
        if test_exp && test_log
            # Lie group exp / log
            k1 = exp(G, X)
            Y1 = log(G, k1)
            @test isapprox(𝔤, X, Y1)
            if test_mutating
                k2 = copy(G, g)
                exp!(G, k2, X)
                Y2 = copy(G, g)
                log!(G, Y2, k2)
                @test isapprox(𝔤, Y1, Y2)
            end
            # exp & log
            k1 = exp(G, g, X)
            Y1 = log(G, g, k1)
            @test isapprox(𝔤, X, Y1)
            if test_mutating
                k2 = copy(G, g)
                exp!(G, k2, g, X)
                Y2 = copy(G, g)
                log!(G, Y2, g, k2)
                @test isapprox(𝔤, Y1, Y2)
            end
        end
    end
    return nothing
end

#
#
# --- H
"""
    test_hat_vee(G::LieGroup, g, X; kwargs...)

Test `hat` and `vee` for given Lie group element `g` and a Lie Algebra vector `X`.

# Keyword arguments

* `expected_value=missing`: the expected value of `vee(G,X)`
* `test_mutating::Bool=true`: test the mutating functions
* `test_vee::Bool=true`: test the vee function
* `test_hat::Bool=true`: test the hat function
"""
function test_hat_vee(
    G::LieGroup,
    g,
    X;
    test_mutating::Bool=true,
    test_vee::Bool=true,
    test_hat::Bool=true,
    expected_value=missing,
)
    @testset "hat & vee" begin
        𝔤 = LieAlgebra(G)
        if test_hat
            c = ismissing(expected_value) ? zeros(manifold_dimension(G)) : expected_value
            Y1 = hat(𝔤, c, typeof(X))
            @test is_vector(G, g, Y1; error=:error)
            !ismissing(expected_value) && @test isapprox(𝔤, X, Y1)
            if test_mutating
                Y2 = zero_vector(𝔤, typeof(Y1))
                hat!(𝔤, Y2, c)
                @test isapprox(𝔤, Y1, Y2)
            end
        end
        if test_vee
            c1 = vee(𝔤, X)
            if test_mutating
                c2 = zero(c1)
                vee!(𝔤, c2, X)
                @test c1 ≈ c2
            end
            if !ismissing(expected_value)
                @test c1 ≈ expected_value
            end
        end
        if test_hat && test_vee
            Y1 = hat(𝔤, vee(𝔤, X), typeof(X))
            @test isapprox(𝔤, X, Y1)
            if test_mutating
                Y2 = zero_vector(𝔤, typeof(Y1))
                c = zeros(manifold_dimension(G))
                vee!(𝔤, c, X)
                hat!(𝔤, Y2, c)
                @test isapprox(𝔤, Y1, Y2)
            end
        end
    end
    return nothing
end

#
#
# --- `I`

"""
    test_injectivity_radius(G::LieGroup; kwargs...)

Test the function `injectivity_radius`.

# Keyword arguments

* `expected=missing`: expected value for global injectivity radius.
"""
function test_injectivity_radius(G::LieGroup; expected=missing)
    @testset "injectivity radius" begin
        if ismissing(expected)
            @test injectivity_radius(G) isa Real
            @test injectivity_radius(G) >= 0
        else
            @test injectivity_radius(G) == expected
        end
    end
    return nothing
end

"""
    test_inv_compose(G::LieGroup, g, h, X; kwargs...)

Test the special functions combining inv and compose, `inv_left_compose` and `inv_right_compose`.
For these tests both `compose` and `inv` are required.

# Keyword arguments

* `test_left::Bool=true`: test ``g^{-1}∘h``
* `test_mutating::Bool=true`: test the mutating functions
* `test_right::Bool=true`: test ``g∘h^{-1}``
"""
function test_inv_compose(
    G::LieGroup,
    g,
    h;
    expected_left=missing,
    expected_right=missing,
    test_left::Bool=true,
    test_mutating::Bool=true,
    test_right::Bool=true,
)
    @testset "test compose inv combinations" begin
        if test_left
            v = if ismissing(expected_left)
                compose(G, inv(G, g), h)
            else
                expected_left
            end
            @testset "g^{-1}∘h" begin
                k1 = inv_left_compose(G, g, h)
                @test isapprox(G, k1, v)
                if test_mutating
                    k2 = copy(G, g)
                    inv_left_compose!(G, k2, g, h)
                    @test isapprox(G, k1, k2)
                end
            end
        end
        if test_right
            v = if ismissing(expected_right)
                compose(G, g, inv(G, h))
            else
                expected_right
            end
            @testset "g∘h^{-1}" begin
                k1 = inv_right_compose(G, g, h)
                @test isapprox(G, k1, v)
                if test_mutating
                    k2 = copy(G, g)
                    inv_right_compose!(G, k2, g, h)
                    @test isapprox(G, k1, k2)
                end
            end
        end
    end
    return nothing
end

"""
    test_inv(G::LieGroup, g)

Test the inverse function, both the allocating and the in-place variant,
and that the double inverse is the identity.

# Keyword arguments

* `test_mutating::Bool=true`: test the mutating functions
* `test_identity::Bool=true`: test that `inv(e) == e`
"""
function test_inv(G::LieGroup, g; test_mutating::Bool=true, test_identity::Bool=true)
    @testset "inv" begin
        k1 = inv(G, g)
        @test is_point(G, k1)
        g1 = inv(G, k1)
        @test isapprox(G, g, g1)
        if test_mutating
            k2 = copy(G, g)
            inv!(G, k2, g)
            @test isapprox(G, k1, k2)
            # continue in-place
            inv!(G, k2, k2)
            @test isapprox(G, k2, g)
        end
        if test_identity
            e = Identity(G)
            @test inv(G, e) === e
            if test_mutating
                e2 = copy(G, g)
                inv!(G, e2, e)
                @test is_identity(G, e2)
                e3 = copy(G, g)
                inv!(G, e3, e) # materialize identity
                @test is_identity(G, e3)
            end
        end
    end
    return nothing
end

"""
    test_is_identity(G::LieGroup, g)

Test that the `Identity` returns that `is_identity` is true and that it is a point
"""
function test_identity(G::LieGroup)
    @testset "Identity" begin
        e = Identity(G)
        @test is_point(G, e; error=:error)
        @test is_identity(G, e)
        e2 = Identity(DummyOperation)
        @test !is_point(G, e2)
        @test_throws DomainError !is_point(G, e2; error=:error)
        @test !is_identity(G, e2)
    end
    return nothing
end
#
#
# --- L
"""
    test_lie_bracket(G::LieGroup, X, Y; expected=missing)

Test  `lie_bracket`.

# Keyword arguments
* `expected=missing`: the result of the lie bracket
  if not provided, only consistency between the allocating and the in-place variant is checked.
* `test_mutating::Bool=true`: test the mutating functions
"""
function test_lie_bracket(G::LieGroup, X, Y; expected=missing, test_mutating::Bool=true)
    @testset "lie_bracket" begin
        𝔤 = LieAlgebra(G)
        Z1 = lie_bracket(𝔤, X, Y)
        if test_mutating
            Z2 = copy(𝔤, X)
            lie_bracket!(𝔤, Z2, X, Y)
            @test isapprox(𝔤, Z1, Z2)
        end
        if !ismissing(expected)
            @test isapprox(𝔤, Z1, expected)
        end
    end
end

#
#
# --- R
"""
    test_rand(G::LieGroup)

Test the random function, both the allocating and the in-place variant,
as well as the variant with an `rng`, if one is provided.

both the random point and the random tangent vector variants are tested.

# Keyword arguments

* `test_mutating::Bool=true`: test the mutating functions
* `rng=missing`: test with a specific random number generator
"""
function test_rand(
    G::LieGroup, g; test_mutating::Bool=true, rng::Union{Missing,AbstractRNG}=missing
)
    @testset "rand" begin
        g1 = rand(G)
        @test is_point(G, g1; error=:error)
        if test_mutating
            g2 = copy(G, g)
            rand!(G, g2)
            @test is_point(G, g2; error=:error)
        end
        X1 = rand(G; vector_at=g1)
        @test is_vector(G, g1, X1; error=:error)
        if test_mutating
            X2 = zero_vector(LieAlgebra(G), typeof(g1))
            rand!(G, X2; vector_at=g1)
            @test is_vector(G, g1, X2; error=:error)
        end
        if !ismissing(rng)
            g1 = rand(rng, G)
            @test is_point(G, g1; error=:error)
            if test_mutating
                g2 = copy(G, g)
                rand!(rng, G, g2)
                @test is_point(G, g2; error=:error)
            end
            X1 = rand(rng, G; vector_at=g1)
            @test is_vector(G, g1, X1; error=:error)
            if test_mutating
                X2 = zero_vector(LieAlgebra(G), typeof(X1))
                rand!(rng, G, X2; vector_at=g1)
                @test is_vector(G, g1, X2; error=:error)
            end
        end
    end
    return nothing
end

#
#
# --- S
"""
    test_show(G, repr_string::AbstractString)

Test that show methods work as expected.
For now this (only) checks that `"\$G"` yields the `repr_string`.

Requires `show` (or `repr`) to be implemented.
"""
function test_show(G::Union{GroupAction,LieGroup}, repr_string::AbstractString)
    @testset "repr(G, g, h)" begin
        @test repr(G) == repr_string
    end
    return nothing
end

# The global test function for a Lie group
#
#
"""
    test_lie_group(G::LieGroup, properties::Dict, expectations::Dict)

Test the Lie group ``G`` based on a `Dict` of properties and a `Dict` of `expectations

Possible properties are

* `:Functions` is a vector of all defined functions for `G`
  Note that if `f` is in `:Functions`, and `f!` makes sense, for example for `compose`,
  it is assumed that both are defined.
* `:Points` is a vector of at least 2 points on `G`, the first is not allowed to be the identity numerically
* `:Vectors` is a vector of at least 2 elements from the Lie algebra `𝔤` og `G`
* `:Mutating` is a boolean (`true` by default) whether to test the mutating variants of functions or not.
* `:Name` is a name of the test. If not provided, defaults to `"\$G"`
* `:Rng` is a random number generator, if provided, the random functions are tested with this generator as well

Possible `expectations` are

* `:adjoint` for the result of `conjgate` in the case where `diff_conjugate` is not implemented
* `:atol => 0.0` a global absolute tolerance
* `:atols -> Dict()` a dictionary `function -> atol` for specific function tested.
* `:conjugate` for the result of `conjgate
* `:conjugate_default => false` to activate the test of the default implementation of `conjugate`
* `:diff_inv` for the result of `diff_inv` with respect to the first point and the first vector.
* `:diff_left_compose` for the result of `diff_left_compose` with respect to the first two points and the first vector.
* `:diff_right_compose` for the result of `diff_right_compose` with respect to the first two points and the first vector.
* `:inv_left_compose` for the result of `inv_left_right_compose` with respect to the first two points
* `:inv_right_compose` for the result of `inv_right_compose` with respect to the first two points
* `:repr` is a sting one gets from `repr(G)`
* `:vee` for the result of `vee(G, X)` where `X` is the first of the vectors
"""
function test_lie_group(G::LieGroup, properties::Dict, expectations::Dict=Dict())
    atol = get(expectations, :atol, 0.0)
    mutating = get(properties, :Mutating, true)
    functions = get(properties, :Functions, Function[])
    points = get(properties, :Points, [])
    vectors = get(properties, :Vectors, [])
    test_name = get(properties, :Name, "$G")
    function_atols = get(expectations, :atols, Dict())
    @testset "$(test_name)" begin
        # Call function tests based on their presence in alphabetical order
        #
        #
        # --- A
        if (adjoint in functions)
            v = get(expectations, :adjoint, missing)
            test_adjoint(G, points[1], vectors[1]; expected=v, test_mutating=mutating)
        end
        #
        #
        # --- C
        if (compose in functions)
            ti = all(in.([inv, is_identity], Ref(functions)))
            compose_atol = get(function_atols, :compose, 0.0)
            identity_atol = get(function_atols, :is_identity, 0.0)
            local_atol = max(compose_atol, identity_atol, atol)
            test_compose(
                G,
                points[1],
                points[2];
                test_inverse=ti,
                test_mutating=mutating,
                atol=local_atol,
            )
        end
        # since there is a default, also providing compose&inv suffices
        if (conjugate in functions) || (all(in.([compose, inv], Ref(functions))))
            v = get(expectations, :conjugate, missing)
            td = get(expectations, :conjugate_default, false)
            test_conjugate(
                G, points[1], points[2]; expected=v, test_mutating=mutating, test_default=td
            )
        end
        # Either `copyto` or the default with `identity_element` available
        if any(in.([copyto!, identity_element], Ref(functions))) && (mutating)
            test_copyto(G, points[1])
        end
        #
        #
        # --- D
        if (diff_conjugate in functions)
            v = get(expectations, :diff_conjugate, missing)
            test_diff_conjugate(
                G, points[1], points[2], vectors[1]; expected=v, test_mutating=mutating
            )
        end

        if (diff_inv in functions)
            v = get(expectations, :diff_inv, missing)
            test_diff_inv(G, points[1], vectors[1]; expected=v, test_mutating=mutating)
        end

        if (diff_left_compose in functions)
            v = get(expectations, :diff_left_compose, missing)
            test_diff_left_compose(
                G, points[1], points[2], vectors[1]; expected=v, test_mutating=mutating
            )
        end

        if (diff_right_compose in functions)
            v = get(expectations, :diff_right_compose, missing)
            test_diff_right_compose(
                G, points[1], points[2], vectors[1]; expected=v, test_mutating=mutating
            )
        end

        #
        #
        # --- E
        if any(in.([exp, log], Ref(functions)))
            exp_atol = get(function_atols, :exp, 0.0)
            log_atol = get(function_atols, :log, 0.0)
            local_atol = max(exp_atol, log_atol, atol)
            test_exp_log(
                G,
                points[1],
                points[2],
                vectors[1];
                atol=local_atol,
                test_exp=(exp in functions),
                test_log=(log in functions),
                test_mutating=mutating,
            )
        end

        #
        #
        # --- H
        if any(in.([hat, vee], Ref(functions)))
            v = get(expectations, :vee, missing)
            test_hat_vee(
                G,
                points[1],
                vectors[1];
                expected_value=v,
                test_hat=(hat in functions),
                test_mutating=mutating,
                test_vee=(vee in functions),
            )
        end

        #
        #
        # --- `I`
        if any(in.([inv_left_compose, inv_right_compose], Ref(functions)))
            vl = get(expectations, :inv_left_compose, missing)
            vr = get(expectations, :inv_right_compose, missing)
            test_inv_compose(
                G,
                points[1],
                points[2];
                expected_left=vl,
                expected_right=vr,
                test_left=(inv_left_compose in functions),
                test_mutating=mutating,
                test_right=(inv_right_compose in functions),
            )
        end
        if (injectivity_radius in functions)
            ir = get(expectations, :injectivity_radius, missing)
            test_injectivity_radius(G; expected=ir)
        end
        if (inv in functions)
            test_inv(G, points[1]; test_mutating=mutating)
        end
        if (is_identity in functions)
            test_identity(G)
        end
        #
        #
        # --- L
        if (lie_bracket in functions)
            v = get(expectations, :lie_bracket, missing)
            test_lie_bracket(G, vectors[1], vectors[2]; expected=v, test_mutating=mutating)
        end

        #
        #
        # --- R
        if (rand in functions)
            v = get(properties, :Rng, missing)
            test_rand(G, points[1]; rng=v, test_mutating=mutating)
        end

        #
        #
        # --- S
        if (any(in.([show, repr], Ref(functions)))) && haskey(expectations, :repr)
            test_show(G, expectations[:repr])
        end
    end
end

"""
    test_group_action(G::LieGroup, properties::Dict, expectations::Dict)

Test the Lie group ``G`` based on a `Dict` of properties and a `Dict` of `expectations`.

Possible properties are

* `:AlgebraVectors` is a vector of at least 3 elements from the Lie algebra `𝔤` of `G`
* `:Functions` is a vector of all defined functions for `G`
  Note that if `f` is in `:Functions`, and `f!` makes sense, for example for `compose`,
  it is assumed that both are defined.
* `:GroupPoints` is a vector of at least three points on `G`, the first is not allowed to be the identity numerically
* `:ManifoldPoints` is a vector of at least three points on `M`
* `:TangentVectors` is a vector of at least three tangent vectors on `M`, each in the tangent space of the corresponding `:ManifoldPoint`
* `:Mutating` is a boolean (`true` by default) whether to test the mutating variants of functions or not.
* `:Name` is a name of the test. If not provided, defaults to `"\$G"`

Possible `expectations` are

* `:apply` for the result of `apply` on the first group and manifold point
* `:diff_apply` for the result of `apply` on the first group and manifold point together with the first tangent vector
* `:atol` a global absolute tolerance, defaults to `1e-8`
* `:group` is the `LieGroup` describing the action
* `:manifold` is the `AbstractManifold` the action acts upon
* `:repr` is a sting one gets from `repr(G)`
"""
function test_group_action(A::GroupAction, properties::Dict, expectations::Dict=Dict())
    a_tol = get(expectations, :atol, 1e-8)
    mutating = get(properties, :Mutating, true)
    functions = get(properties, :Functions, Function[])
    group_points = get(properties, :GroupPoints, [])
    @assert length(group_points) > 2
    algebra_vectors = get(properties, :AlgebraVectors, [])
    @assert length(algebra_vectors) > 2
    manifold_points = get(properties, :ManifoldPoints, [])
    @assert length(manifold_points) > 2
    tangent_vectors = get(properties, :TangentVectors, [])
    @assert length(tangent_vectors) > 2
    test_name = get(properties, :Name, "$A")
    @testset "$(test_name)" begin
        # Call function tests based on their presence in alphabetical order
        #
        #
        # --- A
        if (apply in functions)
            v = get(expectations, :apply, missing)
            test_apply(
                A, group_points[1], manifold_points[1]; expected=v, test_mutating=mutating
            )
        end
        if (diff_apply in functions)
            v = get(expectations, :diff_apply, missing)
            test_diff_apply(
                A,
                group_points[1],
                manifold_points[1],
                tangent_vectors[1];
                expected=v,
                test_mutating=mutating,
            )
        end
        if (diff_group_apply in functions)
            v = get(expectations, :diff_group_apply, missing)
            test_diff_group_apply(
                A,
                group_points[1],
                manifold_points[1],
                algebra_vectors[1];
                expected=v,
                test_mutating=mutating,
            )
        end
        if (base_lie_group in functions)
            v = get(expectations, :group, missing)
            !ismissing(v) && @testset "base_lie_group" begin
                @test base_lie_group(A) == v
            end
        end
        if (base_manifold in functions)
            v = get(expectations, :manifold, missing)
            !ismissing(v) && @testset "base_manifold" begin
                @test base_manifold(A) == v
            end
        end
        if (any(in.([show, repr], Ref(functions)))) && haskey(expectations, :repr)
            test_show(A, expectations[:repr])
        end
    end
end

export test_lie_group, test_group_action
end # module
