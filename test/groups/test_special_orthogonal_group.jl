using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite
using LieGroupsTestSuite: rotation_matrix

begin
    G = SpecialOrthogonalGroup(2)
    g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
    g2 = [0.0 -1.0; 1.0 0.0]
    g3 = [1.0 0.0; 0.0 1.0]
    X1, X2, X3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 0.2 0.0], [0.0 0.0; 0.0 0.0]
    properties = Dict(
        :Name => "The special orthogonal group SO(2)",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            adjoint,
            compose,
            conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            hat,
            identity_element,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectations = Dict(
        :repr => "SpecialOrthogonalGroup(2)", :lie_bracket => zero(X1), :atol => 1e-10
    )
    test_lie_group(G, properties, expectations)
    G2 = SpecialOrthogonalGroup(2; parameter=:field)
    test_lie_group(G2, properties, expectations)
    #
    #
    # SO(3)
    H = SpecialOrthogonalGroup(3)
    h1 = rotation_matrix(3, 2, 1, π / 4)
    h2 = rotation_matrix(3, 2, 1, π / 8) * rotation_matrix(3, 3, 1, π / 4)
    h3 = rotation_matrix(3, 3, 1, π / 4) * rotation_matrix(3, 3, 2, π / 8)
    Y1 = [0.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
    Y2 = [0.0 0.0 -0.2; 0.0 0.0 0.0; 0.2 0.0 0.0]
    Y3 = [0.0 0.3 0.0; -0.3 0.0 0.4; 0.0 -0.4 0.0]
    # Test only specialized functions here
    properties2 = Dict(
        :Name => "The special orthogonal group SO(3) – specialized functions",
        :Points => [h1, h2, h3],
        :Vectors => [Y1, Y2, Y3],
        :Functions => [exp, hat, log, show, vee],
    )
    expectations2 = Dict(
        :repr => "SpecialOrthogonalGroup(3)", :atols => Dict(:exp => 1e-15)
    )
    test_lie_group(H, properties2, expectations2)
    H2 = SpecialOrthogonalGroup(3; parameter=:field)
    test_lie_group(H2, properties2, expectations2)
    # One special case – maybe add to test suite?
    # Hit cosθ = -1
    e = Identity(H)
    g_e = [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
    X_e = [0.0 0.0 0.0; 0.0 0.0 -π; 0.0 π 0.0]
    @test isapprox(H, e, log(H, e, g_e), X_e)
    @test isapprox(H, exp(H, e, X_e), g_e)
    #
    #
    # SO(4)
    J = SpecialOrthogonalGroup(4)
    j1 = rotation_matrix(4, 3, 1, π / 4)
    j2 = rotation_matrix(4, 4, 1, π / 8) * rotation_matrix(4, 3, 2, π / 4)
    j3 = rotation_matrix(4, 4, 1, π / 4) * rotation_matrix(4, 4, 2, π / 8)
    Z1 = [0.0 0.1 0.0 0.0; -0.1 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    Z2 = [0.0 0.0 0.2 0.0; 0.0 0.0 0.0 0.0; -0.2 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    Z3 = [0.0 0.1 0.0 0.3; 0.0 0.0 -0.4 0.0; 0.0 0.4 0.0 0.0; -0.3 0.0 0.0 0.0]
    # Test only specialized functions here
    properties3 = Dict(
        :Name => "The orthogonal group SO(4) – specialized functions",
        :Points => [j1, j2, j3],
        :Vectors => [Z1, Z2, Z3],
        :Functions => [exp, hat, log, show, vee],
    )
    expectations3 = Dict(
        :repr => "SpecialOrthogonalGroup(4)", :atols => Dict(:exp => 1e-15)
    )
    test_lie_group(J, properties3, expectations3)
    @testset "𝔰𝔬(4) edge cases" begin
        e = Identity(J)
        X = zero_vector(J, e)
        d = vee(J, X)
        p = exp(J, e, X)
        for c in LieGroupsTestSuite.𝔰𝔬4_edges_cases_explog
            @testset "$c on $J" begin
                hat!(J, X, c)
                vee!(J, d, X)
                @test isapprox(c, d)
                l = exp(X)
                exp!(J, p, e, X)
                @test l ≈ p
                p2 = exp(J, e, log(J, e, p))
                @test isapprox(J, p, p2; atol=1e-6)
            end
        end
    end
end
