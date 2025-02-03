using LieGroups, Random, Test, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special Euclidean" begin
    G = SpecialEuclideanGroup(2)
    g1 = 1 / sqrt(2) .* [1.0 1.0 sqrt(2); -1.0 1.0 0.0; 0.0 0.0 sqrt(2)]
    g2 = [0.0 -1.0 0.0; 1.0 0.0 1.0; 0.0 0.0 1.0]
    g3 = [1.0 0.0 1.0; 0.0 1.0 1.0; 0.0 0.0 1.0]
    X1 = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
    X2 = [0.0 0.30 1.0; -0.30 0.0 1.0; 0.0 0.0 0.0]
    X3 = [0.0 0.1 1.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
    properties = Dict(
        :Name => "The special Euclidean group (affine matrices)",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            # adjoint,
            compose,
            # conjugate,
            # diff_inv,
            # diff_left_compose,
            # diff_right_compose,
            exp,
            hat,
            identity_element,
            inv,
            # inv_left_compose,
            # inv_right_compose,
            is_identity,
            # lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectations = Dict(
        :repr => "SpecialEuclideanGroup(2)",
        :atol => 1e-14,
        #:diff_inv => -X1,
        #:diff_left_compose => X1,
        #:diff_right_compose => X1,
        #:lie_bracket => zero(X1),
    )
    test_lie_group(G, properties, expectations)

    GL = SpecialEuclideanGroup(2)
    gL1 = ArrayPartition(1 / sqrt(2) * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
    gL2 = ArrayPartition([0.0 -1.0; 1.0 0.0], [0.0, 1.0])
    gL3 = ArrayPartition([1.0 0.0; 0.0 1.0], [1.0, 1.0])
    XL1 = ArrayPartition([0.0 -0.23; 0.23 0.0], [0.0, 1.0])
    XL2 = ArrayPartition([0.0 0.30; -0.30 0.0], [1.0, 1.0])
    XL3 = ArrayPartition([9.0 0.1; -0.1 0.0], [1.0, 0.0])
    propertiesL = Dict(
        :Name => "The special Euclidean group (left variant)",
        :Points => [gL1, gL2, gL3],
        :Vectors => [XL1, XL2, XL3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            # adjoint,
            compose,
            # conjugate,
            # diff_inv,
            # diff_left_compose,
            # diff_right_compose,
            exp,
            hat,
            identity_element,
            inv,
            # inv_left_compose,
            # inv_right_compose,
            is_identity,
            # lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectationsL = Dict(
        :repr => "SpecialEuclideanGroup(2)",
        :atol => 1e-14,
        #:diff_inv => -X1,
        #:diff_left_compose => X1,
        #:diff_right_compose => X1,
        #:lie_bracket => zero(X1),
    )
    test_lie_group(GL, propertiesL, expectationsL)
    #
    #
    # Right variant
    GR = TranslationGroup(2) ⋊ SpecialOrthogonalGroup(2)
    gR1 = ArrayPartition([1.0, 0.0], 1 / sqrt(2) * [1.0 1.0; -1.0 1.0])
    gR2 = ArrayPartition([0.0, 1.0], [0.0 -1.0; 1.0 0.0])
    gR3 = ArrayPartition([1.0, 1.0], [1.0 0.0; 0.0 1.0])
    XR1 = ArrayPartition([0.0, 1.0], [0.0 -0.23; 0.23 0.0])
    XR2 = ArrayPartition([1.0, 1.0], [0.0 0.3; -0.3 0.0])
    XR3 = ArrayPartition([1.0, 0.0], [9.0 -0.1; 0.1 0.0])
    propertiesR = Dict(
        :Name => "The special Euclidean group (right variant)",
        :Points => [gR1, gR2, gR3],
        :Vectors => [XR1, XR2, XR3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            # adjoint,
            compose,
            conjugate,
            # diff_inv,
            # diff_left_compose,
            # diff_right_compose,
            exp,
            hat,
            identity_element,
            inv,
            # inv_left_compose,
            # inv_right_compose,
            is_identity,
            # lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectationsR = Dict(
        :repr => "SpecialEuclideanGroup(2; variant=:right)",
        :atol => 1e-14,
        #:diff_inv => -X1,
        #:diff_left_compose => X1,
        #:diff_right_compose => X1,
        #:lie_bracket => zero(X1),
    )
    test_lie_group(GR, propertiesR, expectationsR)
end
