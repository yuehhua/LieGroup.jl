using LieGroups, Random, Test, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    G = SpecialEuclideanGroup(2)
    g1 = ArrayPartition(1 / sqrt(2) * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
    g2 = ArrayPartition([0.0 -1.0; 1.0 0.0], [0.0, 1.0])
    g3 = ArrayPartition([1.0 0.0; 0.0 1.0], [1.0, 1.0])
    X1 = ArrayPartition([1.0 0.0; 0.0 0.0], [0.0, 1.0])
    X2 = ArrayPartition([0.0 0.0; 0.0 1.0], [1.0, 1.0])
    X3 = ArrayPartition([9.0 0.5; 0.5 0.0], [1.0, 0.0])
    properties = Dict(
        :Name => "The special Euclidean group",
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
            # exp,
            # hat,
            identity_element,
            inv,
            # inv_left_compose,
            # inv_right_compose,
            # is_identity,
            # lie_bracket,
            # log,
            rand,
            show,
            # vee,
        ],
    )
    expectations = Dict(
        :repr => "SpecialEuclideanGroup(2)",
        #:diff_inv => -X1,
        #:diff_left_compose => X1,
        #:diff_right_compose => X1,
        #:lie_bracket => zero(X1),
    )
    test_lie_group(G, properties, expectations)
end