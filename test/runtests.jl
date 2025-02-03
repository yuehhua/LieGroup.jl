#!/usr/bin/env julia --optimize=0
#
#

using LieGroups, Test

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

function include_test(path)
    @info "Testing $path"
    @time include(path)  # show basic timing, (this prints a newline at end)
end

@testset "LieGroups.jl" begin
    @testset "Lie Group Interface" begin
        include_test("test_interface.jl")
    end
    @testset "Generic Group Operations" begin
        include_test("operations/test_addition_operation.jl")
        include_test("operations/test_multiplication_operation.jl")
    end
    @testset "Generic Group Actions" begin
        include_test("actions/test_action_interface.jl")
        include_test("actions/test_operation_action.jl")
    end
    @testset "Lie Groups" begin
        @testset "Meta Lie Groups" begin
            include_test("groups/test_power_group.jl")
            include_test("groups/test_product_group.jl")
            include_test("groups/test_semidirect_product_group.jl")
        end
        include_test("groups/test_general_linear_group.jl")
        include_test("groups/test_heisenberg_group.jl")
        include_test("groups/test_orthogonal_group.jl")
        include_test("groups/test_unitary_group.jl")
        include_test("groups/test_special_orthogonal_group.jl")
        include_test("groups/test_special_euclidean_group.jl")
        include_test("groups/test_special_unitary_group.jl")
        include_test("groups/test_translation_group.jl")
    end
    include("test_aqua.jl")
end
