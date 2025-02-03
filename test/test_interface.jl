using LieGroups, Test, ManifoldsBase

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic Lie Group Interface functions" begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    rs = "LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation())"
    @test repr(G) == rs
    𝔤 = LieAlgebra(G)
    op2 = LieGroupsTestSuite.DummySecondOperation()
    rs2 = "LieAlgebra(LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()))"
    @test repr(𝔤) == rs2
    @test is_identity(G, Identity(op))
    @test !is_identity(G, Identity(op2))
    @test base_manifold(G) === M
    e = Identity(op)
    @test compose(G, e, e) == e
    @test compose!(G, e, e, e) === e
    @test isapprox(G, e, Identity(op))
    @test !isapprox(G, e, Identity(op2))
    @test is_point(G, e)
    @test !is_point(G, Identity(op2))
    @test_throws DomainError is_point(G, Identity(op2); error=:error)
    @testset "Methoderrors for the non-implemented mutating variants" begin
        g = :none
        h = :alsonone
        X = :nonetoo
        begin # locally define identity element
            LieGroups.identity_element(::typeof(G)) = :id
            LieGroups.identity_element(::typeof(G), T::Type) = :id
            LieGroups.identity_element!(::typeof(G), g) = (g[] = :id)
            LieGroups.exp!(::typeof(G), h, X) = :id
            @test exp(G, X) === :id
            #
            # same for log
            ManifoldsBase.allocate_result(::typeof(G), ::typeof(LieGroups.log), g) = :g
            LieGroups.log!(::typeof(G), X, g) = :g
            @test log(G, g) === :g
            g2 = Ref(:g)
            inv!(G, g2, e)
            @test g2[] == :id
            # delete methods again
            Base.delete_method(which(identity_element, (typeof(G),)))
            Base.delete_method(which(identity_element, (typeof(G), Type)))
            Base.delete_method(which(identity_element!, typeof.([G, g2])))
            Base.delete_method(which(LieGroups.exp!, typeof.([G, h, X])))
            Base.delete_method(which(ManifoldsBase.allocate_result, typeof.([G, log, g])))
            Base.delete_method(which(LieGroups.log!, typeof.([G, X, g])))
        end
        # so they are undefined here again but we checked the exp fallback
        @test_throws MethodError exp!(G, g, X)
        @test_throws MethodError log!(G, X, g)
    end
end
@testset "Generic Lie Algebra Interface functions" begin
    @testset "Generic get_coordinates/get_vector passthrough on 𝔤" begin
        M = ManifoldsBase.DefaultManifold(2)
        op = AdditionGroupOperation()
        G = LieGroup(M, op)
        𝔤 = LieAlgebra(G)
        B = DefaultOrthonormalBasis()
        p = [1.0, 2.0]
        q = [0.0, 0.0]
        # coordinates and vector on 𝔤 are here the same as the ones on M at 0
        X = [1.0, 0.0]
        @test get_coordinates(𝔤, X, B) == get_coordinates(M, q, X, B)
        Y = copy(X)
        @test get_coordinates!(𝔤, Y, X, B) == get_coordinates!(M, Y, q, X, B)
        @test X == Y
        c = [0.0, 1.0]
        @test get_vector(𝔤, c, B) == get_vector(M, q, c, B)
        d = copy(c)
        @test get_vector!(𝔤, d, c, B) == get_vector!(M, d, q, c, B)
        @test c == d
    end
end
