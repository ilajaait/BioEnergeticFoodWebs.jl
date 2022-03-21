module TestRandom
using Test
using BioEnergeticFoodWebs

@testset "Niche model." begin
    @test_throws AssertionError nichemodel(10, 150)

    A = nichemodel(10, 20)
    @test size(A, 1) == 10
    @test size(A, 2) == 10

    A = nichemodel(10, 0.12)
    @test size(A, 1) == 10
    @test size(A, 2) == 10

    A = nichemodel(10, 0.12, toltype=:rel)
    A = nichemodel(10, 0.12, toltype=:abs)

    A = [0 1; 1 0]
    @test BioEnergeticFoodWebs.connectance(A) == 0.5
end

@testset "Compute number of species." begin
    A = zeros(Int, 3, 3)
    @test BioEnergeticFoodWebs.numberspecies(A) == 3
    A = zeros(Int, 10, 10)
    @test BioEnergeticFoodWebs.numberspecies(A) == 10
end

@testset "Identifying producers correctly." begin
    A = [0 0; 0 0]
    @test BioEnergeticFoodWebs.whoisproducer(A) == [1, 1]
    A = [0 1; 0 0]
    @test BioEnergeticFoodWebs.whoisproducer(A) == [0, 1]
    A = [0 1; 1 0]
    @test BioEnergeticFoodWebs.whoisproducer(A) == [0, 0]
end

end
