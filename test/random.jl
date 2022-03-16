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

@testset "Finding candidates for facilitation links." begin

    # Only producer.
    A = [0 0 0; 0 0 0; 0 0 0]
    @test BioEnergeticFoodWebs.wherecanfacilitate(A) == [[2, 1], [3, 1], [1, 2], [3, 2],
        [1, 3], [2, 3]]

    # Producers = species 1 & 2.
    A = [0 0 0; 0 0 0; 1 0 0]
    @test BioEnergeticFoodWebs.wherecanfacilitate(A) == [[2, 1], [3, 1], [1, 2], [3, 2]]

    # Producers = species 1.
    A = [0 0 0; 0 0 1; 1 0 0]
    @test BioEnergeticFoodWebs.wherecanfacilitate(A) == [[2, 1], [3, 1]]

    # No producer.
    A = [0 1 0; 0 0 1; 1 0 0]
    @test BioEnergeticFoodWebs.wherecanfacilitate(A) == []

end

@testset "Creating facilitation matrix from number of links." begin

    # No producer.
    A = [0 1 0; 0 0 1; 1 0 0]
    @test_throws ArgumentError BioEnergeticFoodWebs.facilitationlinks(A, 1)

    # Only producer.
    A = [0 0 0; 0 0 0; 0 0 0]
    @test BioEnergeticFoodWebs.facilitationlinks(A, 0) == [0 0 0; 0 0 0; 0 0 0]
    @test BioEnergeticFoodWebs.facilitationlinks(A, 6) == [0 1 1; 1 0 1; 1 1 0]
    for L in 1:6
        @test count(==(1), BioEnergeticFoodWebs.facilitationlinks(A, L)) == L
    end

end

@testset "Creating facilitation matrix from connectance." begin

    # No producer.
    A = [0 1 0; 0 0 1; 1 0 0]
    @test_throws ArgumentError BioEnergeticFoodWebs.facilitationlinks(A, 1.0)

    # Connectance outside [0,1].
    @test_throws AssertionError BioEnergeticFoodWebs.facilitationlinks(A, 1.1)
    @test_throws AssertionError BioEnergeticFoodWebs.facilitationlinks(A, -0.1)

    # Only producer.
    A = zeros(Int, 4, 4)
    for L in 1:12
        C = L / 16
        @test count(==(1), BioEnergeticFoodWebs.facilitationlinks(A, C)) == L
    end

end

end
