module TestSimulate
  using Base.Test
  using BioEnergeticFoodWebs

  food_chain = [0 0 0;0 0 1; 1 0 0]
  p = model_parameters(food_chain)
  # Fail if the biomass vector has bad length
  @test_throws AssertionError simulate(p, rand(5))
  # Fail if stop < start
  @test_throws AssertionError simulate(p, rand(3), start=0, stop=-1)
end

module TestSimulateHandChecked
  using Base.Test
  using BioEnergeticFoodWebs

  A = [0 1 0; 0 0 0; 0 1 0]
  p = model_parameters(A)
  b0 = vec([0.2 0.4 0.1])
  der = similar(b0)
  BioEnergeticFoodWebs.dBdt(der, b0, p, 0.0)
  # 0.1604888888888889,-0.4,0.08024444444444445
  @test der[1] ≈ 0.160 atol = 0.01
  @test der[2] ≈ -0.4 atol = 0.01
  @test der[3] ≈ 0.080 atol = 0.01
end

module TestSimulateSanityCheck
  using Base.Test
  using BioEnergeticFoodWebs

  # A producer with no predation reaches K
  free_producers = [0 1 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0]
  p = free_producers |> model_parameters
  n = rand(4)

  s = simulate(p, n, start=0, stop=25)
  @test s[:B][26,4] ≈ p[:K] atol=0.002

  # Using system-wide regulation, producers with no consumption reach K / n
  A = zeros(Int64, (4, 4))
  p = model_parameters(A, productivity=:system)
  n = rand(4)
  s = simulate(p, n, start=0, stop=15)
  @test s[:B][16,4] ≈ p[:K]/4 atol=0.001

  # A consumer with a resource with 0 biomass goes extinct
  A = [0 1; 0 0]
  p = A |> model_parameters
  n = vec([1.0, 0.0])

  # s = simulate(p, n, start=0, stop=50, use=:stiff)
  # @test_approx_eq_eps s[:B][end,1] 0.0 0.001
  # @test_approx_eq_eps s[:B][end,2] 0.0 0.001

  s = simulate(p, n, start=0, stop=50)
  @test s[:B][end,1] ≈ .0 atol=0.001
  @test s[:B][end,2] ≈ 0.0 atol=0.001

  #test that stiff and nonstiff methods yield similar results
  A = [0 0 1 ; 0 0 1 ; 0 0 0]
  p = model_parameters(A)
  b = [.5, .5, .5]
  s_stiff = simulate(p, b, use = :stiff)
  meanB_stiff = mean(s_stiff[:B], 1)
  s_nonstiff = simulate(p, b, use = :nonstiff)
  meanB_nonstiff = mean(s_stiff[:B], 1)
  @test meanB_stiff == meanB_nonstiff
end


module TestSimulateProductivity
  using Base.Test
  using BioEnergeticFoodWebs

  A = zeros(Int64, (4, 4))
  n = rand(4)

  # Using system-wide regulation, producers with no consumption reach K / n
  p = model_parameters(A, productivity=:system)
  s = simulate(p, n, start=0, stop=15)
  @test s[:B][16,4] ≈ p[:K]/4 atol=0.001

  # Using species-wide regulation, producers with no consumption reach K
  p = model_parameters(A, productivity=:species)
  s = simulate(p, n, start=0, stop=15)
  @test s[:B][16,4] ≈ p[:K] atol=0.001

  # NOTE The following tests start with n = [1 1 1 1]

  # Using competitive regulation with α = 1 is neutral
  p = model_parameters(A, productivity=:competitive, α=1.0)
  s = simulate(p, ones(4), start=0, stop=15)
  @test s[:B][16,4] ≈ p[:K]/4 atol=0.001

  # Using competitive regulation with α > 1 is exclusive
  p = model_parameters(A, productivity=:competitive, α=1.05)
  s = simulate(p, ones(4), start=0, stop=15)
  @test s[:B][16,4] < p[:K]/4

  # Using competitive regulation with α < 1 is overyielding
  p = model_parameters(A, productivity=:competitive, α=0.95)
  s = simulate(p, ones(4), start=0, stop=15)
  @test s[:B][16,4] > p[:K]/4

  # Using competitive regulation with α = 0 is species-level regulation
  p = model_parameters(A, productivity=:competitive, α=0.0)
  s = simulate(p, ones(4), start=0, stop=15)
  @test s[:B][16,4] ≈ p[:K] atol=0.001

  # Using nutrient intake
  p = model_parameters(A, productivity=:nutrients)
  G = BioEnergeticFoodWebs.growthrate(p, ones(4), 1 ; c = p[:supply])
  g = BioEnergeticFoodWebs.get_growth(ones(4), p ; c = p[:supply])[:growth]
  @test G[1] ≈ 0.964 atol=0.001
  @test g[1] == g[2] == g[3] == g[4]
end

module TestSimulateRewiring
  using Base.Test
  using BioEnergeticFoodWebs


  #print info message when extinction occurs
  A = [0 0 0; 0 1 0; 0 0 0]
  A_after = Int.(zeros(A))
  b = rand(3)

  # p = model_parameters(A, rewire_method = :Gilljam)
  # @test_warn "INFO: extinction of species 2" simulate(p, b)
  # @test p[:A] == A_after
  # p = model_parameters(A, rewire_method = :stan)
  # @test_warn "INFO: extinction of species 2" simulate(p, b)
  # @test_warn "INFO: extinction of species 2" simulate(p, b)
  # @test p[:A] == A_after
  # p = model_parameters(A, rewire_method = :ADBM)
  # @test_warn "INFO: extinction of species 2" simulate(p, b)
  # @test_warn "INFO: extinction of species 2" simulate(p, b)
  # @test p[:A] == A_after

  #don't print info message when rewire_method = :none, even when an extinction occurs
  p = model_parameters(A)
  @test_nowarn simulate(p, b)
  #don't print info message when rewire_method != :none but no extinction occurs
  A = [0 0 0 ; 0 0 0 ; 0 0 0] #free producers
  p = model_parameters(A, rewire_method = :stan)
  b = rand(3)
  @test_nowarn simulate(p, b)

end

module TestSimulateNP
  using Base.Test
  using BioEnergeticFoodWebs

  A = [0 1 1 ; 0 0 0 ; 0 0 0]
  b0 = [0.5, 0.5, 0.5]
  k1 = [0, 0.1, 0.2]
  k2 = [0, 0.15, 0.15]
  # When nutrient concentration is 0, then producers growth is 0 and nutrient growth is 1
  c0 = [0.0, 0.0]
  p = model_parameters(A, productivity = :nutrients, K1 = k1, K2 = k2)
  G = zeros(3)
  for i in 1:3
    if p[:is_producer][i]
      G[i] = BioEnergeticFoodWebs.growthrate(p, b0, i, c = c0)[1]
    else
      G[i] = 0.0
    end
  end
  n = BioEnergeticFoodWebs.nutrientuptake(c0, b0, p, G)
  @test n ≈ [1.0, 1.0] atol=0.001
  @test G ≈ [0.0, 0.0, 0.0] atol=0.001

  c0 = [0.0, 1.0]
  G = zeros(3)
  for i in 1:3
    if p[:is_producer][i]
      G[i] = BioEnergeticFoodWebs.growthrate(p, b0, i, c = c0)[1]
    else
      G[i] = 0.0
    end
  end
  n = BioEnergeticFoodWebs.nutrientuptake(c0, b0, p, G)
  @test n[1] ≈ 1.0 atol=0.001
  @test G ≈ [0.0, 0.0, 0.0] atol=0.001

  #sanity check
  c0 = [3.0, 3.0]
  G = zeros(3)
  for i in 1:3
    if p[:is_producer][i]
      G[i] = BioEnergeticFoodWebs.growthrate(p, b0, i, c = c0)[1]
    else
      G[i] = 0.0
    end
  end
  n = BioEnergeticFoodWebs.nutrientuptake(c0, b0, p, G)
  @test G[1] ≈ 0.0 atol=0.001
  @test G[2] ≈ 0.9524 atol=0.001
  @test G[3] ≈ 0.9375 atol=0.001
  @test n[1] ≈ -0.69494 atol=0.001
  @test n[2] ≈ -0.2225 atol=0.001

end
