module TestDefault
  using Base.Test
  using BioEnergeticFoodWebs
  using NamedTuples

  food_chain = [0 1 0 ; 0 0 1; 0 0 0]
  metab_status = [true, false, false]
  p = model_parameters(food_chain, vertebrates = metab_status)

  # DEFAULT ARG.
  # growth rates
  @test p[:r] == [1, 1, 1]
  # metabolic rates
  a_invertebrate = 0.3141
  a_vertebrate = 0.88
  a_producer = 0.138
  @test p[:x] == [a_vertebrate, a_invertebrate, a_producer]
  # maximum consumption rate
  y_vertebrate = 4.0
  y_invertebrate = 8.0
  y_producer = 0.0
  @test p[:y] == [y_vertebrate, y_invertebrate, y_producer]
  # handling time
  @test p[:ht] == 1 ./ p[:y]
  # attack rate
  @test p[:ar] == 1 ./ (0.5 * p[:ht])
  # half saturation constant
  hsc = 1 ./ (p[:ar] .* p[:ht])
  hsc[isnan.(hsc)] .= 0.0
  @test p[:Γ] == hsc

  # PASSED ARG.
  p = model_parameters(food_chain, vertebrates = metab_status,
                       handlingtime = NoEffectTemperature(:handlingtime, parameters_tuple = @NT(y_vertebrate = 3.0, y_invertebrate = 7.0)),
                       attackrate = NoEffectTemperature(:attackrate, parameters_tuple = @NT(Γ = 0.8)),
                       metabolicrate = NoEffectTemperature(:metabolism, parameters_tuple = @NT(a_vertebrate = 0.8, a_invertebrate = 0.3, a_producer = 0.1)),
                       growthrate = NoEffectTemperature(:growth, parameters_tuple = @NT(r = 2.0)))
  @test p[:r] == [2.0, 2.0, 2.0]
  # metabolic rates
  a_invertebrate = 0.3
  a_vertebrate = 0.8
  a_producer = 0.1
  @test p[:x] == [a_vertebrate, a_invertebrate, a_producer]
  # maximum consumption rate
  y_vertebrate = 3.0
  y_invertebrate = 7.0
  y_producer = 0.0
  @test p[:y] == [y_vertebrate, y_invertebrate, y_producer]
  # handling time
  @test p[:ht] == 1 ./ p[:y]
  # attack rate
  @test p[:ar] == 1 ./ (0.8 * p[:ht])
  # half saturation constant
  hsc = 1 ./ (p[:ar] .* p[:ht])
  hsc[isnan.(hsc)] .= 0.0
  @test p[:Γ] == hsc
  @test p[:Γ] == [0.8, 0.8, 0.0]

  # ERRORS

  @test_throws Exception model_parameters(food_chain, vertebrates = metab_status, handlingtime = NoEffectTemperature(:y))
  @test_throws Exception model_parameters(food_chain, vertebrates = metab_status, attackrate = NoEffectTemperature(:y))
  @test_throws Exception model_parameters(food_chain, vertebrates = metab_status, metabolicrate = NoEffectTemperature(:y))
  @test_throws Exception model_parameters(food_chain, vertebrates = metab_status, growthrate = NoEffectTemperature(:y))

end

module TestEppley
  using Base.Test
  using BioEnergeticFoodWebs
  using NamedTuples

  omnivory = [0 1 1 ; 0 0 1 ; 0 0 0]
  metabolic_status = [:true, :false, :false]
  bmass = [100.0, 10.0, 1.0]
  temp = 270.0

  # GROWTH - defaults
  p_g_def = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedEppley(:r))
  r_g_def = bmass .^ -0.25 .* 0.81 .* exp(0.0631 .* (temp.-273.15)) * (1 .- (((temp.-273.15) .- (298.15-273.15)) ./ (35/2)).^2)
  @test p_g_def[:r] == r_g_def

  # GROWTH - changing temperature
  temp_2 = 250.0
  p_g_temp = model_parameters(omnivory, T = temp_2, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedEppley(:r))
  r_g_temp = bmass .^ -0.25 .* 0.81 .* exp(0.0631 .* (temp_2.-273.15)) * (1 .- (((temp_2.-273.15) .- (298.15-273.15)) ./ (35/2)).^2)
  @test p_g_temp[:r] == r_g_temp

  # GROWTH - passed
  pt = @NT(maxrate_0 = 0.6, eppley_exponent = 0.1, T_opt = 215, β = -0.2, range = 20)
  p_g = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedEppley(:r, parameters_tuple = pt))
  r_g = bmass .^ -0.2 .* 0.6 .* exp(0.1 .* (temp.-273.15)) * (1 .- (((temp.-273.15) .- (215-273.15)) ./ (20/2)).^2)
  @test p_g[:r] == r_g

  # METABOLISM - default
  p_x_def = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedEppley(:x))
  @test p_x_def[:x] == p_g_def[:r]

  # METABOLISM - changing temperature
  p_x_temp = model_parameters(omnivory, T = temp_2, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedEppley(:x))
  @test p_x_temp[:x] == p_g_temp[:r]

  # METABOLISM - passed arguments
  pt_x = @NT(maxrate_0_producer = 0.8, maxrate_0_invertebrate = 0.6, maxrate_0_vertebrate = 5,
      eppley_exponent_producer = 0.1, eppley_exponent_invertebrate = 0.06, eppley_exponent_vertebrate = 0.09,
      T_opt_producer = 298, T_opt_invertebrate = 250, T_opt_vertebrate = 270,
      range_producer = 30, range_invertebrate = 32, range_vertebrate = 33,
      β_producer = -0.2, β_invertebrate = -0.3, β_vertebrate = -0.27)
  p_x = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedEppley(:x, parameters_tuple = pt_x))
  maxrate_0_all = [5, 0.6, 0.8]
  eppley_exponent_all = [0.09, 0.06, 0.1]
  T_opt_all = [270, 250, 298]
  T_opt_all = T_opt_all .- 273.15
  range_all = [33, 32, 30]
  β_all = [-0.27, -0.3, -0.2]
  x_x = bmass.^β_all .* maxrate_0_all .* exp.(eppley_exponent_all .* (temp.-273.15)) .* (1 .- (((temp.-273.15) .- T_opt_all) ./ (range_all./2)).^2)
  @test x_x == p_x[:x]

  # ERRORS
  @test_throws Exception model_parameters(omnivory, metabolicrate = ExtendedEppley(:handlingtime))
  @test_throws Exception model_parameters(omnivory, metabolicrate = ExtendedEppley(:attackrate))
  @test_throws Exception model_parameters(omnivory, metabolicrate = ExtendedEppley(:y))

end

module TestExponentialBA
  using Base.Test
  using BioEnergeticFoodWebs
  using NamedTuples

  omnivory = [0 1 1 ; 0 0 1 ; 0 0 0]
  metabolic_status = [:true, :false, :false]
  bmass = [100.0, 10.0, 1.0]
  temp = 270.0

  #GROWTH
  #defaults
  p_r_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExponentialBA(:r))
  k = 8.617e-5
  T0K = 273.15
  r_d = exp(-16.54) .* (bmass .^-0.31) .* exp.(-0.69 .* (293.15 .- (temp + T0K)) ./ (k * (temp + T0K) .* 293.15))
  @test p_r_d[:r] == r_d
  #change temperature
  temp2 = 250.0
  p_r_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExponentialBA(:r))
  r_t = exp(-16.54) .* (bmass .^-0.31) .* exp.(-0.69 .* (293.15 .- (temp2 + T0K)) ./ (k * (temp2 + T0K) .* 293.15))
  @test p_r_t[:r] == r_t
  #passed arguments
  pt_r = @NT(norm_constant = -18, activation_energy = -0.8, T0 = 290, β = -0.25)
  p_r_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExponentialBA(:r, parameters_tuple = pt_r))
  r_2 = exp(-18) .* (bmass .^-0.25) .* exp.(-0.8 .* (290 .- (temp + T0K)) ./ (k * (temp + T0K) .* 290))
  @test p_r_2[:r] == r_2

  #METABOLISM
  #defaults
  p_x_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExponentialBA(:x))
  @test p_x_d[:x] == r_d
  #change temperature
  p_x_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExponentialBA(:x))
  @test p_x_t[:x] == r_t
  #passed arguments
  pt_x = @NT(norm_constant_producer = -16, norm_constant_invertebrate = -17, norm_constant_vertebrate = -18,
             activation_energy_producer = -0.6, activation_energy_invertebrate = -0.7, activation_energy_vertebrate = -0.8,
             T0_producer = 270, T0_invertebrate = 280, T0_vertebrate = 290,
             β_producer = -0.2, β_invertebrate = -0.3, β_vertebrate = -0.4)
  norm_constant_all = [-18, -17, -16]
  activation_energy_all = [-0.8, -0.7, -0.6]
  T0_all = [290, 280, 270]
  β_all = [-0.4, -0.3, -0.2]
  p_x_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExponentialBA(:x, parameters_tuple = pt_x))
  x_2 = exp.(norm_constant_all) .* (bmass .^β_all) .* exp.(activation_energy_all .* (T0_all .- (temp + T0K)) ./ (k * (temp + T0K) .* T0_all))
  @test p_x_2[:x] == x_2

  #ATTACK
  #defaults
  p_ar_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExponentialBA(:attackrate))
  ar_d = exp.([-13.1, -13.1, 0.0]) .* (bmass .^[-0.8, -0.8, 0.0]) .* (bmass' .^[0 -0.8 0.25; 0 0 0.25 ; 0 0 0]) .* exp.([-0.38, -0.38, 0.0] .* ([293.15, 293.15, 0.0] .- (temp + T0K)) ./ (k * (temp + T0K) .* [293.15, 293.15, 0.0]))
  ar_d[isnan.(ar_d)] = 0
  @test p_ar_d[:ar] == ar_d
  #change temperature
  p_ar_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExponentialBA(:attackrate))
  ar_t = exp.([-13.1, -13.1, 0.0]) .* (bmass .^[-0.8, -0.8, 0.0]) .* (bmass' .^[0 -0.8 0.25; 0 0 0.25 ; 0 0 0]) .* exp.([-0.38, -0.38, 0.0] .* ([293.15, 293.15, 0.0] .- (temp2 + T0K)) ./ (k * (temp2 + T0K) .* [293.15, 293.15, 0.0]))
  ar_t[isnan.(ar_t)] = 0
  @test p_ar_t[:ar] == ar_t
  #passed arguments
  pt_ar = @NT(norm_constant_vertebrate = -12, norm_constant_invertebrate = -14,
  						activation_energy_vertebrate = -0.3, activation_energy_invertebrate = -0.4,
  						T0_vertebrate = 290, T0_invertebrate = 270,
  						β_producer = 0.2, β_vertebrate = -0.9, β_invertebrate = -0.7)
  p_ar_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExponentialBA(:attackrate, parameters_tuple = pt_ar))
  ar_2 = exp.([-12, -14, 0.0]) .* (bmass .^[-0.9, -0.7, 0.0]) .* (bmass' .^[0 -0.7 0.2; 0 0 0.2 ; 0 0 0]) .* exp.([-0.3, -0.4, 0.0] .* ([290, 270, 0.0] .- (temp + T0K)) ./ (k * (temp + T0K) .* [290, 270, 0.0]))
  ar_2[isnan.(ar_2)] = 0
  @test p_ar_2[:ar] == ar_2

  #HANDLING
  #defaults
  p_ht_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, handlingtime = ExponentialBA(:handlingtime))
  ht_d = exp.([9.66, 9.66, 0.0]) .* (bmass .^[0.47, 0.47, 0.0]) .* (bmass' .^[0.0 0.47 -0.45 ; 0 0 -0.45; 0 0 0]) .* exp.([0.26, 0.26, 0.0] .* ([293.15, 293.15, 0.0] .- (temp + T0K)) ./ (k * (temp + T0K) .* [293.15, 293.15, 0.0]))
  ht_d[isnan.(ht_d)] .= 0
  @test p_ht_d[:ht] == ht_d
  #change temperature
  p_ht_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, handlingtime = ExponentialBA(:handlingtime))
  ht_t = exp.([9.66, 9.66, 0.0]) .* (bmass .^[0.47, 0.47, 0.0]) .* (bmass' .^[0.0 0.47 -0.45 ; 0 0 -0.45; 0 0 0]) .* exp.([0.26, 0.26, 0.0] .* ([293.15, 293.15, 0.0] .- (temp2 + T0K)) ./ (k * (temp2 + T0K) .* [293.15, 293.15, 0.0]))
  ht_t[isnan.(ht_t)] = 0
  @test p_ht_t[:ht] == ht_t
  #passed arguments
  pt_ht = @NT(norm_constant_vertebrate = 9, norm_constant_invertebrate = 10,
  						activation_energy_vertebrate = 0.2, activation_energy_invertebrate = 0.3,
  						T0_vertebrate = 290, T0_invertebrate = 270,
  						β_producer = -0.4, β_vertebrate = 0.3, β_invertebrate = 0.5)
  p_ht_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, handlingtime = ExponentialBA(:handlingtime, parameters_tuple = pt_ht))
  ht_2 = exp.([9, 10, 0.0]) .* (bmass .^[0.3, 0.5, 0.0]) .* (bmass' .^[0.0 0.5 -0.4 ; 0 0 -0.4; 0 0 0]) .* exp.([0.2, 0.3, 0.0] .* ([290, 270, 0.0] .- (temp + T0K)) ./ (k * (temp + T0K) .* [290, 270, 0.0]))
  ht_2[isnan.(ht_2)] = 0
  @test p_ht_2[:ht] == ht_2

  #ERRORS
  @test_throws Exception model_parameters(omnivory, metabolicrate = ExponentialBA(:y))

end

module TestExtendedBA
  using Base.Test
  using BioEnergeticFoodWebs
  using NamedTuples

  omnivory = [0 1 1 ; 0 0 1 ; 0 0 0]
  metabolic_status = [:true, :false, :false]
  bmass = [100.0, 10.0, 1.0]
  temp = 270.0
  temp2 = 250.0
  k = 8.617e-5 # Boltzmann constant

  #GROWTH
  #defaults
  p_r_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedBA(:r))
  Δenergy = 1.15 .- 0.53
  r_d = 3e8 .* bmass .^(-0.25) .* exp.(.-0.53 ./ (k * temp)) .* (1 ./ (1 + exp.(-1 / (k * temp) .* (1.15 .- (1.15 ./ 298.15 .+ k .* log(0.53 ./ Δenergy)).*temp))))
  @test p_r_d[:r] == r_d
  #change temperature
  p_r_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedBA(:r))
  @test p_r_d[:r] != p_r_t[:r]
  #passed arguments
  pt_r = @NT(norm_constant = 4000)
  p_r_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = ExtendedBA(:r, parameters_tuple = pt_r))
  r_2 = (r_d ./ 3e8) .* 4000
  @test p_r_2[:r] ≈ r_2 atol=1e-6

  #METABOLISM
  #defaults
  p_x_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedBA(:x))
  x_d = 3e8 .* bmass .^(-0.25) .* exp.(.-0.53 ./ (k * temp)) .* (1 ./ (1 + exp.(-1 / (k * temp) .* (1.15 .- (1.15 ./ 298.15 .+ k .* log.(0.53 ./ Δenergy)).* temp))))
  @test p_x_d[:x] == x_d
  #change temperature
  p_x_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedBA(:x))
  @test p_x_t[:x] != p_x_d[:x]
  #passed arguments
  pt_x = @NT(β_vertebrate = -0.5)
  p_x_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = ExtendedBA(:x, parameters_tuple = pt_x))
  p_x_2[:x] != p_x_d[:x]
  x_2 = 3e8 .* bmass .^([-0.5, -0.25, -0.25]) .* exp.(.-0.53 ./ (k * temp)) .* (1 ./ (1 + exp.(-1 / (k * temp) .* (1.15 .- (1.15 ./ 298.15 .+ k .* log.(0.53 ./ Δenergy)).* temp))))
  @test p_x_2[:x] == x_2

  #ATTACK
  #defaults
  p_ar_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExtendedBA(:attackrate))
  ar_d = [3e8, 3e8, 3e8] .* bmass .^([-0.25, -0.25, 0.0]) .* bmass' .^([0.0 -0.25 -0.25 ; 0.0 0.0 -0.25 ; 0.0 0.0 0.0]) .* exp.(.-[0.53, 0.53, 0.53] ./ (k * temp)) .* (1 ./ (1 .+ exp.(-1 / (k * T) .* ([1.15, 1.15, 1.15] .- ([1.15, 1.15, 1.15] ./ [298.15, 298.15, 0.0] .+ k .* log.([0.53, 0.53, 0.53] ./ Δenergy)) .* temp))))
  ar_d[isnan.(ar_d)] = 0
  @test p_ar_d[:ar] == ar_d
  #change temperature
  p_ar_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExtendedBA(:attackrate))
  @test p_ar_d[:ar] != p_ar_t[:ar]
  #passed arguments
  pt_ar = @NT(norm_constant_invertebrate = 3e7)
  p_ar_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExtendedBA(:attackrate, parameters_tuple = pt_ar))
  @test p_ar_2[:ar] != p_ar_d[:ar]
  ar_2 = [3e8, 3e7, 3e8] .* bmass .^([-0.25, -0.25, 0.0]) .* bmass' .^([0.0 -0.25 -0.25 ; 0.0 0.0 -0.25 ; 0.0 0.0 0.0]) .* exp.(.-[0.53, 0.53, 0.53] ./ (k * temp)) .* (1 ./ (1 .+ exp.(-1 / (k * T) .* ([1.15, 1.15, 1.15] .- ([1.15, 1.15, 1.15] ./ [298.15, 298.15, 0.0] .+ k .* log.([0.53, 0.53, 0.53] ./ Δenergy)) .* temp))))
  @test p_ar_2[:ar] == ar_2

  #ERRORS
  @test_throws Exception model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExtendedBA(:handlingtime))
  @test_throws Exception model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = ExtendedBA(:y))

end

module TestGaussian
  using Base.Test
  using BioEnergeticFoodWebs
  using NamedTuples

  omnivory = [0 1 1 ; 0 0 1 ; 0 0 0]
  metabolic_status = [:true, :false, :false]
  bmass = [100.0, 10.0, 1.0]
  temp = 270.0
  temp2 = 250.0

  #GROWTH
  #defaults
  p_r_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = Gaussian(:r))
  r_d = bmass .^ -0.25 .* 0.5 .* exp( .- (temp .- 298.15) .^ 2 ./ (2 .* 20 .^ 2))
  @test p_r_d[:r] == r_d
  #change temperature
  p_r_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, growthrate = Gaussian(:r))
  @test p_r_t[:r] != p_r_d[:r]
  #passed arguments
  pt_r = @NT(norm_constant = 0.7)
  p_r_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, growthrate = Gaussian(:r, parameters_tuple = pt_r))
  @test p_r_2[:r] != p_r_d[:r]
  r_2 = bmass .^ -0.25 .* 0.7 .* exp( .- (temp .- 298.15) .^ 2 ./ (2 .* 20 .^ 2))
  @test p_r_2[:r] == r_2

  #METABOLISM
  #defaults
  p_x_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = Gaussian(:x))
  x_d = bmass .^ -0.25 .* 0.5 .* exp.( .- (temp .- 298.15) .^ 2 ./ (2 .* 20 .^ 2))
  @test p_x_d[:x] == x_d
  #change temperature
  p_x_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = Gaussian(:x))
  @test p_x_t[:x] != p_x_d[:x]
  #passed arguments
  pt_x = @NT(range_producer = 30)
  p_x_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, metabolicrate = Gaussian(:x, parameters_tuple = pt_x))
  @test p_x_2[:x] != p_x_d[:x]
  x_2 = bmass .^ -0.25 .* 0.5 .* exp.( .- (temp .- 298.15) .^ 2 ./ (2 .* [20,20,30] .^ 2))
  @test p_x_2[:x] == x_2

  #ATTACK
  #defaults
  p_ar_d = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = Gaussian(:attackrate))
  ar_d = bmass .^ [-0.25,-0.25,0.0] .* bmass' .^ [0.0 -0.25 -0.25 ; 0.0 0.0 -0.25 ; 0.0 0.0 0.0] .* [0.5, 0.5, 0.0] .* exp.( .- (temp .- [295,295,0]) .^ 2 ./ (2 .* [20,20,0] .^ 2))
  ar_d[isnan.(ar_d)] .= 0
  @test p_ar_d[:ar] == ar_d
  #change temperature
  p_ar_t = model_parameters(omnivory, T = temp2, bodymass = bmass, vertebrates = metabolic_status, attackrate = Gaussian(:attackrate))
  @test p_ar_t[:ar] != p_ar_d[:ar]
  #passed arguments
  pt_ar = @NT(T_opt_invertebrate = 270)
  p_ar_2 = model_parameters(omnivory, T = temp, bodymass = bmass, vertebrates = metabolic_status, attackrate = Gaussian(:attackrate, parameters_tuple = pt_ar))
  @test p_ar_2[:ar] != p_ar_d[:ar]
  ar_2 = bmass .^ [-0.25,-0.25,0.0] .* bmass' .^ [0.0 -0.25 -0.25 ; 0.0 0.0 -0.25 ; 0.0 0.0 0.0] .* [0.5, 0.5, 0.0] .* exp.( .- (temp .- [295,270,0]) .^ 2 ./ (2 .* [20,20,0] .^ 2))
  @test p_ar_2[:ar] == ar_2

  #HANDLING
  #defaults
  #change temperature
  #passed arguments
  #ERRORS

end