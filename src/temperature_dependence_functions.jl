#=
**Functions of thermal performance curve for model parameters**

We included different functions for temperature dependence :
1) Extended Eppley function
2) Exponential Boltzmann-Arrhenius function
3) Extended Boltzmann-Arrhenius function (Johnson-Lewin)
4) Gaussian (inverted Gaussian) function

In each case, the function returns the biological rate value at a given temperature.
=#

T_param_expba = @NT(norm_constant = 0.2e11, activation_energy = 0.65, β = -0.25)
T_param_extba = @NT(norm_constant = 0.2e11, activation_energy = 0.65, β = -0.25, deactivation_energy = 1.15, T_opt = 293)
T_param_gauss = @NT(norm_constant = 0.5, range = 20, β = -0.25, T_opt = 293)
T_param_eppley = @NT(maxrate_0 = 0.81, eppley_exponent = 0.0631, T_opt = 308.185, range = 35, β = -0.25)

"""
**Option 1 : Extended Eppley function**

This function can be called as an argument in `model_parameters` to define an extended Eppley funtion (ref) for one of:
    - metabolic rate
    - producers growth rate
    - handling time
    - attack rate

| Parameter       | Meaning                                                           | Default values| Reference |
|:----------------|:------------------------------------------------------------------|:--------------|:----------|
| maxrate_0       | Maximum rate at 273.15 degrees Kelvin                             | 0.81          |    TODO   |
| eppley_exponent | Exponential rate of increase                                      | 0.0631        |    TODO   |
| T_opt           | location of the maximum of the quadratic portion of the function  | 308.185       |    TODO   |
| range           | thermal breadth                                                   | 35            |    TODO   |
| β               | allometric exponent                                               | -0.25         |    TODO   |

Example:
TODO
"""

function extended_eppley(T_param)
    return (bodymass, T) -> bodymass^T_param.β * T_param.a * exp.(T_param.b * T) * (1 - ((T - T_param.T_opt) / (T_param.range/2)).^2)
end

"""
**Option 3 : Exponential Boltzmann-Arrhenius function**

| Parameter    | Meaning                               |
|:-------------|:--------------------------------------|
| temp         | temperature range (Kelvin)            |
| p0           | scaling coefficient                   |
| E            | activation energy                     |
| p[:bodymass] | body mass                             |
| beta         | allometric exponent                   |
| k            | Boltzmann constant (k=8.617e-5)       |

Parameters are for instance:

p0=0.2e11
E=0.65
m=1
beta=-0.25

"""

function exponentialBA(T_param)
    k=8.617e-5
    return (bodymass, T) -> T_param.norm_constant.*((bodymass.^T_param.β).*exp.(-T_param.activation_energy./(k*T)))
end

"""
**Option 4 : Extended Boltzmann-Arrhenius function**

| Parameter     | Meaning                                               |
|:--------------|:------------------------------------------------------|
| temp          | temperature range (Kelvin)                            |
| p0            | scaling coefficient                                   |
| E             | activation energy                                     |
| Ed            | deactivation energy                                   |
| topt          | temperature at which trait value is maximal           |
| p[:bodymass]  | body mass                                             |
| beta          | allometric exponent                                   |
| k             | Boltzmann constant (k=8.617e-5)                       |

Parameters are for instance:

p0=0.2e12
E=0.65
Ed=0.72
topt=295
p[:bodymass]=1
beta=-0.25


"""

function extended_BA(T_param)
    kt = 8.617e-5 # Boltzmann constant
    kt = k * T
    pwr = T_param.β*exp(-T_param.activation_energy/kt)
    Δenergy = T_param.deactivation_energy - T_param.activation_energy
    lt = 1 / (1 + exp(-1 / kt * (T_param.deactivation_energy - (T_param.deactivation_energy / T_param.T_opt + k * log(T_param.activation_energy / Δenergy))*T)))
    return(bodymass, T) -> T_param.norm_const * (bodymass^pwr) * lt
end


"""
**Option 5 : Gaussian function**

| Parameter    | Meaning                                        |
|:-------------|:-----------------------------------------------|
| temp         | temperature range (Kelvin)                     |
| p0           | minimal/maximal trait value                    |
| s            | performance breath (width of function)         |
| topt         | temperature at which trait value is maximal    |
| p[:bodymass] | body mass                                      |
| beta         | allometric exponent                            |

Parameters are for instance:

p0=0.5
s=20
topt=295
p[:bodymass]=1
beta=-0.25

"""

function gaussian(T_param)
    if T_param.shape = :hump
        return(bodymass, T) -> bodymass^T_param.β * T_param.norm_const * exp(-(T-T_param.T_opt)^2/(2*T_param.range^2))
    elseif T_param.shape = :U
        return(bodymass, T) -> bodymass^T_param.β * T_param.norm_const * exp((T-T_param.T_opt)^2/(2*T_param.range^2))
    end
end
