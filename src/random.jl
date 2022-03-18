"""
** Connectance of a network**

Returns the connectance of a square matrix, defined as ``S/L^2``.

"""
function connectance(S::Int64, L::Int64)
    @assert S > 1
    C = L / S^2
    @assert C < 1.0
    @assert C > 0.0
    return C
end

function connectance(A::Array{Int64,2})
    return connectance(size(A, 1), sum(A))
end

"""
**Niche model of food webs**

Takes a number of species `S` and a number of interactions `L`, and returns
a food web with predators in rows, and preys in columns. This function is
used internally by `nichemodel` called with a connectance.

"""
function nichemodel(S::Int64, L::Int64)

    C = connectance(S, L)

    # Beta distribution parameter
    β = 1 / (2 * C) - 1.0

    # Pre-allocate the network
    A = zeros(Int64, (S, S))

    # Generate body size
    n = sort(rand(Uniform(0.0, 1.0), S))

    # Pre-allocate centroids
    c = zeros(Float64, S)

    # Generate random ranges
    r = n .* rand(Beta(1.0, β), S)

    # Generate random centroids
    for s in 1:S
        c[s] = rand(Uniform(r[s] / 2, n[s]))
    end

    # The smallest species has a body size and range of 0
    n[n.==minimum(n)] .= 0.0
    r[n.==minimum(n)] .= 0.0

    for consumer in 1:S
        for resource in 1:S
            if n[resource] < c[consumer] + r[consumer]
                if n[resource] > c[consumer] - r[consumer]
                    A[consumer, resource] = 1
                end
            end
        end
    end

    return A

end

"""
**Niche model of food webs**

Takes a number of species `S` and a connectance `C`, and returns a food web
with predators in rows, and preys in columns. Note that the connectance is
first transformed into an integer number of interactions.

This function has two keyword arguments:

1. `tolerance` is the allowed error on tolerance (see below)

2. `toltype` is the type or error, and can be `:abs` (absolute) and `:rel`
(relative). Relative tolerance is the amount of error allowed, relative to the
desired connectance value. If the simulated network has a tolerance x, the
target connectance is c, then the relative error is |1-x/c|.

"""
function nichemodel(S::Int64, C::Float64; tolerance::Float64=0.05, toltype::Symbol=:abs)
    @assert C < 1.0
    @assert C > 0.0
    @assert tolerance > 0.0
    @assert toltype ∈ [:abs, :rel]
    L = round(Int64, C * S^2)
    A = nichemodel(S, L)
    if toltype == :abs
        tolfunc = (x) -> abs(x - C) < tolerance
    else
        tolfunc = (x) -> abs(1 - x / C) < tolerance
    end
    while !(tolfunc(connectance(A)))
        A = nichemodel(S, L)
    end
    return A
end

"""Compute number of species in the system from the trophic matrix."""
numberspecies(A) = size(A, 1)

"""
**Who is a primary producer?**

Return BitVector of primary producers.
1 - The species is a primary producer.
0 - The species is a consumer (i.e. consumes at least one other species).
"""
function whoisproducer(A)
    vec(.!any(A, dims=2))
end

"""
**Create facilitation links**

Takes trophic matrix (A) and number or connectance of facilitation links (resp. L and C).
Returns matrix of facilitation links.
Facilitating species are in rows, facilitated species are in columns.
"""
function facilitationlinks(A, L::Int64)

    # Test.
    @assert L >= 0

    # Find candidate links.
    S = numberspecies(A)
    facilitation_matrix = zeros(Int, S, S) # initialize
    candidates = wherecanfacilitate(A)
    size(candidates, 1) < L ? throw(ArgumentError("L is too large.")) : nothing

    # Choose randomly between candidates.
    choosen = sample(candidates, L, replace=false)
    for l in 1:L
        i, j = choosen[l]
        facilitation_matrix[i, j] = 1
    end

    facilitation_matrix
end

function facilitationlinks(A, C::Float64)

    # Test. Connectance is in [0,1]?
    @assert C >= 0.0
    @assert C <= 1.0

    # Create facilitation matrix.
    S = numberspecies(A)
    L = C * (S^2) # number of facilitation links
    L = Int64(round(L)) # formatting
    #? Creal = L / S^2 # realized connectance
    facilitationlinks(A, L) # facilitation matrix
end

"""
**Where facilitation links can be added?**

Find directed pair of species between which facilitation could be added.
Returns the list of indexes of candidate pairs.
"""
function wherecanfacilitate(A)

    # Set up.
    S = numberspecies(A)
    candidates = []
    isproducer = whoisproducer(A)
    producers = (1:S)[isproducer] # indexes of producers

    # Loop over facilitation matrix and add candidates.
    for j in producers
        for i in 1:S
            i != j ? push!(candidates, [i, j]) : nothing
        end
    end

    candidates
end

""" Find sessile species. """
function whoissessile(A)
    issessile = whoisproducer(A) #? sessile ⟺ producer (can change)
    issessile
end

"""
**Where competition for space links can be added?**

Find directed pair of species between which facilitation could be added.
Returns the list of indexes of candidate pairs.
"""
function wherecancompete(A)

    # Set up.
    S = numberspecies(A)
    candidates = []
    issessile = whoissessile(A)
    sessile = (1:S)[issessile] # indexes sessile species

    # Loop to find candidates for competition.
    for i in sessile
        for j in sessile
            i != j ? push!(candidates, [i, j]) : nothing
        end
    end

    candidates
end

"""
**Create competition links**

Takes trophic matrix (A) and number or connectance of competition links (resp. L and C).
Returns matrix of competition links.
Facilitating species are in rows, facilitated species are in columns.
"""
function competitionlinks(A, L::Int64)

    # Test.
    @assert L >= 0

    # Set up.
    S = numberspecies(A)
    competition_matrix = zeros(Int, S, S) # initialize
    candidates = wherecancompete(A)
    size(candidates, 1) < L ? throw(ArgumentError("Too many links. Decrease L.")) : nothing

    # Choose randomly between candidates.
    choosen = sample(candidates, L, replace=false)
    for l in 1:L
        i, j = choosen[l] # pick index of selected link
        competition_matrix[i, j] = 1 # create link
    end

    competition_matrix
end

function competitionlinks(A, C::Float64)
    # Test. Connectance is in [0,1]?
    @assert C >= 0.0
    @assert C <= 1.0

    # Create competition matrix.
    S = numberspecies(A)
    L = C * (S^2) # number of competition links
    L = Int64(round(L)) # formatting
    #? Creal = L / S^2 # realized connectance
    competitionlinks(A, L) # competition matrix
end

# TODO: Create predator interference matrix

""" Find predators of a given prey, given the trophic network (A). """
function whoispredator(A, prey)
    S = numberspecies(A)
    ispredator = A[:, prey] .== 1
    (1:S)[ispredator]
end

function wherecaninterfere(A)

    #  Set up.
    S = numberspecies(A)
    candidates = []

    # Loop over species to find candidates.
    for i in 1:S
        predators = whoispredator(A, i)
        length(predators) >= 2 ? append!(candidates, push!.(allpairs(predators), i)) : nothing
    end

    candidates
end

""" Find all pairs of values of a given vector """
function allpairs(vector)
    n = length(vector)
    pairs = []
    for i in 1:n, j in (i+1):n
        push!(pairs, [vector[i], vector[j]])
        push!(pairs, [vector[j], vector[i]])
    end
    pairs
end

# TODO: Create refuge matrix
