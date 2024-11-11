using Dates

struct System
    L::Float64
    μ::Float64
    β::Float64
    Vext::Function
    ϕ::Function
    particles::Vector{Float64}
    amount::Vector{Float64}
    System(L::Number, μ::Number, T::Number, Vext::Function, ϕ::Function) = new(L, μ, 1 / T, Vext, ϕ, [], [])
end

mutable struct Histograms
    bins::Int
    dx::Float64
    ρ::Vector{Float64}
    count::Int
    function Histograms(system::System; bins=1000)
        new(bins, system.L / bins, zeros(bins), 0)
    end
end

mutable struct gHistograms
    bins::Int
    dx::Float64
    g::Vector{Float64}
    count::Int
    function gHistograms(system::System; bins=500)
        new(bins, system.L / (2*bins), zeros(bins), 0)
    end
end

function pbc!(system::System, i)
    system.particles[i] -= floor(system.particles[i] / system.L) * system.L
end

function dist(xi, xj, L)
    result = xj - xi
    result -= round(result / L) * L
    abs(result)
end

function add_particle!(system::System, x)
    push!(system.particles, x)
end

function remove_particle!(system::System, i)
    deleteat!(system.particles, i)
end

function calc_energy(system::System, i)
    xi = system.particles[i] #xi=particle position
    E = system.Vext(xi)
    for xj in system.particles
        if xi == xj
            continue
        end
        E += system.ϕ(dist(xi, xj, system.L)) # ϕ only dependent from particle distance
    end
    E
end

function trial_insert(system::System)
    add_particle!(system, rand() * system.L)
    i = length(system.particles)
    ΔE = calc_energy(system, i)
    if rand() > system.L / length(system.particles) * exp(system.β * (system.μ - ΔE))
        # Rejected, undo trial insert
        remove_particle!(system, i)
    end
end

function trial_delete(system::System)
    if isempty(system.particles)
        return
    end
    i = rand(1:length(system.particles))
    ΔE = calc_energy(system, i)
    if rand() < length(system.particles) / system.L * exp(system.β * (ΔE - system.μ))
        # Accepted, do the actual removal
        remove_particle!(system, i)
    end
end

function get_results(system::System, histograms::Histograms, ghist::gHistograms)
    dx = histograms.dx
    xs = collect(dx/2:dx:system.L-dx/2)
    leng = 500
    ρ = histograms.ρ / (histograms.count * dx)
    avN = sum(system.amount) / length(system.amount)
    ρb =  avN / system.L
    xs, ρ, xs[1:leng], 2*ghist.g / (avN* ghist.count * dx * ρb), ρb   
end

function trial_move(system::System; Δxmax=0.1)
    if isempty(system.particles)  # If no particles are in the system, do nothing
        return
    end
    i = rand(1:length(system.particles))  # Select random particle with index i, collection in brackets
    xbefore = system.particles[i]  # Save its initial position
    Ebefore = calc_energy(system, i)  # Calculate the initial potential energy of particle i
    system.particles[i] += Δxmax * (2 * rand() - 1)  # Move particle to a new position
    pbc!(system, i)  # Apply periodic boundary conditions (this places the particle back in the box if it has moved outside of the valid range)
    Eafter = calc_energy(system, i)  # Calculate the potential energy of particle i after it has been moved
    if rand() > exp(-system.β * (Eafter - Ebefore))
        system.particles[i] = xbefore  # Trial move rejected. Reset particle to previous state
    end
end

function sweep(system::System; transitions=10, insert_delete_probability=0.2)
    for _ in 1:transitions
        if rand() < insert_delete_probability  # Randomly select trial transition: either move or insert/delete
            rand() < 0.5 ? trial_insert(system) : trial_delete(system)  # Do trial insertions and removals with equal probability
        else
            trial_move(system)
        end
    end
end

function simulate(L::Number, μ::Number, T::Number, Vext::Function, ϕ::Function; equilibration_time=Dates.Second(1), production_time=Dates.Second(2), sweep_transitions=10)
    system = System(L, μ, T, Vext, ϕ)  # The state of the system is encapsulated in this struct
    histograms = Histograms(system)
    ghist = gHistograms(system)
    equilibration_start = now()
    while now() - equilibration_start < equilibration_time  # Equilibration stage, no sampling
        sweep(system; transitions=sweep_transitions)
    end
    production_start = now()
    while now() - production_start < production_time
        sweep(system; transitions=sweep_transitions)
        sample(system, histograms)
        RDF(system, ghist)
    end
    get_results(system, histograms, ghist)  # Normalizes histograms and returns (xs, ρ)
end

function sample(system::System, histograms::Histograms)
    for x in system.particles
        bin = ceil(Int, x / L * histograms.bins)  # Calculate the bin index from the given particle position (x/L)*bins
        histograms.ρ[bin] += 1
    end
    push!(system.amount,length(system.particles))
    histograms.count += 1  # Needed later for normalization
end

function RDF(system::System, ghist::gHistograms)
    for i=1:length(system.particles)
        x = system.particles[i]
        for j = i+1:length(system.particles)
            dx = (x-system.particles[j])
            dx -= floor(dx/system.L)*system.L
            r = abs(dx)
            if r < L/2
                bin = ceil(Int, (2*r) / L * ghist.bins)
                ghist.g[bin] += 1
            end
        end           
    end
    ghist.count += 1
end
