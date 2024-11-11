using Flux.Zygote

function read_sim_data(dir)
    ρ_profiles = Vector{Vector{Float64}}()
    c1_profiles = Vector{Vector{Float64}}()
    ϕ_profiles = Vector{Vector{Float64}}()
    for sim in readdir(dir, join=true)
        if !endswith(sim, ".dat") #check if data is of .dat type
            continue        # if file not of .dat type, do not execute loop for this file
        end
        xs, μloc, ρ, phi = eachcol(readdlm(sim))
        ϕ = phi[1:150]
        c1 = log.(ρ) .- μloc
        push!(ρ_profiles, ρ)
        push!(c1_profiles, c1)
        push!(ϕ_profiles, ϕ)
    end
    ρ_profiles, c1_profiles, ϕ_profiles
end

function generate_windows(ρ; window_bins)
    ρ_windows = Zygote.Buffer(zeros(Float32, window_bins, length(ρ))) # We use a Zygote Buffer here to keep autodifferentiability
    pad = window_bins ÷ 2 - 1 # a number
    ρpad = vcat(ρ[end-pad:end], ρ, ρ[1:1+pad]) 
    for i in 1:length(ρ)
        ρ_windows[:,i] = ρpad[i:i+window_bins-1]
    end
    copy(ρ_windows)  # copy needed due to Zygote.Buffer
end

function generate_phi(ϕ, ρ)
    ϕ_func = Zygote.Buffer(zeros(Float32, length(ϕ), length(ρ)))
    for i in 1:length(ρ)
        ϕ_func[:,i] = ϕ
    end
    copy(ϕ_func)
end

function generate_inout(ρ_profiles, c1_profiles, ϕ_profiles; window_width, dx)
    window_bins = 2 * round(Int, window_width / dx) + 1
    ρ_windows_all = Vector{Vector{Float32}}()
    c1_values_all = Vector{Float32}()
    ϕ_functions_all = Vector{Vector{Float32}}()
    for (ρ, c1, ϕ) in zip(ρ_profiles, c1_profiles, ϕ_profiles)
        ρ_windows = generate_windows(ρ; window_bins)
        ϕ_func = generate_phi(ϕ,ρ)
        for i in collect(1:5:length(c1)) 
            if !isfinite(c1[i])
                continue
            end
            push!(ρ_windows_all, ρ_windows[:,i])
            push!(c1_values_all, c1[i])
            push!(ϕ_functions_all, ϕ_func[:,i])
        end
    end
    reduce(hcat, ρ_windows_all), c1_values_all', reduce(hcat, ϕ_functions_all)
end

function ϕ(x, params) 
    l = length(params)
    Δ = 1.5/(l-1)
    if 0 <= x < 1.5
        bin = x ÷ Δ
        bin = Integer(bin)
        ϵ1 = params[bin+1]
        ϵ2 = params[bin+2]
        return (ϵ2-ϵ1)*(x-bin*Δ)/(Δ) + ϵ1
    else
        return 0
    end
end

function get_params(f::Function)
    grid = collect(0:0.01:1.5-0.01)
    params = f.(grid)
    for i=1:length(params)
        if params[i] > 9 || isinf(params[i]) == true 
            params[i] = 9
        end
    end 
    for i=1:length(params)
        if isnan(params[i]) == true 
            params[i] = 9
        end
    end 
    return params
end



###############################################################################################################
# DFT

function minimize(L::Number, μ::Number, T::Number, ϕ::Vector{Float32}, Vext::Function, get_c1::Function; α::Number=0.03, maxiter::Int=10000, dx::Number=0.01, floattype::Type=Float32, tol::Number=max(eps(floattype(1e3)), 1e-8))
    L, μ, T = floattype.((L, μ, T))
    ϕ = floattype.(ϕ)
    xs = collect(floattype, dx/2:dx:L)  # Construct the numerical grid
    Vext = Vext.(xs)  # Evaluate the external potential on the grid
    infiniteVext = isinf.(Vext)  # Check where Vext is infinite to set ρ = 0 there
    ρ, ρEL = zero(xs), zero(xs)  # Preallocate the density profile and an intermediate buffer for iteration
    fill!(ρ, 0.5)  # Start with a bulk density of 0.5
    c1 = get_c1(xs)  # Obtain the c1 functional for the given numerical grid
    i = 0
    while true
        ρEL .= exp.((μ .- Vext) ./ T .+ c1(ρ,ϕ))  # Evaluate the RHS of the Euler-Lagrange equation
        ρ .= (1 - α) .* ρ .+ α .* ρEL  # Do a Picard iteration step to update ρ
        ρ[infiniteVext] .= 0  # Set ρ to 0 where Vext = ∞
        clamp!(ρ, 0, Inf)  # Make sure that ρ does not become negative
        Δρmax = maximum(abs.(ρ - ρEL)[.!infiniteVext])  # Calculate the remaining discrepancy to check convergence
        i += 1
        if Δρmax < tol
            println("Converged (step: $(i), ‖Δρ‖ = $(Δρmax) < $(tol) = tolerance)")
            break  # The remaining discrepancy is below the tolerance: break out of the loop and return the result
        end
        if !isfinite(Δρmax) || i >= maxiter
            println("Did not converge (step: $(i) of $(maxiter), ‖Δρ‖: $(Δρmax), tolerance: $(tol))")
            return nothing  # The iteration did not converge, there is no valid result
        end
    end
    xs, ρ
end

function MINIMIZE(L::Number, μ::Number, T::Number, ϕ::Function, Vext::Function, model)
    s = 1/T
    G(x) = s*Vext(x)
    Φ(x) = s*ϕ(x)
    params = get_params(Φ)
    params = Float32.(params)
    μ = s*μ
    T = 1
    minimize(L, μ, T, params, G, xs -> get_c1_neural(model, params))
end








