Vext_sin(x; n::Int, A::Number, φ::Number, L::Number) = A * sin(2π * x * n / L + φ)

Vext_lin(x; x1::Number, x2::Number, E1::Number, E2::Number) = x > x1 && x < x2 ? E1 + (x - x1) * (E2 - E1) / (x2 - x1) : 0

Vext_wall(x; xw::Number, L::Number) = x < xw || x > L - xw ? Inf : 0

function generate_Vext(L::Number; num_sin=4, num_lin=rand(1:5), wall=true)
    Avar = 1.0
    sin_parameters = []
    for n in 1:num_sin  # Generate random parameters for periodic sine functions with increasing frequency
        push!(sin_parameters, (n = n, A = randn() * Avar, φ = rand() * 2π, L = L))
    end
    Evar = 1.0
    lin_parameters = []
    for _ in 1:num_lin  # Generate random parameters for discontinuous linear segments
        push!(lin_parameters, (x1 = round(rand() * L, digits=2), x2 = round(rand() * L, digits=2), E1 = randn() * Evar, E2 = randn() * Evar))
    end
    xwmax = 1.0
    wall_params = (xw = round(rand() * xwmax, digits=2), L = L)  # Set a random wall width
    function (x)  # Return a method which evaluates a combination of all functions with the chosen parameters above
        result = 0.0
        for sin_params in sin_parameters
            result += Vext_sin(x; sin_params...)
        end
        for lin_params in lin_parameters
            result += Vext_lin(x; lin_params...)
        end
        if wall
            result += Vext_wall(x; wall_params...)
        end
        result
    end
end

function Generateparams(L::Number; num_sin=4, num_lin=rand(1:5), wall=true)
    Avar = 1.0
    sin_parameters = []
    for n in 1:num_sin  # Generate random parameters for periodic sine functions with increasing frequency
        push!(sin_parameters, (n = n, A = randn() * Avar, φ = rand() * 2π, L = L))
    end
    Evar = 1.0
    lin_parameters = []
    for _ in 1:num_lin  # Generate random parameters for discontinuous linear segments
        push!(lin_parameters, (x1 = round(rand() * L, digits=2), x2 = round(rand() * L, digits=2), E1 = randn() * Evar, E2 = randn() * Evar))
    end
    xwmax = 1.0
    wall_params = (xw = round(rand() * xwmax, digits=2), L = L)  # Set a random wall width
    return sin_parameters, lin_parameters, wall_params
end

function GenerateV(x, sin_parameters, lin_parameters, wall_params, wall=true)
    result = 0.0
    for sin_params in sin_parameters
        result += Vext_sin(x; sin_params...)
    end
    for lin_params in lin_parameters
        result += Vext_lin(x; lin_params...)
    end
    if wall
        result += Vext_wall(x; wall_params...)
    end
    result
end
