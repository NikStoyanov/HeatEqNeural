using CSV
using Flux
using Plots
using DiffEqFlux
using DifferentialEquations

N_tot = 20
tot_t = 1.0 # TODO: get plate thickness.
dx = tot_t / N_tot

# TODO: confirm sink and ambient.
sinkT = -190.0
ambT = 20.0
u20 = ambT * ones(N_tot)
h = 40 * ones(180) # TODO: check if closer initial guess is possible.
tspan = (0.0, 180.0)

mutable struct hprob
    h # Vector with size N = number of test data points.
    c # Specific heat.
    ρ # Density.
    k # Thermal conductivity.
    sinkT # Cryogenic liquid temperature.
    ambT # Initial condition.
    htcAmb # Atmospheric HTC. TODO: Refer to study from PPG document to justify value.
end

pr = hprob(h, 450.0, 7850.0, 44.0, sinkT, ambT, 4.0)

# Calculate h as the interpolation for the current time point.
function htc(ti)
    y1 = y2 = 0.0
    x1 = x2 = 0

    for (t, h) in enumerate(pr.h)
        if ti >= t
            x1 = t
            y1 = h
        end
    end

    if x1 == length(pr.h)
        return y1
    else
        x2 = x1 + 1
        y2 = pr.h[x2]
        return (y2 - y1) / (x2 - x1) * (ti - x2) + y1
    end
end

# TODO: Parametrize solve to use p.
function heat_transfer(du, u, p, t)
    ΔT = (pr.sinkT + u[1]) / 2 - u[1]
    α = pr.k / (pr.ρ * pr.c)

    # Node exposed to cryogenic liquid.
    du[1] = 2 * (ΔT * htc(t) /
                 (pr.c * pr.ρ) - α * (u[1] - u[2]) / dx) / dx

    # Inner nodes.
    for i in 2:N_tot - 1
        du[i] = α * (u[i - 1] - 2u[i] + u[i + 1]) / dx^2
    end

    # Node exposed to ambient atmosphere. TC is located here.
    du[end] = 2 * (α * (u[end-1] - u[end]) /
                   dx - (u[end] - pr.ambT) * pr.htcAmb /
                   (pr.ρ * pr.c)) / dx
end

# TODO: Setup https://github.com/FluxML/model-zoo/blob/da4156b4a9fb0d5907dcb6e21d0e78c72b6122e0/other/diffeq/ode.jl
# Solve ODE.
function direct_problem()
    val_prob = ODEProblem(heat_transfer, u20, tspan)
    val_sol = solve(val_prob, Tsit5())
    plot(val_sol)
end

# Read test data.
function read_data()
    filename = "Fig_67.csv_Result.csv"
    CSV.read(filename; skipto = 1)
end

direct_problem()
