using CSV
using Flux
using Plots
using DataFrames
using DiffEqFlux
using DifferentialEquations

N_tot = 20
tot_t = 1.0 # TODO: get plate thickness.
dx = tot_t / N_tot

# TODO: confirm sink and ambient.
sinkT = -190.0
ambT = 20.0
u20 = ambT * ones(N_tot)
tspan = (0.0, 180.0)

mutable struct hprob
    c # Specific heat.
    ρ # Density.
    k # Thermal conductivity.
    sinkT # Cryogenic liquid temperature.
    ambT # Initial condition.
    htcAmb # Atmospheric HTC. TODO: Refer to study from PPG document to justify value.
end

pr = hprob(450.0, 7850.0, 44.0, sinkT, ambT, 4.0)

# Calculate h as the interpolation for the current time point.
function htc(h, ti)
    y1 = y2 = 0.0
    x1 = x2 = 0

    for (t, h) in enumerate(h)
        if ti >= t
            x1 = t
            y1 = h
        end
    end

    if x1 == length(h)
        return y1
    else
        x2 = x1 + 1
        y2 = h[x2]
        return (y2 - y1) / (x2 - x1) * (ti - x2) + y1
    end
end

# Discretise ODE.
function heat_transfer(du, u, h, t)
    ΔT = (pr.sinkT + u[1]) / 2 - u[1]
    α = pr.k / (pr.ρ * pr.c)

    # Node exposed to cryogenic liquid.
    du[1] = 2 * (ΔT * htc(h, t) /
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
    h = 40 * ones(180) # TODO: check if closer initial guess is possible.
    val_prob = ODEProblem(heat_transfer, u20, tspan, h)
    val_sol = solve(val_prob, Tsit5())
    plot(val_sol)
end

# Read test data.
function read_data()
    filename = "Fig_67.csv_Result.csv"
    data = CSV.read(filename)
    table = DataFrame(data)

    # Drop duplicates.
    nrows, ncols = size(table)
    for row in 2:nrows
        if table[row, 1] == table[row - 1, 1]
            println("Row ", row, " is a duplicate. Deleting...")
            table[row, 1] = 0
        end
    end

    table = table[table[1].!=0,:]

    return table
end

data = read_data()
