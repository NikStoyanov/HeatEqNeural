using CSV
using Flux
using Plots
using DiffEqFlux
using DifferentialEquations

N_tot = 20
tot_t = 1.0 # TODO: get plate thickness.
dx = tot_t / N_tot
val_thick_20 = range(0, stop=tot_t, length=N_tot)

# TODO: confirm sink and ambient.
sintT = -190.0
ambT = 20.0
u20 = ambT * ones(N_tot)
h = zeros(N_tot) # TODO: check if closer initial guess is possible.

mutable struct hprob
    h # Vector with size N = number of test data points.
    c # Specific heat.
    ρ # Density.
    k # Thermal conductivity.
    sinkT # Cryogenic liquid temperature.
    ambT # Initial condition.
    htcAmb # Atmospheric HTC. TODO: Refer to study from PPG document to justify value.
end

# Calculate h as the interpolation for the current time point.
function htc(ti)
    y1, y2, x1, x2 = 0.0
    for (t, h) in enumerate(hprob.h)
        if ti >= t
            x1 = ti
            y1 = h
        end
    end

    # TODO: handle exception.
    if target != 0.0
        x2 = x1 + 1.0
        y2 = hprob.h[Int(x2)]
    end

    return (y2 - y1) / (x2 - x1) * (ti - x2) + y1
end

# TODO: Solve ODE.
function heat_transfer(du, u, p, t)
    ΔT = (hprob.sinkT + u[1]) / 2 - u[1]

    # Node exposed to cryogenic liquid.
    du[1] = 2 * (ΔT * htc(u.t) /
                 (prot_c * prot_ρ) - α1 * (u[1] - u[2]) / dx) / dx

    # Inner nodes.
    for i in 2:N_prot - 1
        du[i] = α1 * (u[i - 1] - 2u[i] + u[i + 1]) / dx^2
    end

    # Node exposed to ambient atmosphere. TC is located here.
    du[end] = 2 * (α2 * (u[end-1] - u[end]) /
                   dx - (u[end] - ambT) * htcAmbient /
                   (steel_c * steel_ρ)) / dx
end

# TODO: Setup https://github.com/FluxML/model-zoo/blob/da4156b4a9fb0d5907dcb6e21d0e78c72b6122e0/other/diffeq/ode.jl
# Solve ODE.
function direct_problem()
    val_prob = ODEProblem(heat_transfer, u20, tspan, p)
    val_sol = solve(val_prob, Tsit5())
end

# Read test data.
function read_data()
    filename = "Fig_67.csv_Result.csv"
    CSV.read(filename; skipto = 1)
end
