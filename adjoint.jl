using CSV
using Flux
using Plots
using DataFrames
using DiffEqFlux
using DifferentialEquations

mutable struct hprob
    c # Specific heat.
    ρ # Density.
    k # Thermal conductivity.
    sinkT # Cryogenic liquid temperature.
    ambT # Initial condition.
    htcAmb # Atmospheric HTC. TODO: Refer to study from PPG document to justify value.
end

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
    for i in 2:N - 1
        du[i] = α * (u[i - 1] - 2u[i] + u[i + 1]) / dx^2
    end

    # Node exposed to ambient atmosphere. TC is located here.
    du[end] = 2 * (α * (u[end-1] - u[end]) /
                   dx - (u[end] - pr.ambT) * pr.htcAmb /
                   (pr.ρ * pr.c)) / dx
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

# Flux reverse-mode AD through the differential equation solver.
function predict_rd()
    diffeq_adjoint(p, prob, Tsit5(), saveat = 1.0)[end, :]
end

# Least squares error.
function loss_rd()
    sum(abs2, test_data[:, 2] .- predict_rd())
end

test_data = read_data()

N = 20
t = 0.00635
dx = t / N

sinkT = test_data[end, 2]
ambT = test_data[1, 2]
u20 = ambT * ones(N)
tspan = (Float64(test_data[1, 1]), Float64(test_data[end, 1]))

pr = hprob(450.0, 7850.0, 44.0, sinkT, ambT, 4.0)

# Solve ODE.
p = param(550 * ones(test_data[end, 1])) # TODO: check if closer initial guess is possible.
#p = 550 * ones(test_data[end, 1])
prob = ODEProblem(heat_transfer, u20, tspan, p, saveat = 1.0)

#sol = solve(prob)
#plot(sol.t[:], sol[end, :])
#plot!(test_data[:, 1], test_data[:, 2])

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cb = function ()
    display(show(loss_rd()))
    sol = solve(remake(prob, p = Flux.data(p)), Tsit5(), saveat = 1.0)
    display(plot(sol.t[:], sol[end, :]))
    display(plot!(test_data[:, 1], test_data[:, 2]))
end

cb()
Flux.train!(loss_rd, [p], data, opt, cb = cb)
