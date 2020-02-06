# Determine the mesh sensitivity for the default unexposed HTC and the recovered value of
# the exposed HTC. This is ran post-factum of the recovery calculation to justify it.
# In case the correct mesh was not picked the adjoint calculation muts be repeated with
# the converged value of the mesh size.

using CSV
using Plots
using DataFrames
using DifferentialEquations

struct hprob
    c # Specific heat.
    ρ # Density.
    k # Thermal conductivity.
    sinkT # Cryogenic liquid temperature.
    ambT # Initial condition.
    htcAmb # Atmospheric HTC.
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
    for i in 2:size(u)[1] - 1
        du[i] = α * (u[i - 1] - 2u[i] + u[i + 1]) / dx^2
    end

    # Node exposed to ambient atmosphere. TC is located here.
    du[end] = 2 * (α * (u[end-1] - u[end]) /
                   dx - (u[end] - pr.ambT) * pr.htcAmb /
                   (pr.ρ * pr.c)) / dx
end

# Read test data.
function read_data(filename)
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

# Read recovered htc.
function read_htc(filename)
    data = CSV.read(filename)
    htc_recovered = DataFrame(data)

    # Overwrite non-physical values i.e. < 0
    nrows, ncols = size(htc_recovered)
    for row in 1:nrows
        if htc_recovered[row, 2] <= 0.0
            println("Row ", row, " is < 0. Overwrite.")
            htc_recovered[row, 2] = 1.0
        end
    end

    return htc_recovered
end

const testfile = "Fig_67.csv_Result.csv"
const test_data = read_data(testfile)

# This is either the raw HTC or a spline fit.
const htc_file = "Experiment4/checkpoint_htc_360.csv"
const htc_recovered = read_htc(htc_file)

const sinkT = test_data[end, 2]
const ambT = test_data[1, 2]

tspan = (Float64(test_data[1, 1]), Float64(test_data[end, 1]))

const pr = hprob(450.0, 7850.0, 42.0, sinkT, ambT, 4.0)

const p = Vector(htc_recovered[2])

# Pick correct value.
const t = 0.00635

# Study the mesh sensitivity.
const N1 = 10
const N2 = 20
const N3 = 40
const N4 = 80

dx = t / N1
u = [ambT for i in 1:N1]
prob = ODEProblem(heat_transfer, u, tspan, p, saveat = 1.0)
sol1 = solve(prob, Tsit5(), saveat = 1.0)

dx = t / N2
u = [ambT for i in 1:N2]
prob = ODEProblem(heat_transfer, u, tspan, p, saveat = 1.0)
sol2 = solve(prob, Tsit5(), saveat = 1.0)

dx = t / N3
u = [ambT for i in 1:N3]
prob = ODEProblem(heat_transfer, u, tspan, p, saveat = 1.0)
sol3 = solve(prob, Tsit5(), saveat = 1.0)

dx = t / N4
u = [ambT for i in 1:N4]
prob = ODEProblem(heat_transfer, u, tspan, p, saveat = 1.0)
sol4 = solve(prob, Tsit5(), saveat = 1.0)

plot(sol1.t[:], sol1[end, :], label = "N = 10", color = :black, linestyle = :dash)
plot!(sol2.t[:], sol2[end, :], label = "N = 20", color = :black, linestyle = :dot)
plot!(sol3.t[:], sol3[end, :], label = "N = 40", color = :black, linestyle = :dashdot)
plot!(sol4.t[:], sol4[end, :], label = "N = 80", color = :black)
plot!(fmt = :svg,
      grid = false,
      xlabel = "Time (min)",
      ylabel = "Temperature (°C)",
      xlims = (0, 200),
      ylims = (-200, 50),
      xticks = 0:25:200,
      yticks = -200:25:50,
      minorticks = true)

savefig("mesh_sensitivity.svg")
