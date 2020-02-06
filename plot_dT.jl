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

# This is either the raw HTC or a spline fit.
htc_file = "Experiment4/checkpoint_htc_360.csv"
htc_recovered = read_htc(htc_file)

# Experimental LN2 HTC
exp_file = "ln2.csv"
exp_recovered = read_htc(exp_file)

const testfile = "Fig_67.csv_Result.csv"
const test_data = read_data(testfile)

const N = 10
const t = 0.00635
const dx = t / N

const sinkT = test_data[end, 2]
const ambT = test_data[1, 2]
const u20 = [ambT for i in 1:N]

tspan = (Float64(test_data[1, 1]), Float64(test_data[end, 1]))

const pr = hprob(450.0, 7850.0, 42.0, sinkT, ambT, 4.0)

# Solve ODE.
const p = Vector(htc_recovered[2])

prob = ODEProblem(heat_transfer, u20, tspan, p, saveat = 1.0)
sol = solve(prob, Tsit5(), saveat = 1.0)
#plot(sol.t[:], sol[1, :], label = "Exposed")
#plot!(sol.t[:], sol[end, :], label = "Unexposed")

exposedT = sol[1, :]
vapourT = (exposedT .+ sinkT) / 2
ΔT = abs.(sinkT .- vapourT)
plot(xlims = (1, 120))
#plot!(xscale = :log, yscale = :log)
plot!(xlabel = "\\Delta T", ylabel = "HTC W/mK")
plot!(ΔT, p, label = "Recovered")
plot!(exp_recovered[1], exp_recovered[2], label = "Small scale experiment")
