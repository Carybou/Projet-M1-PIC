using Plots
using ProgressMeter


include("pic1d.jl")

kernel_d(pi, 1)

(x0, v0, wei) = samples(1000, 0.5, 0.5)
plot(v0)
plot(x0)

# +
function run(nstep)
    xmin =-pi
    xmax = pi
    vmin = -1
    vmax = 1
    k = 1
    nx = 64
    dt = 0.001
    mesh = Mesh1D(xmin, xmax, nx)
    (x0, y0, wei) = samples(nx, 0.5, 0.5)
    histogram(x0)
    p = Particles(x0, y0, wei)
    @gif for istep = 1:nstep # Loop over time
        update_velocities!(p, k, dt)
        update_positions!(p, mesh, dt)
        histogram(p.x)
    end every (nstep ÷ 100)
end

# -

run(2000)
