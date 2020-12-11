using Plots
using ProgressMeter
using SparseArrays


include("pic1d.jl")

"""xmin =0
xmax = 4pi
k = 1
nx = 30

dt = 0.001
m = Mesh1D(xmin, xmax, nx)
(x0, y0, wei) = samples(1000, 0.5, 0.5)
p = Particules(x0, y0, wei, 1000)
r, rt = compute_rho(p, m)
plot(r)
plot(y0)
plot(size  = (500,500),
        legend=false, widen  = false)

histogram!(x0, bins = 50, normalize = true)"""

# +
function run(nstep)
    xmin = 0
    xmax = 4pi
    k = 1
    nx = 1000

    dt = 0.001
    mesh = Mesh1D(xmin, xmax, nx)
    (x0, y0, wei) = samples(nx, 0.5, 0.5)
    p = Particules(x0, y0, wei, nx)

    matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                                0 => -2*ones(Float64,nx+1),
                                1 => ones(Float64,nx-1))

    E = []

    anime = @animate for istep = 1:nstep # Loop over time
        update_velocities!(p, k, dt)
        update_positions!(p, mesh, dt)
        plot(size  = (500,500), xlim = (xmin, xmax), ylim = (0,.5),
                legend=false, widen = false)

        histogram!(real(p.x), bins = Int(nx/100), normalize = true)
        push!(E, compute_int_E(p, mesh, matrix_poisson))
    end
    #gif(anime)
    plot(E)
end

# -

run(2000)
