using Plots
using ProgressMeter
using SparseArrays
using LinearAlgebra


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

histogram!(x0, bins = 50, normalize = true)

nx = 6
matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                            0 => -2*ones(Float64,nx+1),
                            1 => ones(Float64,nx-1))
det(matrix_poisson)"""

# +
function run(nstep)
    xmin = 0
    xmax = 4pi
    k = 5 #paramètre de troncature du noyau 
    np = 10000 #nb de particules
    nx = 300   #nb de mailles pour le calcul de l energie

    dt = 0.1
    mesh = Mesh1D(xmin, xmax, nx)
    (x0, y0, wei) = samples(np, 0.5, 0.5)
    p = Particules(x0, y0, wei, np)


    matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                                0 => -2*ones(Float64,nx+1),
                                1 => ones(Float64,nx-1))

    E = []
    rho_i, rho_ti = compute_rho(p, mesh)
    phi_i = compute_phi(rho_i, rho_ti, mesh, matrix_poisson)
    e = compute_E(phi_i, mesh)
    potential = []
    anime = @animate for istep = 1:nstep # Loop over time
        push!(potential, update_velocities_S_C!(p, k, dt))
        update_positions!(p, mesh, dt)
        plot(size  = (500,500), xlim = (xmin, xmax), ylim = (0,.5),
                legend=false, widen = false)
        histogram!(real(p.x), bins = nx, normalize = true)
        push!(E, compute_int_E(p, mesh, matrix_poisson))
    end
    #plot(potential)
    """pot = pop!(potential)
    plot(size  = (500,500),
            legend=false, widen = false)
    plot!(real(p.x), real(pot), seriestype = :scatter) #sin x /0.4 pi + dephasage, phi
    rho, rho_t = compute_rho(p, mesh)
    phi = compute_phi(rho, mesh, matrix_poisson)
    plot(phi_i)
    plot!(phi)
    plot(E, yaxis=:log)"""
    plot(p.x, pop!(potential), seriestype = :scatter)
    #plot!(phi_i)
    #gif(anime)
end

# -

run(100)
