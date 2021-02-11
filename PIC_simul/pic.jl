using Plots
using ProgressMeter
using SparseArrays
using LinearAlgebra

include("pic1d.jl")
"""
xmin =0
xmax = 4pi
k = 1
nx = 30

dt = 0.001
m = Mesh1D(xmin, xmax, nx)
(x0, y0, wei) = samples_2(1000, xmin, xmax)
p = Particules(x0, y0, wei, 1000)
r, rt = compute_rho(p, m)
plot(r)
plot(y0)
plot(size  = (500,500), widen  = false)

#fit(Histogram, (x0, y0), bins = 50)
histogram2d (x0[1:1:end], y0[1:1:end], bins = 100)
plot(x0, y0)

nx = 6
matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                            0 => -2*ones(Float64,nx+1),
                            1 => ones(Float64,nx-1))
det(matrix_poisson)"""

# +
function run(nstep)
    xmin = 0
    xmax = 4pi
    k = 1 #paramètre de troncature du noyau
    np = 10000 #nb de particules
    nx = 100   #nb de mailles pour le calcul de l energie

    dt = 0.1
    alpha = 0.5
    kx = 2pi /(xmax -xmin)
    mesh = Mesh1D(xmin, xmax, nx)
    (x0, y0, wei) = samples_2(np, alpha, kx)
    p = Particules(x0, y0, wei, np)


    matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                                0 => -2*ones(Float64,nx+1),
                                1 => ones(Float64,nx-1))

    rho_i, rho_ti = compute_rho(p, mesh)
    phi_i = compute_phi(rho_i, rho_ti, mesh, matrix_poisson)
    e = compute_E(phi_i, mesh)
    potential = []
    energy_pic = []
    energy_hamil_elec = []
    energy_hamil_cine = []
    energy_tot = []
    anime = @animate for istep = 1:nstep # Loop over time
        plot(size  = (500,500), ylim = (0,.5),
            widen = false, st=:surface, camera=(-30,30))
        histogram2d(real(p.x)[1:1:end], real(p.v)[1:1:end], bins = nx, normalize = true)

        energy_cine = 1/2 * sum(p.wei .* (p.v .* p.v))

        phi_t = update_velocities!(p, k, dt, kx)
        push!(energy_pic, compute_int_E(p, mesh, matrix_poisson))
        update_positions!(p, mesh, dt)
        energy_elec = sum(p.wei .* phi_t)

        push!(potential, phi_t)
        push!(energy_hamil_elec, energy_elec)
        push!(energy_hamil_cine, energy_cine)
        push!(energy_tot, energy_elec + energy_cine)
    end
    #plot(potential)
    """rho, rho_t = compute_rho(p, mesh)
    plot(real(p.x), real(pot), seriestype = :scatter)
    phi = compute_phi(rho, rho_t, mesh, matrix_poisson)
    plot(phi)
    plot!(energy_hamil)
    """
    plot(real(rho_i))
    plot!(phi_i)
    #plot!(phi_i)
    #plot(p.x, pop!(potential), seriestypes = :scatter)
    plot(energy_pic)
    plot!(energy_hamil_elec)
    plot!(energy_tot)
    gif(anime)
end

# -

run(300)
# faire distribution / check 3 méthodes de calcul de phi (PIC (ok), hamiltonien (ok), exact et th)
# coeff dans calcul hamiltonien avec TFD N(x) = \sum 1/(2pi*(k*kx)^2) exp(i*kx*k*x) x \in [0, 2pi/kx],
# qui est C ? (1/2pi ?)
