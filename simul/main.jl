using Plots
using ProgressMeter
using SparseArrays
using LinearAlgebra

include("package.jl")

function run(nstep, dt, L, k, np, nx, alpha)

    kx = 2pi / L
    mesh = Mesh1D(L, nx)
    (x0, y0, wei) = samples(np, alpha, kx)
    p = Particules(x0, y0, wei, np)

    potential = []
    energy_elec_from_phi = []
    energy_pic = []
    energy_hamil_elec = []
    energy_hamil_cine = []
    energy_tot = []
    anime = @animate for istep = 1:nstep # Loop over time
        plot(size  = (500,500), ylim = (0,.5),
            widen = false, st=:surface, camera=(-30,30))
        histogram2d(real(p.x)[1:1:end], real(p.v)[1:1:end], bins = nx, normalize = true,
            xlabel = "position", ylabel = "velocity")
        energy_cine = 1/2 * sum(p.wei .* (p.v .* p.v))

        phi_t, der_phi_t = update_velocities!(p, k, dt, kx)
        push!(energy_pic, compute_int_E(p, mesh))
        update_positions!(p, mesh, dt)
        energy_elec = sum(p.wei .* phi_t)

        push!(potential, phi_t)
        push!(energy_elec_from_phi, sum(der_phi_t.^2)*L/np)
        push!(energy_hamil_elec, energy_elec)
        push!(energy_hamil_cine, energy_cine)
        push!(energy_tot, 1/2 * energy_elec + energy_cine) #pb avec le 1/2
    end
    return energy_hamil_elec, energy_pic, energy_elec_from_phi, energy_tot, anime
end

nstep = 100
dt = 0.2
L = 2pi/0.5
k = 1 #paramètre de troncature du noyau
np = 100000 #nb de particules
nx = 100   #nb de mailles pour le calcul de l energie
alpha = 0.1 #perturbation initiale

energy_hamil_elec, energy_pic, energy_elec_from_phi, energy_tot, anime = run(nstep, dt, L, k, np, nx, alpha)

gif(anime)

time = [i*dt for i =1:nstep]
plot(time, energy_tot)
plot(time, log.(sqrt.(energy_pic)), label="PIC", xlabel = "time")
plot!(time, log.(sqrt.(max.(energy_hamil_elec, 0.0000000001))), label = "Hamiltonian")
plot!(time, log.(sqrt.(energy_elec_from_phi)), label = "Potential")
