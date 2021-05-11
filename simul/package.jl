using Sobol
using Random
using Distributions

"""
Describe the universe of the plasma : one dimension, between 2 bounds 0 and L
"""
struct Mesh1D

   L
   nx #number of samples in the uniserve
   dx #distance between 2 samples

   function Mesh1D( L, nx)
       dx = L / nx
       new(L, nx, dx)
   end

end

"""
Describe meta particules, represented by a Dirac disbtibution in (x, v), with a weight (~a high) of wei
"""
struct Particules

    x #list of the positions
    v #list of the velocities
    wei #list of the weights of the particules
    nbpart #nulber of particules

    function Particules(x0, v0, wei, nbpart)
        new(x0, v0, wei, nbpart)
    end
end

"""
Define the kernel of the model, is n = 1 kernel = cos (HMF model)
"""
function kernel(x, n)
    return(sum(1/k^2 * exp(im * k * x) for k in -1*n:n if k != 0))
end

"""
Define the derivate kernel of the model, is n = 1 kernel = cos (HMF model)
"""
function kernel_d(x, n)
    return(sum(im/k * exp(im * k * x) for k in -1*n:n if k != 0))
end
""" attention ffttd donc facteur dans la somme  a voir """

""" S[k] = \\sum_l=1^n {\\beta_l * sin(k kx x_l)} et C[k] = \\sum_l=1^n {\\beta_l * cos(k kx x_l)}
utile pour le calcul de la vitesse et du potentiel electrique phi"""
function compute_S_C(p, n, kx)
    S = []
    C = []
    for k in 1:n
        push!(S, sum(p.wei .* sin.((k * kx) .* p.x)))
        push!(C, sum(p.wei .* cos.((k * kx) .* p.x)))
    end
    return S, C
end

"""
update particle position xp (phi_T)
"""
function mod_float(p, L)
    new_p = []
    for x in p
        if Real(x) > L
            push!(new_p,x - L)
        elseif real(x) < 0
            push!(new_p,L + x)
        else push!(new_p,x)
        end
    end
    return new_p
end

function update_positions!(p, mesh, dt)
    L = mesh.L

    p.x .+= p.v .* dt
    p.x .= mod_float(p.x, L)
end

"""
update particle velocities vp (phi_v)
"""
function update_velocities!(p, n, dt, kx)
    S, C = compute_S_C(p, n, kx)
    der_phi = []
    phi = []
    for position in p.x
        der_phi_x = - sum([1/(pi*k) * (sin(k*kx*position) * C[k] -
                                        cos(k*kx*position) * S[k]) for k in 1:n])
        phi_x = sum([1/(pi*k^2*kx) * (cos(k*kx*position) * C[k] -
                                    sin(k*kx*position) * S[k]) for k in 1:n])
        push!(der_phi, der_phi_x)
        push!(phi, phi_x)
    end
    p.v .+=  - dt .* der_phi
    return phi, der_phi
end

"""
```math
f_0(x,v,t) = \\frac{1}{2π} ( 1 + \\alpha cos(k_x x))exp( - \\frac{v^2}{2})
```
The newton function solves the equation ``P(x)-r=0`` with Newton’s method
```math
x^{n+1} = x^n – (P(x)-(2\\pi r / k)/f(x)
```
with
```math
P(x) = \\int_0^x (1 + \\alpha cos(k_x y)) dy = x + \\frac{\\alpha}{k_x} sin(k_x x)
```
"""

function samples(nsamples, alpha, kx)
    function newton(r)
        x1, x2 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x2-x1) > 1e-12)
            p = x1 + alpha * sin( kx * x1) / kx #primitive 1 + alpha cos(kx x)
            f = 1 + alpha * cos( kx * x1)
            x1, x2 = x2, x1 - (p - r) / f
        end
        return x2
    end
    x0 = []
    s = SobolSeq(1)
    for i=1:nsamples
        x = pop!(Sobol.next!(s))
        x = newton(x)
        push!(x0, x)
    end
    distr_v = Normal()
    v0 = rand(distr_v, nsamples)
    wei = [2*pi/kx/nsamples for i=1:nsamples]
    return x0, v0, wei
end

"""
compute rho, dcharge density (ie int f dv)
"""
function compute_rho( p, m )

   L = m.L
   nx = m.nx
   dx = m.dx
   rho = zeros(nx + 1)

   for ipart=1:p.nbpart
      xp = p.x[ipart]
      i = Int(floor(Real(xp) / dx) + 1) #number of the cell with the particule
      dum = p.wei[ipart] / dx
      mx_i = (i-1) * dx
      mx_j = i * dx
      a1 = (mx_j-xp) * dum
      a2 = (xp-mx_i) * dum
      rho[i]   += a1
      rho[i+1] += a2
   end

   rho[nx+1] += rho[1]
   rho[1] = rho[nx + 1]

   rho ./= dx

   rho_total  = sum(rho[1:nx]) * dx / L
   return rho, rho_total
end

"""
int E^2(t,x) dx entre 0 et L (energie elec, dont on connait l evolution temp)
    avec E = - der_x phi, der_x E = rho_tot, calculer avec poisson
    - laplacien phi = rho -> difference finie
"""

function compute_phi(rho, rho_total, mesh)
    dx = mesh.dx
    L = mesh.L
    nx = mesh.nx
    matrix_poisson = spdiagm(-1 => ones(Float64,nx-1),
                                0 => -2*ones(Float64,nx+1),
                                1 => ones(Float64,nx-1))
    phi = matrix_poisson \ ((rho .- rho_total) .* -dx^2 )
        # - rho_total car le monde est circulaire, sol périodique
    phi_total = sum(phi) * dx / L
    return phi .- phi_total # - phi_total densité car centré en 0
end

function compute_E(phi, mesh)
    dx = mesh.dx
    e = 0 .* phi
    e[1:end] .= phi[1:end]
    e .= (circshift(e, 1) .- circshift(e, -1)) ./ (2dx)
    return e
end

function compute_int_E(p, mesh)
    dx = mesh.dx
    rho, rho_t = compute_rho(p, mesh)
    phi = compute_phi(rho, rho_t, mesh)
    E = compute_E(phi, mesh)
    return sum(E .^2) * dx
end

"""
calcul energie elec potentielle et cinétique avec le hamiltonien

debug pas de temps de 10 ité pour voir l evolution de E et H
"""
