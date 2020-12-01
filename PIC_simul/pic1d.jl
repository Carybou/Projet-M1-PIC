using Sobol

export Mesh1D

"""
Describe the universe of the plasma : one dimension, between 2 bounds xmin and xmax
"""
struct Mesh1D

   xmin
   xmax
   nx #number of samples in the uniserve
   dx #distance between 2 samples

   function Mesh1D( xmin, xmax, nx)
       dx = (xmax - xmin) / nx
       new(xmin, xmax, nx, dx)
   end

end

export Particles

"""
Describe meta particules, represented by a Dirac disbtibution in (x, v), with a weight (~a high) of wei
"""
struct Particles

    x #list of the positions
    v #list of the velocities
    wei #list of the weights of the particules

    function Particles(x0, v0, wei)
        new(x0, v0, wei)
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

export update_positions!

"""
update particle position xp (phi_T)
"""
function update_positions!(p, mesh, dt)

    xmin = mesh.xmin
    xmax = mesh.xmax

    p.x .+= p.v .* dt
    #faire un modulo sans int

end

export update_velocities!

"""
update particle velocities vp (phi_v)
"""
function update_velocities!(p, k, dt)

    to_add=[]
    for position in p.x
        to_push = - dt * sum(p.wei .* kernel_d.(position .- p.x, k))
        push!(to_add, to_push)
    end
    p.v .+= to_add #reverse ?

end

"""
generated random samples
"""
function random_generation(nsamples, mesh, vmin, vmax, f0)
    x0 = []
    v0 = []
    wei = []
    s = SobolSeq(1)

    for i=1:nsamples
        x = Sobol.next!(s)
        x = pop!(x)
        x = (1-x)* mesh.xmin + x*mesh.xmax
        v = Sobol.next!(s)
        v = pop!(v)
        v = (1-v) * vmin + v*vmax
        push!(x0, f0(x))
        push!(v0, v)
        push!(wei, 1/nsamples)
    end
    return (x0, v0, wei)
end

"""
```math
f_0(x,v,t) = \\frac{n_0}{2π v_{th}^2} ( 1 + \\alpha cos(k_x x))
 exp( - \\frac{v^2}{2 v_{th}^2})
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
            p = x1 + alpha * sin( kx * x1) / kx
            f = 1 + alpha * cos( kx * x1)
            x1, x2 = x2, x1 - (p - r) / f
        end
        return x2
    end

    x0 = []
    v0 = []
    wei = []
    s = SobolSeq(1)

    for i=1:nsamples
        v = sqrt(-2 * log( (i-0.5)/nsamples))
        x = pop!(Sobol.next!(s))
        x = newton(x)
        push!(x0, x)
        push!(v0, v)
        push!(wei, 2*pi/kx/nsamples)
    end
    return x0, v0, wei
end
