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

export Particles1D

"""
Describe meta particules, represented by a Dirac disbtibution in (x, v), with a weight (~a high) of wei
"""
struct Particles1D

    wei #list of the weights of the particules
    x #list of the positions
    v #list of the velocities

end

"""
Define the kernel of the model, is n = 1 kernel = cos (HMF model)
"""
function kernel(x, n)
    sum(1/k^2 * exp(im * k * x) for k in -1*n:n if k != 0)
end

"""
Define the derivate kernel of the model, is n = 1 kernel = cos (HMF model)
"""
function kernel_d(x, n)
    sum(im/k * exp(im * k * x) for k in -1*n:n if k != 0)
end

export update_positions!

"""
update particle position xp
"""
function update_positions!(p, mesh, dt)

    xmin = mesh.xmin
    xmax = mesh.xmax

    p.xp .+= p.vp .* dt
    p.xp .= xmin .+ mod.( p.xp .- xmin , xmax - xmin)

end

export update_velocities!

"""
update particle velocities vp
"""
function update_velocities!(p, e, coeffs, dt)

    p.vp .+= coeffs * e .* p.qm .* dt

end
