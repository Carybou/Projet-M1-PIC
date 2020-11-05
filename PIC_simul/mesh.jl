export Mesh

"""
    Mesh( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [0:nx]
- `ny` : indices are in [0:ny]
- `dimx` x `dimy`: mesh area
- `x, y` : node positions
- `dx, dy` : step size
"""
struct Mesh

    nx
    dimx
    x
    dx

    function Mesh(dimx, nx)

        x = LinRange(0, dimx, nx+1) |> collect

        dx = dimx / nx

        new( nx, dimx, x, dx)

    end

end
