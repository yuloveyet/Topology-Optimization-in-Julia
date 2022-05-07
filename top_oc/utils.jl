# export loads!, bcs!, passive_domain!
export prj, deta, dprj, cnt
export meshgrid
export Visualization

function prj(v::Array{Float64}, eta::Float64, beta::Int)
    return (tanh(beta * eta) .+ tanh.(beta * (v[:] .- eta))) ./ (tanh(beta * eta) + tanh(beta * (1 - eta)))
end

function deta(v::Array{Float64}, eta::Float64, beta::Int)
    return -beta * csch(beta) .* sech.(beta * (v[:] .- eta)) .^ 2 .* sinh.(v[:] * beta) .* sinh.((1 .- v[:]) * beta)
end

function dprj(v::Array{Float64}, eta::Float64, beta::Int)
    return beta * (1 .- tanh.(beta * (v .- eta)) .^ 2) ./ (tanh(beta * eta) + tanh(beta * (1 - eta)))
end

function cnt(v::Int, vCnt::NTuple{4,Int}, l::Int, ch::Float64, maxchang::Float64)
    return v + (l >= vCnt[1]) * (v < vCnt[2]) * (mod(l, vCnt[3]) == 0 || ch <= maxchang) * vCnt[4]
end


"""
meshgrid(vx)
Computes an (x,y)-grid from the vectors (vx,vx).
For more information, see the MATLAB documentation.
"""
meshgrid(v::AbstractVector) = meshgrid(v, v)

"""
meshgrid(vx,vy)
Computes an (x,y)-grid from the vectors (vx,vy).
For more information, see the MATLAB documentation.
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repeat(vx, m, 1), repeat(vy, 1, n))
end

"""
meshgrid(vx,vy,vz)
Computes an (x,y,z)-grid from the vectors (vx,vy,vz).
For more information, see the MATLAB documentation.
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
    vz::AbstractVector{T}) where {T}
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

function Visualization(setup::SetUp, x::Array{Float64}, loop::Int)
    nelx, nely = setup.nelx, setup.nely
    # cmap = cgrad(:Blues_9, rev=false)
    plot = heatmap(reshape(x, nely, nelx), c=:Blues_9, aspect_ratio=:equal, yflip=true, grid=false, axis=:off, tick=false, colorbar=false, border=nothing, dpi=300, size=(400, nely / nelx * 400), legend=:none, display_type=:gui)
    # display(plot)
    # if mod(loop, 10) == 0
    savefig(plot, "./top_mma/res_$loop.pdf")
    # end

    # PLOT FINAL DESIGN
    # heatmap(1.0 .- x[end:-1:1, :], yaxis=false, xaxis=false, legend=:none,color=:greys, grid=false, border=nothing, aspect_ratio=:equal)
    nothing
end