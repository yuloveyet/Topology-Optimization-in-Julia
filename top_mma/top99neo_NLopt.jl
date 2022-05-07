# __precompile__()
# module TopOpt99neo_MMA

using LinearAlgebra, SparseArrays
using Plots
      #             : heatmap, savefig, @animate
using ImageFiltering: imfilter
using Statistics: mean

using BenchmarkTools

using NLopt

BLAS.set_num_threads(1)

# abstract type top99neo end

include("utils.jl")
include("MMA.jl")

export SetUp, Mat, DiscretizationFeature, LoadsSupportsBCs, Initialization, Filter
export Optimization, Visualization

mutable struct SetUp
    nelx:: Int
    nely:: Int

    eta:: Float64
    beta:: Int
    penalCnt:: Array
    betaCnt:: Array

    move:: Float64
    maxit:: Int
    maxchang:: Float64

    pasS:: Array{Int}
    pasV:: Array{Int}

    function SetUp()
          nelx                                = 120
          nely                                = 40
        # volfrac                             = 0.3
          maxit                               = 500
          move                                = 0.1
          beta                                = 2
          eta                                 = 0.5
          maxchang                            = 1e-6
          penalCnt                            = [1, 3, 20, 0.2]
          betaCnt                             = [1, 32, 20, 2]
        # continuation scheme on beta parCont = { istart,maxPar,steps,deltaPar }

        #    elNrs = reshape(1:nely*nelx, nely, nelx)
        #    a1    = elNrs[Int.(nely/4:nely/2), Int.(nelx/4:nelx/2)]
        #    pasS,  pasV       = Array([]), a1[:]
        pasS, pasV = Array([]), Array([])

        new(nelx, nely, eta, beta, penalCnt, betaCnt, move, maxit, maxchang, pasS, pasV)
    end

end

mutable struct Mat
    E0:: Float64
    Emin:: Float64
    ν:: Float64
    penal:: Float64
    volfrac:: Float64

    function Mat()
        E0      = 1.0
        Emin    = 1e-9 * E0
        ν       = 0.3
        penal   = 3.0
        volfrac = 0.3
        new(E0, Emin, ν, penal, volfrac)
    end
end

struct DiscretizationFeature
    nEl:: Int
    nodeNrs:: Array{Int,2}
    nDof:: Int
    act

    cMat:: Array{Int}
    Iar:: Array{Int}
    Ke:: Array{Float64,1}
    Ke0:: Array{Float64,2}

    function DiscretizationFeature(setup::SetUp, mat::Mat)
             nelx  = setup.nelx
             nely  = setup.nely
        pasS, pasV = setup.pasS, setup.pasV
             ν     = mat.ν

        nEl     = nely * nelx
        nodeNrs = reshape(1:(1+nelx)*(1+nely), nely + 1, nelx + 1)
        nDof    = (nely + 1) * (nelx + 1) * 2
        act     = setdiff(collect(1:nEl), union(pasS, pasV))

        cVec = reshape(2 * nodeNrs[1:end-1, 1:end-1] .+ 1, nEl, 1)
        cMat = Int.(repeat(cVec, 1, 8) + repeat([0 1 2 * nely .+ [2 3 0 1] -2 -1], nelx * nely, 1))

            FuckRow    = [1 2 3 4 5 6 7 8]
        sI::Array{Int}, sII::Array{Int} = copy(FuckRow), fill(1, 1, 8)
        @inbounds for j in 2: 8
              sI  = cat(sI, FuckRow[j:8]'; dims=2)
            # sII = cat(2, sII, repmat(j, 1, 8 - j + 1))
              sII = cat(sII, fill(j, 1, 8 - j + 1); dims=2)
        end
        iK::Array{Int,2}, jK::Array{Int,2} = cMat[:, sI][:, 1, :]', cMat[:, sII][:, 1, :]'
            Iar          = sort([iK[:] jK[:]]; dims=2, rev=true) # comma is a newline
        # iK[:], jK[:] .= 0.0, 0.0
            c1           = [12, 3, -6, -3, -6, -3, 0, 3, 12, 3, 0, -3, -6, -3, -6, 12, -3, 0, -3, -6, 3, 12, 3, -6, 3, -6, 12, 3, -6, -3, 12, 3, 0, 12, -3, 12]
            c2           = [-4, 3, -2, 9, 2, -3, 4, -9, -4, -9, 4, -3, 2, 9, -2, -4, -3, 4, 9, 2, 3, -4, -9, -2, 3, 2, -4, 3, -2, 9, -4, -9, 4, -4, -3, -4]
            Ke           = 1 / (1 - ν^2) / 24 .* (c1 .+ ν .* c2) # half-KE vector
        # full KE
        # Ke0::Array{Float64} = zeros(8, 8)
        # start_id,            end_id = 1, 8
        # for i in 1: 8
        # Ke0[i:8, i] = Ke[start_id:end_id]
        # start_id,    end_id = end_id + 1, 2 * end_id - start_id
        # end

        Ke0::Array{Float64}    = zeros(8, 8)
        #    Index::Array{Int} = [sI' sII']
        #    Ke0[sI, sII]      = Ke
             Index             = findall(isequal(1), tril(ones(8, 8)))
        Ke0[Index]            = Ke'
        # Ke0 = reshape(Ke0, 8, 8)
          Ke0 = Ke0 + Ke0' - diagm(diag(Ke0))  
        # diag return the diag of a matrix
        # diagm build a matirx with a diag
        new(nEl, nodeNrs, nDof, act, cMat, Iar, Ke, Ke0)
    end
end

mutable struct LoadsSupportsBCs
    lcDof:: Array{Int}
    F
    free:: Array{Int}
    fixed:: Array{Int}

    function LoadsSupportsBCs(setup::SetUp, disfeature::DiscretizationFeature)
        nelx    = setup.nelx
        nely    = setup.nely
        nodeNrs = disfeature.nodeNrs
        nDof    = disfeature.nDof

        load_position::Symbol          = :half_MBB
        if             load_position  == :half_MBB
        load_nodey    , load_nodex     = 1, 1
        #              fixed           = union([1:2:2*(nely+1)], 2 * nodeNrs[nely+1, nelx+1])
                       fixed           = union(collect(1:2:2*(nely+1)), 2 * nodeNrs[end, end])
        elseif         load_position  == :cantilever
                       load_nodey      = nely + 1
                       load_nodex      = nelx / 2 + 1
                       fixed           = 1:2*(nely+1)
        end

        F = spzeros(nDof, 1)

        load_type::Symbol                          = :pin
        if         load_type == :pin # 1 point
                   lcDof                           = collect(2 * nodeNrs[load_nodey, load_nodex])
        F[2, 1]                           = -1.0
        elseif     load_type == :points # 5 points
        #          lcDof                           = collect(2 * nodeNrs[load_nodey, load_nodex], nodeNrs[load_nodey, load_nodex-1], nodeNrs[load_nodey, load_nodex-2], nodeNrs[load_nodey, load_nodex+1], nodeNrs[load_nodey, load_nodex+2])
        F[lcDof', ones(length(lcDof'))] .= -1.0
        elseif     load_type == :line
                   lcDof                           = [2:2*(nely+1):nDof]
                   F                               = spzeros(nDof, 1)
        F[lcDof', ones(length(lcDof'))] .= -1.0
        #          F                               = sparse(lcDof', ones(length(lcDof')), -1.0)
        end
        all  = collect(1:nDof)
        free = setdiff(all, fixed)

        new(lcDof, F, free, fixed)
    end
end

struct Initialization
    x:: Array{Float64}
    xPhys:: Array{Float64}
    xOld:: Array{Float64}
    ch:: Float64
    loop:: Int

    function Initialization(setup::SetUp, disfeature::DiscretizationFeature, mat::Mat)
          pasV    = setup.pasV
          pasS    = setup.pasS
          act     = disfeature.act
          volfrac = mat.volfrac
          nEl     = disfeature.nEl
        # nDof    = disfeature.nDof
        # column vectors
              x       = zeros(nEl, 1)
        x[act] .= volfrac
        x[pasS] .= 1.0
        xPhys, xOld,   ch, loop = copy(x), ones(nEl, 1), 1.0, 0
        # x̅  x̃

        new(x, xPhys, xOld, ch, loop)
    end
end

mutable struct Filter
    rmin:: Float64
    ft:: Int
    h:: Array{Float64}
    Hs:: Array{Float64}
    dHs:: Array{Float64}

    function Filter(setup::SetUp)
           nelx    = setup.nelx
           nely    = setup.nely
        #  volfrac = mat.volfrac
        #  bcF     = setup.bcF
           rmin    = 6.5
           ft      = 3 # 1 density 2 projection 3 volume preserving
        dy, dx     = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1)
           h       = max.(0, rmin .- sqrt.(dx .^ 2 + dy .^ 2))
           Hs      = imfilter(ones(nely, nelx), h, "symmetric")
           dHs     = Hs

        new(rmin, ft, h, Hs, dHs)
    end
end


function FiniteElementAnalasys(mat::Mat, disfeature::DiscretizationFeature, load::LoadsSupportsBCs, xPhys::Array{Float64})
    nEl , nDof,  Iar,     Ke        = disfeature.nEl, disfeature.nDof, disfeature.Iar, disfeature.Ke
         act    = disfeature.act
    Emin, penal, E0      = mat.Emin, mat.penal,       mat.E0
    F   , free  = load.F, load.free
    # Interplotation model
      sK      = Emin .+ xPhys .^ penal .* (E0 - Emin)
      sK      = reshape(Ke[:] * sK', length(Ke) * nEl, 1)
      K       = sparse(Iar[:, 1], Iar[:, 2], vec(sK), nDof, nDof)
      U       = zeros(nDof, 1)
    # U[free] = cholesky(Symmetric(K[free, free], :L), check=false) \ F[free
    U[free] = cholesky(Symmetric(K, :L)[free, free], check=false) \ F[free]
    Obj   = F' * U
    Vf = sum(xPhys[act])/length(act)
    return U, Obj, Vf
end

function SensitivityAnalasys(setup::SetUp, filter::Filter, mat::Mat, disfeature::DiscretizationFeature, U::Array{Float64}, xPhys::Array{Float64})
    nelx, nely, act         = setup.nelx,     setup.nely,      disfeature.act
    dHs , h    = filter.dHs, filter.h
    E0  , Emin, penal       = mat.E0,         mat.Emin,  mat.penal
    nEl , cMat, Ke0         = disfeature.nEl, disfeature.cMat, disfeature.Ke0
         act   = disfeature.act

    dsK, dV    = zeros(nEl, 1), zeros(nEl, 1)
    dV[act]   .= 1.0/length(act)
    dsK[act]   = -penal * (E0 - Emin) .* xPhys[act] .^ (penal - 1)

    dc  = dsK .* sum((U[cMat] * Ke0) .* U[cMat], dims=2)
    dc  = imfilter(reshape(dc, nely, nelx) ./ dHs, h, "symmetric")
    dV0 = imfilter(reshape(dV, nely, nelx) ./ dHs, h, "symmetric")

    return reshape(dc, nEl, 1), reshape(dV0, nEl, 1)
end

function MMAupdate(act, xval::Array{Float64}, low::Array{Float64}, upp::Array{Float64}, 
        xold1:: Array{Float64}, xold2:: Array{Float64},
    Obj, Vf:: Float64,
        dc:: Array{Float64}, dV0:: Array{Float64}, loop:: Int, volfrac:: Float64)
    move = 0.1
    ### Initiation of MMA ###
    m = 1
    # active design variable
    n     = length(act)
    onen  = ones(n, 1)
    onem  = ones(m, 1)
    zeron = zeros(n, 1)
    zerom = zeros(m, 1)
    a_mma = zerom
    c_mma = 1.0e3 * onem
    d_mma = zerom
    a0    = 1.0
    # column vector
    # xval = xval
    xmin = max.(xval .- move, zeron)
    xmax = min.(xval .+ move, onen)

    # low = low
    # upp = upp
    # objective function   
    f0val  = Obj
    df0dx  = dc[act]
    df0dx2 = 0.0 * df0dx
    # constraint function
    fval  = Vf / volfrac - 1.0  # column vector
    dfdx  = reshape(dV0[act], 1, length(act)) ./ volfrac   # (m * n)
    dfdx2 = 0.0 * dfdx

    # The MMA subproblem is solved at the point xval: 
    xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp = 
        mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2,
            f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, a0, a_mma, c_mma, d_mma)
    return xmma, low, upp
end

function Visualization(setup::SetUp, x::Array{Float64}, loop::Int)
    nelx, nely = setup.nelx, setup.nely
    #    cmap  = cgrad(:Blues_9, rev=false)
         plot  = heatmap(reshape(x, nely, nelx), c=:Blues_9, aspect_ratio=:equal, yflip=true, grid=false, axis=:off, tick=false, colorbar=false, border=nothing, dpi=300, size=(400, nely / nelx * 400), legend=:none, display_type=:gui)
    # display(plot)
    # if mod(loop, 10) = = 0
    # savefig(plot, "./res/res_$loop.pdf")
    savefig(plot, "./res/design_$loop.png")
    # end

    # PLOT FINAL DESIGN
    # heatmap(1.0 .- x[end:-1:1, :], yaxis=false, xaxis=false, legend=:none,color=:greys, grid=false, border=nothing, aspect_ratio=:equal)
    nothing
end


function Optimization(setup::SetUp, mat::Mat, load::LoadsSupportsBCs, filter::Filter, ini::Initialization, disfeature::DiscretizationFeature)
    ch, loop  = ini.ch, ini.loop
    x , xPhys, xOld    = ini.x, ini.xPhys, ini.xOld

    maxit,       maxchang = setup.maxit, setup.maxchang
    nely,       nelx     = setup.nely,  setup.nelx
    eta,       beta,     penalCnt,     betaCnt   = setup.eta, setup.beta, setup.penalCnt, setup.betaCnt
    volfrac, penal, nEl      = mat.volfrac, mat.penal, disfeature.nEl
    act    = disfeature.act

    Hs, h, ft = filter.Hs, filter.h, filter.ft

    opt_hist = []
    vf_hist  = []

    ### Initiation of MMA ###
    xval  = copy(x[act])
    xold1 = copy(xval)
    xold2 = copy(xold1)

    low = copy(xval)
    upp = copy(low)
    #########################
    loop_beta = 0
    anim = @animate while ch > maxchang && loop < maxit || beta < betaCnt[2]
        @time begin
            loop = loop + 1
            # COMPUTE PHYSICAL DENSITY FIELD 
                  xTilde = imfilter(reshape(x, nely, nelx), h, "symmetric") ./ Hs
            xPhys[act]   = copy(xTilde[act])
            if ft > 1
                f = (mean(prj(xPhys[act], eta, beta)) .- volfrac) * (ft == 3)
                while abs(f) > maxchang
                    eta = eta - f / mean(deta(xPhys[:], eta, beta))
                    f   = mean(prj(xPhys[act], eta, beta)) - volfrac
                end
                filter.dHs = Hs ./ reshape(dprj(xTilde, eta, beta), nely, nelx)
                xPhys      = prj(xPhys, eta, beta)
            end
            ch   = norm(xPhys - xOld) ./ sqrt(nEl)
            xOld = copy(xPhys)
            #~ SETUP AND SOLVE EQUILIBRIUM EQUATIONS
            U, C, Vf = FiniteElementAnalasys(mat, disfeature, load, xPhys)
            push!(opt_hist, C)
            push!(vf_hist, Vf)
            #~ COMPUTE SENSITIVITIES
            dc, dV0 = SensitivityAnalasys(setup, filter, mat, disfeature, U, xPhys)
            #~ MMA iteration
            xmma, low, upp = MMAupdate(act, xval, low, upp, xold1, xold2, C, Vf, dc, dV0, loop, Mat().volfrac)
            # Some vectors are updated: 
            xold2 = copy(xold1)
            xold1 = copy(xval)
            xval  = copy(xmma)

            x[act] = xval

            #~ CONTINUATION
            # if  mat.penal == 1.0
                mat.penal = cnt(mat.penal, penalCnt,loop, ch, maxchang)
            # end
            if mat.penal < 3.0
                beta      = 2.0
            else
                loop_beta  = loop_beta + 1
                beta       = cnt(beta, betaCnt, loop_beta, ch, maxchang)
                # restart mma
            end

            heatmap(reshape(xPhys, nely, nelx), c=:Blues_9, aspect_ratio=:equal, yflip=true, grid=false, axis=:off, tick=false, colorbar=false, border=nothing, dpi=300, size=(400, nely / nelx * 400), legend=:none)

        end
        if mod(loop, 20) == 0
            println("It.: $loop C.: $C Vf.: $Vf ch.: $ch, p.: $(mat.penal) beta.:$beta eta.: $eta ")
            Visualization(setup, xPhys, loop)
        end
    end
    # gif(anim, "./top/cantilever.gif", fps=8)
    return xPhys, opt_hist, vf_hist, anim
end

# end
