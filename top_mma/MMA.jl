########################################################################################################
### MMA-Julia (1.7.1)         									
########################################################################################################
### MODIFIED BY YULI 2022/5/4                                                            
########################################################################################################
export mmasub
# Loading modules
using LinearAlgebra
using SparseArrays

########################################################################################################
### MMA FUNCTIONS                                                                                    ###
########################################################################################################

# Function for the MMA sub problem
function mmasub(m::Int, n::Int, iter::Int, xval::Array{Float64}, xmin::Array{Float64},
    xmax::Array{Float64}, xold1::Array{Float64}, xold2::Array{Float64}, f0val,
    df0dx::Array{Float64}, df0dx2:: Array{Float64},
    fval::Float64,        dfdx::Array{Float64}, dfdx2::Array{Float64},
    low::Array{Float64}, upp::Array{Float64}, a0::Float64,
    a::Array{Float64}, c::Array{Float64}, d::Array{Float64})

    # """
    # This function mmasub performs one MMA-iteration, aimed at solving the nonlinear programming problem: 
    #
    # Minimize    f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
    # subject to  f_i(x) - a_i*z - y_i < = 0,     i       = 1, ..., m
    # xmin_j <                           = x_j < = xmax_j, j  = 1,  ..., n
    # z >                                = 0,     y_i >   = 0, i   = 1,  ..., m
    # INPUT: 
    #
    # m     = The number of general constraints.
    # n     = The number of variables x_j.
    # iter  = Current iteration number ( =1 the first time mmasub is called).
    # xval  = Column vector with the current values of the variables x_j.
    # xmin  = Column vector with the lower bounds for the variables x_j.
    # xmax  = Column vector with the upper bounds for the variables x_j.
    # xold1 = xval, one iteration ago (provided that iter>1).
    # xold2 = xval, two iterations ago (provided that iter>2).
    # f0val = The value of the objective function f_0 at xval.
    # df0dx = Column vector with the derivatives of the objective function
    #             f_0 with respect to the variables x_j, calculated at xval.
    # fval = Column vector with the values of the constraint functions f_i, calculated at xval.
    # dfdx = (m x n)-matrix with the derivatives of the constraint functions
    #             f_i with respect to the variables x_j, calculated at xval.
    # dfdx(i,j) = the derivative of f_i with respect to x_j.
    # low       = Column vector with the lower asymptotes from the previous
    #             iteration (provided that iter>1).
    # upp = Column vector with the upper asymptotes from the previous
    #             iteration (provided that iter>1).
    # a0 = The constants a_0 in the term a_0*z.
    # a  = Column vector with the constants a_i in the terms a_i*z.
    # c  = Column vector with the constants c_i in the terms c_i*y_i.
    # d  = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    #
    # OUTPUT: 
    #
    # xmma = Column vector with the optimal values of the variables x_j
    #             in the current MMA subproblem.
    # ymma = Column vector with the optimal values of the variables y_i
    #             in the current MMA subproblem.
    # zmma = Scalar with the optimal value of the variable z
    #             in the current MMA subproblem.
    # lam = Lagrange multipliers for the m general MMA constraints.
    # xsi = Lagrange multipliers for the n constraints alfa_j - x_j < = 0.
    # eta = Lagrange multipliers for the n constraints x_j - beta_j < = 0.
    # mu  = Lagrange multipliers for the m constraints -y_i <         = 0.
    # zet = Lagrange multiplier for the single constraint -z <        = 0.
    # s   = Slack variables for the m general MMA constraints.
    # low = Column vector with the lower asymptotes,                   calculated and used
    #             in the current MMA subproblem.
    # upp = Column vector with the upper asymptotes, calculated and used
    #             in the current MMA subproblem.
    # """

    epsimin = sqrt(m + n) * 10^(-9)
    feps    = 0.000001
    asyinit = 0.5
    asyincr = 1.05
    asydecr = 0.65
    albefa  = 0.1
    een     = ones(n, 1)
    zeron   = zeros(n, 1)


    if iter <= 2
       low    = xval - asyinit .* (xmax - xmin)
       upp    = xval + asyinit .* (xmax - xmin)
    else
               zzz                    = (xval - xold1) .* (xold1 - xold2)
               factor                 = copy(een)
        factor[findall(zzz .> 0.0)]  .= asyincr
        factor[findall(zzz .< 0.0)]  .= asydecr
               low                    = xval - factor .* (xold1 - low)
               upp                    = xval + factor .* (upp - xold1)
               lowmin                 = xval - 10.0 .* (xmax - xmin)
               lowmax                 = xval - 0.01 .* (xmax - xmin)
               uppmin                 = xval + 0.01 .* (xmax - xmin)
               uppmax                 = xval + 10.0 .* (xmax - xmin)
               low                    = max.(low, lowmin)
               low                    = min.(low, lowmax)
               upp                    = min.(upp, uppmax)
               upp                    = max.(upp, uppmin)
    end

      zzz  = low + albefa .* (xval - low)
      alfa = max.(zzz, xmin)

      zzz  = upp - albefa .* (upp - xval)
      beta = min.(zzz, xmax)


      ux1    = upp - xval
      ux2    = ux1 .* ux1
      ux3    = ux2 .* ux1
      xl1    = xval - low
      xl2    = xl1 .* xl1
      xl3    = xl2 .* xl1
      ul1    = upp - low
      ulinv1 = een ./ ul1

      uxinv1 = een ./ ux1
      xlinv1 = een ./ xl1
      uxinv3 = een ./ ux3
      xlinv3 = een ./ xl3
      diap   = (ux3 .* xl1) ./ (2 * ul1)
      diaq   = (ux1 .* xl3) ./ (2 * ul1)

            p0                     = copy(zeron)
            q0                     = copy(zeron)
    p0[findall(df0dx .> 0.0)] = df0dx[findall(df0dx .> 0.0)]
            p0                     = p0 + 0.001 .* abs.(df0dx) + feps * ulinv1
            p0                     = p0 .* ux2
    #       q0                     = zeron
    q0[findall(df0dx .< 0.0)] = -df0dx[findall(df0dx .< 0.0)]
            q0                     = q0 + 0.001 * abs.(df0dx) + feps * ulinv1
            q0                     = q0 .* xl2
            dg0dx2                 = 2 * (p0 ./ ux3 + q0 ./ xl3)
            del0                   = df0dx2 - dg0dx2
            delpos0                = copy(zeron)
    delpos0[findall(del0 .> 0.0)]  = del0[findall(del0 .> 0.0)]
            p0                     = p0 + delpos0 .* diap
            q0                     = q0 + delpos0 .* diaq
            P                      = spzeros(m, n)
    P[findall(dfdx .> 0.0)]  = dfdx[findall(dfdx .> 0.0)]
    #       P                      = P * diag(ux2);
    #       P                      = P * spdiagm(n, n, 0 => ux2)
            P                      = P * sparse(Diagonal(ux2[:]))
            Q                      = spzeros(m, n)
    Q[findall(dfdx .< 0.0)]  = -dfdx[findall(dfdx .< 0.0)]
    #       Q                      = Q * diag(xl2);
    #       Q                      = Q * spdiagm(n, n, 0 => xl2)
            Q                      = Q * sparse(Diagonal(xl2[:]))
    #       dgdx2                  = 2.0*(P*diag(uxinv3) + Q*diag(xlinv3));
    #       dgdx2                  = P * spdiagm(n, n, 0 => uxinv3) + Q * spdiagm(n, n, 0 => xlinv3)
            dgdx2                  = P * sparse(Diagonal(uxinv3[:])) + Q * sparse(Diagonal(xlinv3[:]))
            dgdx2                  = 2.0 * dgdx2
            del                    = dfdx2 - dgdx2
            delpos                 = zeros(m, n)
     delpos[findall(del .> 0.0)]   = del[findall(del .> 0.0)]
    #       P                      = P + delpos*diag(diap);
    #       P                      = P + delpos * spdiagm(n, n, 0 => diap)
            P                      = P + delpos * sparse(Diagonal(diap[:]))
    #       Q                      = Q + delpos*diag(diaq);
    #       Q                      = Q + delpos * spdiagm(n, n, 0 => diap)
            Q                      = Q + delpos * sparse(Diagonal(diap[:]))
            b                      = P * uxinv1 + Q * xlinv1 .- fval


    # b  = P * uxinv + Q * xlinv - fval
    # Solving the subproblem by a primal-dual Newton method
    xmma, ymma, zmma, lam, xsi, eta, mu, zet, s = 
        subsolv(m, n, epsimin, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d)
    # Return values
    return xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp
end


# Function for solving the subproblem (can be used for MMA and GCMMA)
function subsolv(m::Int, n::Int, epsimin::Float64, low::Array{Float64}, upp::Array{Float64},
    alfa::Array{Float64}, beta::Array{Float64}, p0::Array{Float64}, q0::Array{Float64},
    P::Array{Float64}, Q::Array{Float64}, a0::Float64,        a::Array{Float64}, b::Array{Float64},
    c::Array{Float64}, d::Array{Float64})

    # """
    # This function subsolv solves the MMA subproblem: 
    #
    # minimize SUM[p0j/(uppj-xj) + q0j/(xj-lowj)] + a0*z + SUM[ci*yi + 0.5*di*(yi)^2],
    #
    # subject to SUM[pij/(uppj-xj) + qij/(xj-lowj)] - ai*z - yi < = bi,
    # alfaj <                                                     = xj < = betaj, yi > = 0, z > = 0.
    #
    # Input : m,    n,    low,  upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
    # Output: xmma, ymma, zmma, slack variables and Lagrange multiplers.
    # """

    een     = ones(Float64, n)
    eem     = ones(Float64, m)
    epsi    = 1.0
    epsvecn = epsi .* een
    epsvecm = epsi .* eem
    x       = 0.5 .* (alfa + beta)
    y       = copy(eem)
    z       = 1.0
    lam     = copy(eem)
    xsi     = een ./ (x - alfa)
    xsi     = max.(xsi, een)
    eta     = een ./ (beta - x)
    eta     = max.(eta, een)
    mu      = max.(eem, 0.5 .* c)
    zet     = 1.0
    s       = copy(eem)
    itera   = 0

    while epsi > epsimin    # Start while epsi>epsimin
        epsvecn = epsi .* een
        epsvecm = epsi .* eem
        ux1     = upp - x
        xl1     = x - low
        ux2     = ux1 .* ux1
        xl2     = xl1 .* xl1
        uxinv1  = een ./ ux1
        xlinv1  = een ./ xl1

        plam = p0 + (P)' * lam
        qlam = q0 + (Q)' * lam
        gvec   = P * uxinv1 + Q * xlinv1
        dpsidx = (plam ./ ux2) - (qlam ./ xl2)

        rex = dpsidx - xsi + eta
        rey = c + d .* y - mu - lam
        rez = a0 .- zet .- (a)' * lam
        relam = gvec - a .* z - y + s - b
        rexsi = xsi .* (x - alfa) - epsvecn
        reeta = eta .* (beta - x) - epsvecn
        remu  = mu .* y - epsvecm
        rezet = zet * z - epsi
        res   = lam .* s - epsvecm

        residu1 = [rex' rey' rez]'
        residu2 = [relam' rexsi' reeta' remu' rezet res']'
        residu = [residu1' residu2']'
        residunorm = norm(residu, 2)
        residumax  = maximum(abs.(residu))

        ittt = 0
        while (residumax > 0.9 * epsi) && (ittt < 100) # Start while (residumax>0.9*epsi) and (ittt<100)
            ittt  = ittt + 1
            itera = itera + 1

            if ittt == 100
                println("max inner iter reached")
            end
            ux1    = upp - x
            xl1    = x - low
            ux2    = ux1 .* ux1
            xl2    = xl1 .* xl1
            ux3    = ux1 .* ux2
            xl3    = xl1 .* xl2
            uxinv1 = een ./ ux1
            xlinv1 = een ./ xl1
            uxinv2 = een ./ ux2
            xlinv2 = een ./ xl2
            plam   = p0 + (P)' * lam
            qlam = q0 + (Q)' * lam
              gvec      = P * uxinv1 + Q * xlinv1
            # GG        = P .* transpose(uxinv2) - Q .* transpose(xlinv2)
            # GG        = P .* spdiagm(n,n, 0 => uxinv2) - Q .* spdiagm(n, n, 0 => xlinv2)
              GG        = P * sparse(Diagonal(uxinv2[:])) - Q * sparse(Diagonal(xlinv2[:]))
              dpsidx    = (plam ./ ux2) - (qlam ./ xl2)
              delx      = dpsidx - epsvecn ./ (x - alfa) + epsvecn ./ (beta - x)
              dely      = c + d .* y - lam - epsvecm ./ y
              delz      = a0 .- transpose(a) * lam .- epsi / z
              dellam    = gvec - a .* z - y - b + epsvecm ./ lam
              diagx     = plam ./ ux3 + qlam ./ xl3
              diagx     = 2.0 .* diagx + xsi ./ (x - alfa) + eta ./ (beta - x)
              diagxinv  = een ./ diagx
              diagy     = d + mu ./ y
              diagyinv  = eem ./ diagy
              diaglam   = s ./ lam
              diaglamyi = diaglam + diagyinv

            if m < n # Start if m < n
                blam = dellam .+ dely ./ diagy .- GG * (delx ./ diagx)
                bb   = [blam' delz]'
                # Alam = spdiagm(0 => diaglamyi) + (GG .* transpose(diagxinv)) * transpose(GG)
                # Alam = spdiagm(m, m, 0 => diaglamyi) + (GG .* spdiagm(n, n, 0 => diagxinv)) * transpose(GG)
                Alam = sparse(Diagonal(diaglamyi[:])) + (GG * sparse(Diagonal(diagxinv[:]))) * (GG)'
                AA    = [Alam a
                            a' -zet/z]
                solut = AA \ bb
                dlam  = solut[1:m]
                dz    = solut[m+1]
                dx    = -delx ./ diagx - ((GG)' * dlam) ./ diagx
            else
                diaglamyiinv = eem ./ diaglamyi
                dellamyi     = dellam + dely ./ diagy
                # Axx          = spdiagm(0 => diagx) + (transpose(GG) .* transpose(diaglamyiinv)) * GG
                # Axx          = spdiagm(n,n, 0 => diagx) + (transpose(GG) .* spdiagm(m, m, 0 => diaglamyiinv)) * GG
                Axx          = sparse(Diagonal(diagx[:])) + ((GG)' .* sparse(Diagonal(diaglamyiinv[:]))) * GG
                azz          = zet / z + (a)' * (a ./ diaglamyi)
                axz = -(GG)' * (a ./ diaglamyi)
                bx = delx + (GG)' * (dellamyi ./ diaglamyi)
                bz = delz - (a)' * (dellamyi ./ diaglamyi)
                # AAr1 = [Axx axz]
                # AAr2 = [transpose(axz) azz]
                # AA   = [AAr1; AAr2]
                AA   = [Axx axz
                        axz' azz]
                bb = [-bx' -bz]'
                solut = AA \ bb
                dx    = solut[1:n]
                dz    = solut[n+1]
                dlam  = (GG * dx) ./ diaglamyi - dz .* (a ./ diaglamyi) + dellamyi ./ diaglamyi
            end # End if m<n

            dy   = -dely ./ diagy + dlam ./ diagy
            dxsi = -xsi .+ epsvecn ./ (x - alfa) - (xsi .* dx) ./ (x - alfa)
            deta = -eta .+ epsvecn ./ (beta - x) + (eta .* dx) ./ (beta - x)
            dmu  = -mu .+ epsvecm ./ y - (mu .* dy) ./ y
            dzet = -zet .+ epsi / z - zet * dz / z
            ds   = -s .+ epsvecm ./ lam - (s .* dlam) ./ lam
            xx   = [y' z lam' xsi' eta' mu' zet s']'
            dxx  = [dy' dz dlam' dxsi' deta' dmu' dzet ds']'
            #
            stepxx    = -1.01 .* dxx ./ xx
            stmxx     = maximum(stepxx)
            stepalfa  = -1.01 .* dx ./ (x - alfa)
            stmalfa   = maximum(stepalfa)
            stepbeta  = 1.01 .* dx ./ (beta - x)
            stmbeta   = maximum(stepbeta)
            stmalbe   = max.(stmalfa, stmbeta)
            stmalbexx = max.(stmalbe, stmxx)
            stminv    = max.(stmalbexx, 1.0)
            steg      = 1.0 / stminv
            #
            xold   = copy(x)
            yold   = copy(y)
            zold   = copy(z)
            lamold = copy(lam)
            xsiold = copy(xsi)
            etaold = copy(eta)
            muold  = copy(mu)
            zetold = copy(zet)
            sold   = copy(s)
            #
            itto    = 0
            resinew = 2.0 * residunorm
            # Start: while (resinew>residunorm) and (itto<50)
            while (resinew > residunorm) && (itto < 50)
                itto = itto + 1

                x      = xold + steg .* dx
                y      = yold + steg .* dy
                z      = zold + steg * dz
                lam    = lamold + steg .* dlam
                xsi    = xsiold + steg .* dxsi
                eta    = etaold + steg .* deta
                mu     = muold + steg .* dmu
                zet    = zetold + steg * dzet
                s      = sold + steg .* ds
                ux1    = upp - x
                xl1    = x - low
                ux2    = ux1 .* ux1
                xl2    = xl1 .* xl1
                uxinv1 = een ./ ux1
                xlinv1 = een ./ xl1
                plam   = p0 + (P)' * lam
                qlam   = q0 + (Q)' * lam
                gvec   = P * uxinv1 + Q * xlinv1
                dpsidx = plam ./ ux2 - qlam ./ xl2

                rex   = dpsidx - xsi + eta
                rey   = c .+ d .* y - mu - lam
                rez   = a0 - zet .- (a)' * lam
                relam = gvec - a .* z - y + s - b
                rexsi = xsi .* (x - alfa) - epsvecn
                reeta = eta .* (beta - x) - epsvecn
                remu  = mu .* y - epsvecm
                rezet = zet * z - epsi
                res   = lam .* s - epsvecm

                residu1 = [rex' rey' rez]'
                residu2 = [relam' rexsi' reeta' remu' rezet res']'
                residu = [residu1' residu2']'
                resinew = norm(residu, 2)
                steg    = steg / 2.0
            end # End: while (resinew>residunorm) and (itto<50)
            residunorm = copy(resinew)
            residumax  = maximum(abs.(residu))
            steg       = 2.0 * steg
        end # End: while (residumax>0.9*epsi) and (ittt<200)
        epsi = 0.1 * epsi
    end # End: while epsi>epsimin
    xmma   = copy(x)
    ymma   = copy(y)
    zmma   = copy(z)
    lamma  = lam
    xsimma = xsi
    etamma = eta
    mumma  = mu
    zetmma = zet
    smma   = s
    # Return values
    return xmma, ymma, zmma, lamma, xsimma, etamma, mumma, zetmma, smma
end

