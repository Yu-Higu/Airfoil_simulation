# FS0 / FS1: temporary velocity field (原Fortranのループ順序を維持)
function _velocity_predictor!(d::CFDData, first_step::Bool)
    pre = 1.0 / d.RE

    # STREAMWISE DIRECTION
    for j in 1:JG
        for i in 1:IG
            if d.BETA[i,j,1] != 1.0
                continue
            end
            ududx = 0.5 * ((d.U[i+1,j,1] + d.U[i,j,1]) *
                (d.U[i+1,j,1] - d.U[i,j,1]) * d.PDX[i+1] * d.DX[i] +
                (d.U[i,j,1] + d.U[i-1,j,1]) *
                (d.U[i,j,1] - d.U[i-1,j,1]) * d.PDX[i] * d.DX[i+1]) * d.PDX2[i]
            vdudy = 0.5 * ((d.U[i+1,j,2]*d.DX[i] + d.U[i,j,2]*d.DX[i+1]) *
                (d.U[i,j+1,1] - d.U[i,j,1]) * d.PDYY[j+1] +
                (d.U[i+1,j-1,2]*d.DX[i] + d.U[i,j-1,2]*d.DX[i+1]) *
                (d.U[i,j,1] - d.U[i,j-1,1]) * d.PDYY[j]) * d.PDX2[i]
            u_con = ududx + vdudy
            uvisx = (d.DX[i+1]*d.U[i-1,j,1] - d.DX2[i]*d.U[i,j,1] +
                d.DX[i]*d.U[i+1,j,1]) * 2.0*d.PDX3[i]*pre
            uvisy = (d.DYY[j+1]*d.U[i,j-1,1] - d.DYY2[j]*d.U[i,j,1] +
                d.DYY[j]*d.U[i,j+1,1]) * 2.0*d.PDYY3[j]*pre
            u_vis = uvisx + uvisy
            dpdx = (d.P[i+1,j] - d.P[i,j]) * d.PDXX[i+1]
            rhs = u_con - u_vis + dpdx
            if first_step
                d.UT[i,j,1] = d.U[i,j,1] - d.DT*rhs
            else
                d.UT[i,j,1] = (d.U[i,j,1] - d.DT*(rhs*1.5 - d.CNX_OLD[i,j]*0.5))*d.BETA[i,j,1]
            end
            d.CNX_OLD[i,j] = rhs
        end
    end

    # VERTICAL DIRECTION
    for j in 1:JG-1
        for i in 1:IG
            if d.BETA[i,j,2] != 1.0
                continue
            end
            udvdx = 0.5 * ((d.U[i,j+1,1]*d.DY[j] + d.U[i,j,1]*d.DY[j+1]) *
                (d.U[i+1,j,2] - d.U[i,j,2])*d.PDXX[i+1] +
                (d.U[i-1,j+1,1]*d.DY[j] + d.U[i-1,j,1]*d.DY[j+1]) *
                (d.U[i,j,2] - d.U[i-1,j,2])*d.PDXX[i]) * d.PDY2[j]
            vdvdy = 0.5 * ((d.U[i,j+1,2] + d.U[i,j,2]) *
                (d.U[i,j+1,2] - d.U[i,j,2])*d.PDY[j+1]*d.DY[j] +
                (d.U[i,j,2] + d.U[i,j-1,2]) *
                (d.U[i,j,2] - d.U[i,j-1,2])*d.PDY[j]*d.DY[j+1]) * d.PDY2[j]
            v_con = udvdx + vdvdy
            vvisx = (d.DXX[i+1]*d.U[i-1,j,2] - d.DXX2[i]*d.U[i,j,2] +
                d.DXX[i]*d.U[i+1,j,2]) * 2.0*d.PDXX3[i]*pre
            vvisy = (d.DY[j+1]*d.U[i,j-1,2] - d.DY2[j]*d.U[i,j,2] +
                d.DY[j]*d.U[i,j+1,2]) * 2.0*d.PDY3[j]*pre
            v_vis = vvisx + vvisy
            dpdy = (d.P[i,j+1] - d.P[i,j]) * d.PDYY[j+1]
            rhs = v_con - v_vis + dpdy
            if first_step
                d.UT[i,j,2] = (d.U[i,j,2] - d.DT*rhs)*d.BETA[i,j,2]
            else
                d.UT[i,j,2] = (d.U[i,j,2] - d.DT*(rhs*1.5 - d.CNY_OLD[i,j]*0.5))*d.BETA[i,j,2]
            end
            d.CNY_OLD[i,j] = rhs
        end
    end

    # CUT-CELL VELOCITY CORRECTION
    for i in 1:IG
        for j in 1:JG
            if d.IDX[i,j] == 3
                d.UT[i,j,1] = d.UT[i,j+1,1]*d.GAMMA[i,j,1]
            end
            if d.IDX[i,j] == 5
                d.UT[i,j,1] = d.UT[i,j-1,1]*d.GAMMA[i,j,1]
            end
        end
    end
    return nothing
end

function FS0!(d::CFDData)
    # Fortran FS0 initializes U and P over 0:IG+1,0:JG+1 but BETA is
    # declared only over 0:IG,0:JG. Ghost-cell BETA is therefore treated
    # as the free-stream value (1) rather than performing an out-of-bounds read.
    for j in 0:JG+1
        for i in 0:IG+1
            d.U[i,j,1] = d.BETA[i,j,1]
            d.U[i,j,2] = 0.0
            d.P[i,j] = 0.0
        end
    end
    return _velocity_predictor!(d, true)
end
FS1!(d::CFDData) = _velocity_predictor!(d, false)

function FS2!(d::CFDData)
    igf = 6
    istart = div(IG*(igf-1), igf) + 1
    # FRINGE CHORDWISE VELOCITY, then SPANWISE VELOCITY
    for j in 1:JG
        for i in istart:IG
            d.UT[i,j,1] = (d.UT[i,j,1]-1.0)*(-tanh(4.5*(d.X[i]-d.XL))) + 1.0
        end
    end
    for j in 1:JG
        for i in istart:IG
            d.UT[i,j,2] = d.UT[i,j,2]*(-tanh(4.5*(d.X[i]-d.XL)))
        end
    end
    for i in 1:IG
        d.UT[i,0,1] = d.UT[i,1,1]
        d.UT[i,JG+1,1] = d.UT[i,JG,1]
        d.UT[i,0,2] = 0.0
        d.UT[i,JG,2] = 0.0
    end
    for j in 1:JG
        d.UT[0,j,1] = 1.0
        d.UT[IG+1,j,1] = d.UT[IG,j,1]
        d.UT[0,j,2] = -d.UT[1,j,2]
        d.UT[IG+1,j,2] = -d.UT[IG,j,2]
    end

    pdt = 1.0/d.DT
    omg = 1.4
    eps = Float64(Float32(1.0e-4))
    ncount = 0
    for n in 1:50000
        ncount = n
        ermax = 0.0
        for i in 1:IG
            for j in 1:JG
                phi_new = (d.D_W[i]*d.PHI[i-1,j] + d.D_E[i]*d.PHI[i+1,j] +
                    d.D_S[j]*d.PHI[i,j-1] + d.D_N[j]*d.PHI[i,j+1] -
                    pdt*((d.UT[i,j,1]-d.UT[i-1,j,1])*d.PDX[i] +
                         (d.UT[i,j,2]-d.UT[i,j-1,2])*d.PDY[j]))*d.PDPHI[i,j]
                phidif = phi_new - d.PHI[i,j]
                d.PHI[i,j] = d.PHI[i,j] + omg*phidif
                ermax = max(ermax, abs(phidif))
            end
        end
        for i in 1:IG
            d.PHI[i,0] = d.PHI[i,1]
            d.PHI[i,JG+1] = d.PHI[i,JG]
        end
        for j in 1:JG
            d.PHI[0,j] = d.PHI[1,j]
            d.PHI[IG+1,j] = d.PHI[IG,j]
        end
        if ermax < eps
            break
        end
    end
    for j in 1:JG
        for i in 1:IG
            d.U[i,j,1] = d.UT[i,j,1] - d.DT*(d.PHI[i+1,j]-d.PHI[i,j])*d.PDXX[i+1]
            d.U[i,j,2] = d.UT[i,j,2] - d.DT*(d.PHI[i,j+1]-d.PHI[i,j])*d.PDYY[j+1]
        end
    end
    p_stn = d.P[1,1] + d.PHI[1,1]
    for j in 1:JG
        for i in 1:IG
            d.P[i,j] = d.P[i,j] + d.PHI[i,j] - p_stn
        end
    end
    return ncount
end
