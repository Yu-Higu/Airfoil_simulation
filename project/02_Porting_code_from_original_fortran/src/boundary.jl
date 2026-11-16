function _write_velocity_field(d::CFDData, outdir::AbstractString)
    nstep = lpad(string(mod(d.NUM, 10000)), 4, '0')
    path = joinpath(outdir, "puv_field$(nstep).dat")
    open(path, "w") do io
        for i in 1:IG
            for j in 1:JG
                ubar = (d.U[i,j,1]+d.U[i-1,j,1])*0.5
                vbar = (d.U[i,j,2]+d.U[i,j-1,2])*0.5
                betav = (d.BETA[i,j,2]+d.BETA[i,j-1,2])*0.5
                betau = (d.BETA[i,j,1]+d.BETA[i-1,j,1])*0.5
                vort = -((d.U[i,j,1]+d.U[i-1,j,1])*0.5 -
                    (d.U[i,j-1,1]+d.U[i-1,j-1,1])*0.5)/d.DY[j] +
                    ((d.U[i,j,2]+d.U[i,j-1,2])*0.5 -
                    (d.U[i-1,j,2]+d.U[i-1,j-1,2])*0.5)/d.DX[i]
                @printf(io, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n",
                    d.XX[i], d.YY[j], d.P[i,j], ubar, vbar, betav, betau, vort)
            end
        end
    end
end

function UBC!(d::CFDData; outdir::AbstractString="TEXTDATA")
    # VERTICAL DIRECTION B.C.
    for i in 1:IG
        d.U[i,0,1] = d.U[i,1,1]
        d.U[i,JG+1,1] = d.U[i,JG,1]
        d.U[i,0,2] = 0.0
        d.U[i,JG,2] = 0.0
    end
    # STREAMWISE DIRECTION B.C.
    for j in 1:JG
        d.U[0,j,1] = 1.0
        d.U[IG+1,j,1] = d.U[IG,j,1]
        d.U[0,j,2] = -d.U[1,j,2]
        d.U[IG+1,j,2] = -d.U[IG,j,2]
    end
    # PRESSURE B.C.
    for i in 0:IG
        d.P[i,0] = d.P[i,1]
        d.P[i,JG+1] = d.P[i,JG]
    end
    for j in 1:JG
        d.P[0,j] = d.P[1,j]
        d.P[IG+1,j] = d.P[IG,j]
    end

    divmax = 0.0
    for j in 1:JG
        for i in 1:IG
            qqqq = (d.U[i,j,1]-d.U[i-1,j,1])*d.PDX[i] +
                   (d.U[i,j,2]-d.U[i,j-1,2])*d.PDY[j]
            divmax = max(divmax, abs(qqqq))
        end
    end
    cfl = 0.0
    for j in 1:JG
        for i in 1:IG
            cfll = (abs(d.U[i,j,1]+d.U[i-1,j,1])*d.PDX[i] +
                    abs(d.U[i,j,2]+d.U[i,j-1,2])*d.PDY[j])*0.5*d.DT
            cfl = max(cfl, cfll)
        end
    end
    usum = 0.0
    for j in 1:JG
        for i in 1:IG
            usum += d.U[i,j,1]*d.DY[j]
        end
    end
    usum = usum/Float64(IG)/d.YL

    if mod(d.NUM, 5) == 0
        _write_velocity_field(d, outdir)
    end

    p_sum = 0.0
    for i in 1:IG
        p_top = 0.0; p_btm = 0.0
        for j in 2:JG
            if d.BETA[i,j,2]-d.BETA[i,j-1,2] >= 0.50
                p_top = d.P[i,j]
            end
            if d.BETA[i,j,2]-d.BETA[i,j-1,2] <= -0.50
                p_btm = d.P[i,j-1]
            end
        end
        p_sum += d.DX[i]*(p_btm-p_top)
    end
    cl = p_sum*2.0

    p_sum = 0.0
    for j in 1:JG
        p_top = 0.0; p_btm = 0.0
        for i in 2:IG
            if d.BETA[i,j,1]-d.BETA[i-1,j,1] >= 0.5
                p_top = d.P[i+1,j]
            end
            if d.BETA[i,j,1]-d.BETA[i-1,j,1] <= -0.5
                p_btm = d.P[i-2,j]
            end
        end
        p_sum += d.DY[j]*(p_btm-p_top)
    end
    cd = p_sum*2.0
    return divmax, cfl, usum, cl, cd
end
