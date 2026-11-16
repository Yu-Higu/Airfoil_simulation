function _read_airfoil(path::AbstractString)
    lines = readlines(path)
    length(lines) >= IGS + 1 || error("翼座標ファイルの点数が不足しています: $path")
    wingname = rstrip(lines[1])
    xt = zeros(Float64, IGS); yt = similar(xt); yb = similar(xt)
    for k in 1:IGS
        line = lines[k+1]
        length(line) >= 27 || error("翼座標データ $k 行目が固定幅3F9.5形式ではありません")
        xt[k] = parse(Float64, strip(line[1:9]))
        yt[k] = parse(Float64, strip(line[10:18]))
        yb[k] = parse(Float64, strip(line[19:27]))
    end
    return wingname, xt, yt, yb
end

function CALBETA!(d::CFDData, airfoil_file::AbstractString)
    wingname, wx_t, wy_t, wy_b = _read_airfoil(airfoil_file)
    theta = PI*d.DEG/180.0
    c = cos(-theta); s = sin(-theta)
    wx_b = zeros(Float64, IGS)
    for ii in 1:IGS
        x0 = wx_t[ii]
        twxt = x0*c - wy_t[ii]*s
        twyt = x0*s + wy_t[ii]*c
        twxb = x0*c - wy_b[ii]*s
        twyb = x0*s + wy_b[ii]*c
        wx_t[ii] = twxt + d.XL*0.5 - cos(theta)*0.5
        wy_t[ii] = twyt + d.YL*0.5 + sin(theta)*0.5 - 0.06
        # Fortran uses the original upper-surface x-coordinate in both transforms.
        # This is also the source expression for the lower-surface x transform.
        # wx_b is kept separately to avoid overwriting upper-surface coordinates.
        wx_b[ii] = twxb + d.XL*0.5 - cos(theta)*0.5
        wy_b[ii] = twyb + d.YL*0.5 + sin(theta)*0.5 - 0.06
    end

    xtmp = z1(-3, LX+5); ytmp = z1(-3, LY+5)
    xalp1 = 1.0; xalp2 = 3.0
    for i in 1:IG
        x1 = 2.0*Float64(i)/Float64(IG)-1.0
        xtmp[i] = d.XL*(((xalp1*(x1^xalp2)+x1)*0.5)/(xalp1+1.0)+0.5)
    end
    yalp1 = 1.0; yalp2 = 5.0
    for j in 1:JG
        y1 = 2.0*Float64(j)/Float64(JG)-1.0
        ytmp[j] = d.YL*(((yalp1*(y1^yalp2)+y1)*0.5)/(yalp1+1.0)+0.5)
    end
    xtmp[0] = 0.0; xtmp[-1] = -xtmp[1]; xtmp[-2] = -xtmp[2]; xtmp[-3] = -xtmp[3]
    ytmp[0] = 0.0; ytmp[-1] = -ytmp[1]; ytmp[-2] = -ytmp[2]; ytmp[-3] = -ytmp[3]
    for k in 1:5
        xtmp[IG+k] = 2.0*d.XL - xtmp[IG-k]
        ytmp[JG+k] = 2.0*d.YL - ytmp[JG-k]
    end
    for i in -2:LX+5
        d.X[i] = xtmp[i]
        d.DX[i] = d.X[i] - d.X[i-1]
        d.XX[i] = (d.X[i] + d.X[i-1])*0.5
        d.DXX[i] = d.XX[i] - d.XX[i-1]
    end
    for j in -2:LY+5
        d.Y[j] = ytmp[j]
        d.DY[j] = d.Y[j] - d.Y[j-1]
        d.YY[j] = (d.Y[j] + d.Y[j-1])*0.5
        d.DYY[j] = d.YY[j] - d.YY[j-1]
    end

    # Extend BETA ghost cells as free-stream fluid. The original COMMON
    # declaration ends at IG,JG, although IDY and FS0 read IG+1/JG+1.
    for i in 0:IG+1, j in 0:JG+1
        d.BETA[i,j,1] = 1.0; d.BETA[i,j,2] = 1.0
    end
    # The source loops to IGS while reading II+1 from 52-point arrays.
    # Restrict to IGS-1 to avoid the original out-of-bounds final iteration.
    for i in 0:IG
        ib = 0; it = 0
        for ii in 1:IGS-1
            if wx_b[ii+1] > d.X[i] && d.X[i] >= wx_b[ii]; ib = ii; end
            if wx_t[ii+1] > d.X[i] && d.X[i] >= wx_t[ii]; it = ii; end
        end
        if ib*it != 0
            xt_a = d.X[i]-wx_t[it]; xb_a = d.X[i]-wx_b[ib]
            xt_b = wx_t[it+1]-d.X[i]; xb_b = wx_b[ib+1]-d.X[i]
            yt = (wy_t[it]*xt_b + wy_t[it+1]*xt_a)/(xt_a+xt_b)
            yb = (wy_b[ib]*xb_b + wy_b[ib+1]*xb_a)/(xb_a+xb_b)
            for j in 0:JG
                if d.YY[j] >= yb && d.YY[j] < yt
                    d.BETA[i,j,1] = 0.0
                end
                if d.YY[j] <= yb && d.YY[j+1] > yb
                    d.GAMMA[i,j,1] = (yb-d.YY[j])/(yb-d.YY[j-1])
                end
                if d.YY[j] >= yt && d.YY[j-1] < yt
                    d.GAMMA[i,j,1] = (d.YY[j]-yt)/(d.YY[j+1]-yt)
                end
            end
        end
    end
    for i in 0:IG
        ib = 0; it = 0
        for ii in 1:IGS-1
            if wx_b[ii+1] > d.XX[i] && d.XX[i] >= wx_b[ii]; ib = ii; end
            if wx_t[ii+1] > d.XX[i] && d.XX[i] >= wx_t[ii]; it = ii; end
        end
        if ib*it != 0
            xt_a = d.XX[i]-wx_t[it]; xb_a = d.XX[i]-wx_b[ib]
            xt_b = wx_t[it+1]-d.XX[i]; xb_b = wx_b[ib+1]-d.XX[i]
            yt = (wy_t[it]*xt_b + wy_t[it+1]*xt_a)/(xt_a+xt_b)
            yb = (wy_b[ib]*xb_b + wy_b[ib+1]*xb_a)/(xb_a+xb_b)
            for j in 0:JG
                if d.Y[j] >= yb && d.Y[j] < yt
                    d.BETA[i,j,2] = 0.0
                end
            end
        end
    end
    println("GRID GENERATION OK!")

    for i in 1:IG
        for j in 1:JG-1
            d.IDX[i,j] = 4*Int(d.BETA[i,j-1,1]) + 2*Int(d.BETA[i,j,1]) + Int(d.BETA[i,j+1,1])
            d.IDY[i,j] = 4*Int(d.BETA[i-1,j,2]) + 2*Int(d.BETA[i,j,2]) + Int(d.BETA[i+1,j,2])
        end
    end
    for i in 0:IG+1
        d.DX2[i] = d.DX[i]+d.DX[i+1]
        d.DXX2[i] = d.DXX[i]+d.DXX[i+1]
        d.PDX[i] = 1.0/d.DX[i]
        d.PDXX[i] = 1.0/d.DXX[i]
        d.PDX2[i] = 1.0/(d.DX[i]+d.DX[i+1])
        d.PDX3[i] = 1.0/(d.DX[i]*d.DX[i+1]*(d.DX[i]+d.DX[i+1]))
        d.PDXX3[i] = 1.0/(d.DXX[i]*d.DXX[i+1]*(d.DXX[i]+d.DXX[i+1]))
    end
    for j in 0:JG+1
        d.DY2[j] = d.DY[j]+d.DY[j+1]
        d.DYY2[j] = d.DYY[j]+d.DYY[j+1]
        d.PDY[j] = 1.0/d.DY[j]
        d.PDYY[j] = 1.0/d.DYY[j]
        d.PDY2[j] = 1.0/(d.DY[j]+d.DY[j+1])
        d.PDY3[j] = 1.0/(d.DY[j]*d.DY[j+1]*(d.DY[j]+d.DY[j+1]))
        d.PDYY3[j] = 1.0/(d.DYY[j]*d.DYY[j+1]*(d.DYY[j]+d.DYY[j+1]))
    end
    for i in 1:IG
        for j in 1:JG
            d.D_W[i] = d.PDX[i]*d.PDXX[i]
            d.D_E[i] = d.PDX[i]*d.PDXX[i+1]
            d.D_S[j] = d.PDY[j]*d.PDYY[j]
            d.D_N[j] = d.PDY[j]*d.PDYY[j+1]
            d.PDPHI[i,j] = 1.0/(d.D_W[i]+d.D_E[i]+d.D_S[j]+d.D_N[j])
        end
    end
    dfn = 0.0
    for i in 1:IG
        for j in 1:JG
            dfnn = d.DT/d.RE*(1.0/(d.DX[i]^2)+1.0/(d.DY[j]^2))
            dfn = max(dfn, dfnn)
        end
    end
    @printf("DIFFUSION NUMBER = %11.4E\n", dfn)
    return wingname
end
