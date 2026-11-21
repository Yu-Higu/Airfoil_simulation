# ======================================================================
# Module AirfoilSimulation
# Two-dimensional airfoil flow simulation (Standard Array + LoopVectorization)
# ======================================================================

module AirfoilSimulation

using Printf
using Polyester
using LoopVectorization

export read_setting, read_wing, allocate_sim, calbeta!, fs_predictor!, fs0!, fs1!,
    sor_row, fs2!, write_puv, ubc!, output, seedout, seedin!, main, run_simulation,
    TS_HEADER, ts_line

const PI_F = 3.141592653589793

# 1始まり標準配列の確保ヘルパー
std_oa(::Type{T}, dims...) where {T} = zeros(T, dims...)

# Fortran Ew.d 編集記述子
const EFMT = Dict(d => Printf.Format("%." * string(d - 1) * "e") for d in 1:12)

function fE(x::Real, w::Int, d::Int)
    xf = Float64(x)
    if !isfinite(xf)
        return lpad(isnan(xf) ? "NaN" : (xf < 0 ? "-Infinity" : "Infinity"), w)
    end
    if xf == 0.0
        body = string("0.", "0"^d, "E+00")
    else
        t = Printf.format(EFMT[d], abs(xf))
        k = findfirst('e', t)
        mant = SubString(t, 1, k - 1)
        e10 = parse(Int, SubString(t, k + 1)) + 1
        digs = d > 1 ? string(mant[1], SubString(mant, 3)) : String(mant)
        ae = abs(e10)
        ex = ae < 100 ? string("E", e10 < 0 ? '-' : '+', lpad(ae, 2, '0')) :
             string(e10 < 0 ? '-' : '+', ae)
        body = string("0.", digs, ex)
    end
    s = xf < 0 ? string('-', body) : body
    return length(s) > w ? "*"^w : lpad(s, w)
end

a71(s::AbstractString) = first(rpad(s, 71), 71)

parse_f(s) = parse(Float64, replace(strip(s), r"[dD]" => "e"))
first_token(line) = split(line, r"[\s,]+"; keepempty=false)[1]

const TS_HEADER = "    NN     NUM    TIME    DIV         CFL         AVE.U        C_L      C_D     NCOUNT"

ts_line(num, time, div, cfl, usum, cl, cd, ncount) =
    @sprintf("%6d%8d%10.4f", num, num, time) *
        fE(div, 12, 4) * fE(cfl, 12, 4) * fE(usum, 12, 4) *
        @sprintf("%9.4f%9.4f%8d", cl, cd, ncount)

# ----------------------------------------------------------------------
# setting.dat 読み込み
# ----------------------------------------------------------------------
function read_setting(file::AbstractString)
    lines = readlines(file)
    ig = parse(Int, first_token(lines[1]))
    jg = parse(Int, first_token(lines[2]))
    igs = parse(Int, first_token(lines[3]))
    re = parse_f(first_token(lines[4]))
    deg = parse_f(first_token(lines[5]))
    dt = parse_f(first_token(lines[6]))
    nc = parse(Int, first_token(lines[7]))
    ng = parse(Int, first_token(lines[8]))
    xl = parse_f(first_token(lines[9]))
    yl = parse_f(first_token(lines[10]))
    nwing = String(strip(first_token(lines[11]), ['\'', '"']))
    return (;
        IG=ig, JG=jg, IGS=igs, RE=re, DEG=deg, DT=dt, NC=nc, NG=ng,
        XL=xl, YL=yl, NWING=nwing,
    )
end

# ----------------------------------------------------------------------
# 配列確保（標準Array：1始まり，オフセットシフト考慮）
# ----------------------------------------------------------------------
function allocate_sim(p)
    ig, jg = p.IG, p.JG
    lx, ly = ig + 1, jg + 1

    beta = std_oa(Float64, ig + 2, jg + 2, 2)
    fill!(beta, 1.0)
    gamma = std_oa(Float64, ig + 2, jg + 2, 2)

    return (;
        IG=ig, JG=jg, LX=lx, LY=ly,
        RE=p.RE, DEG=p.DEG, DT=p.DT, NC=p.NC, NG=p.NG,
        XL=p.XL, YL=p.YL, NWING=p.NWING,
        U=std_oa(Float64, ig + 2, jg + 2, 2),
        UT=std_oa(Float64, ig + 2, jg + 2, 2),
        CNX_OLD=std_oa(Float64, ig + 1, jg + 1),
        CNY_OLD=std_oa(Float64, ig + 1, jg + 1),
        P=std_oa(Float64, ig + 9, jg + 5),
        PHI=std_oa(Float64, ig + 9, jg + 5),
        BETA=beta, GAMMA=gamma,
        IDX=std_oa(Int, ig + 1, jg + 1),
        IDY=std_oa(Int, ig + 1, jg + 1),
        X=std_oa(Float64, lx + 9), XX=std_oa(Float64, lx + 9),
        DX=std_oa(Float64, lx + 9), DXX=std_oa(Float64, lx + 9),
        Y=std_oa(Float64, ly + 9), YY=std_oa(Float64, ly + 9),
        DY=std_oa(Float64, ly + 9), DYY=std_oa(Float64, ly + 9),
        PDX=std_oa(Float64, lx + 3), PDXX=std_oa(Float64, lx + 3),
        DX2=std_oa(Float64, lx + 3), DXX2=std_oa(Float64, lx + 3),
        PDX2=std_oa(Float64, lx + 3), PDX3=std_oa(Float64, lx + 3),
        PDXX3=std_oa(Float64, lx + 3),
        PDY=std_oa(Float64, ly + 3), PDYY=std_oa(Float64, ly + 3),
        DY2=std_oa(Float64, ly + 3), DYY2=std_oa(Float64, ly + 3),
        PDY2=std_oa(Float64, ly + 3), PDY3=std_oa(Float64, ly + 3),
        PDYY3=std_oa(Float64, ly + 3),
        D_N=std_oa(Float64, jg + 1), D_E=std_oa(Float64, ig + 1),
        D_W=std_oa(Float64, ig + 1), D_S=std_oa(Float64, jg + 1),
        PDPHI=std_oa(Float64, ig + 1, jg + 1),
        ROWMAX=zeros(Float64, jg),
        RHS=std_oa(Float64, ig, jg),
    )
end

# ----------------------------------------------------------------------
# 翼ファイル読み込み
# ----------------------------------------------------------------------
field(line, a, b) = a > lastindex(line) ? "" : line[a:min(b, lastindex(line))]

function parse_f9(s0)
    s = replace(strip(s0), r"[dD]" => "e")
    isempty(s) && return 0.0
    (occursin('.', s) || occursin('e', s) || occursin('E', s)) && return parse(Float64, s)
    return parse(Int, s) / 1.0e5
end

function read_wing(file::AbstractString)
    lines = readlines(file)
    data = [rstrip(l, ['\r', '\n']) for l in lines[2:end] if !isempty(strip(l))]
    n = length(data)
    xt = zeros(n)
    yt = zeros(n)
    yb = zeros(n)
    for (k, l) in enumerate(data)
        vals = try
            (parse_f9(field(l, 1, 9)), parse_f9(field(l, 10, 18)), parse_f9(field(l, 19, 27)))
        catch
            v = parse_f.(split(l))
            (v[1], v[2], v[3])
        end
        xt[k], yt[k], yb[k] = vals
    end
    return xt, yt, yb
end

# ----------------------------------------------------------------------
# メッシュ生成 (calbeta!)
# ----------------------------------------------------------------------
function calbeta!(s)
    ig, jg, lx, ly = s.IG, s.JG, s.LX, s.LY
    xl, yl, deg, dt, re = s.XL, s.YL, s.DEG, s.DT, s.RE
    x, xx, dx, dxx = s.X, s.XX, s.DX, s.DXX
    y, yy, dy, dyy = s.Y, s.YY, s.DY, s.DYY
    beta, gamma, idx, idy = s.BETA, s.GAMMA, s.IDX, s.IDY
    pdx, pdxx, dx2, dxx2, pdx2, pdx3, pdxx3 = s.PDX, s.PDXX, s.DX2, s.DXX2, s.PDX2, s.PDX3, s.PDXX3
    pdy, pdyy, dy2, dyy2, pdy2, pdy3, pdyy3 = s.PDY, s.PDYY, s.DY2, s.DYY2, s.PDY2, s.PDY3, s.PDYY3
    d_n, d_e, d_w, d_s, pdphi = s.D_N, s.D_E, s.D_W, s.D_S, s.PDPHI

    nwing_file = s.NWING
    if !isfile(nwing_file)
        println(" ERROR: Cannot open wing file: ", nwing_file)
        return 1
    end

    xt, yt, yb = read_wing(nwing_file)
    igs = length(xt)
    @printf(" Detected wing data points (IGS) = %11d\n", igs)

    wng_xt = zeros(igs + 1)
    wng_xb = zeros(igs + 1)
    wng_yt = zeros(igs + 1)
    wng_yb = zeros(igs + 1)
    wng_xt[1:igs] .= xt
    wng_yt[1:igs] .= yt
    wng_yb[1:igs] .= yb

    x_temp = zeros(lx + 9)
    y_temp = zeros(ly + 9)

    theta = PI_F * deg / 180.0

    for ii in 1:igs
        twng_xt = wng_xt[ii] * cos(-theta) - wng_yt[ii] * sin(-theta)
        twng_yt = wng_xt[ii] * sin(-theta) + wng_yt[ii] * cos(-theta)
        twng_xb = wng_xt[ii] * cos(-theta) - wng_yb[ii] * sin(-theta)
        twng_yb = wng_xt[ii] * sin(-theta) + wng_yb[ii] * cos(-theta)

        wng_xt[ii] = twng_xt + xl * 0.5 - cos(theta) * 0.5
        wng_xb[ii] = twng_xb + xl * 0.5 - cos(theta) * 0.5
        wng_yt[ii] = twng_yt + yl * 0.5 + sin(theta) * 0.5 - 0.06
        wng_yb[ii] = twng_yb + yl * 0.5 + sin(theta) * 0.5 - 0.06
    end

    xalp1 = 1.0
    xalp2 = 3.0
    for i in 1:ig
        x1 = 2.0 * Float64(i) / Float64(ig) - 1.0
        x_temp[i+4] = xl * (((xalp1 * (x1^xalp2) + x1) * 0.5) / (xalp1 + 1.0) + 0.5)
    end
    yalp1 = 1.0
    yalp2 = 5.0
    for j in 1:jg
        y1 = 2.0 * Float64(j) / Float64(jg) - 1.0
        y_temp[j+4] = yl * (((yalp1 * (y1^yalp2) + y1) * 0.5) / (yalp1 + 1.0) + 0.5)
    end

    x_temp[0+4] = 0.0
    x_temp[-1+4] = -x_temp[1+4]
    x_temp[-2+4] = -x_temp[2+4]
    x_temp[-3+4] = -x_temp[3+4]
    x_temp[ig+1+4] = 2.0 * xl - x_temp[ig-1+4]
    x_temp[ig+2+4] = 2.0 * xl - x_temp[ig-2+4]
    x_temp[ig+3+4] = 2.0 * xl - x_temp[ig-3+4]
    x_temp[ig+4+4] = 2.0 * xl - x_temp[ig-4+4]
    x_temp[ig+5+4] = 2.0 * xl - x_temp[ig-5+4]

    y_temp[0+4] = 0.0
    y_temp[-1+4] = -y_temp[1+4]
    y_temp[-2+4] = -y_temp[2+4]
    y_temp[-3+4] = -y_temp[3+4]
    y_temp[jg+1+4] = 2.0 * yl - y_temp[jg-1+4]
    y_temp[jg+2+4] = 2.0 * yl - y_temp[jg-2+4]
    y_temp[jg+3+4] = 2.0 * yl - y_temp[jg-3+4]
    y_temp[jg+4+4] = 2.0 * yl - y_temp[jg-4+4]
    y_temp[jg+5+4] = 2.0 * yl - y_temp[jg-5+4]

    for i in -2:(lx+5)
        x[i+4] = x_temp[i+4]
        dx[i+4] = x[i+4] - x[i-1+4]
        xx[i+4] = (x[i+4] + x[i-1+4]) * 0.5
        dxx[i+4] = xx[i+4] - xx[i-1+4]
    end
    for j in -2:(ly+5)
        y[j+4] = y_temp[j+4]
        dy[j+4] = y[j+4] - y[j-1+4]
        yy[j+4] = (y[j+4] + y[j-1+4]) * 0.5
        dyy[j+4] = yy[j+4] - yy[j-1+4]
    end

    for i in 0:ig, j in 0:jg
        beta[i+1, j+1, 1] = 1.0
        beta[i+1, j+1, 2] = 1.0
    end

    for i in 0:ig
        ib = 0
        it = 0
        xi = x[i+4]
        for ii in 1:igs
            if (xi < wng_xb[ii+1]) && (xi >= wng_xb[ii]); ib = ii; end
            if (xi < wng_xt[ii+1]) && (xi >= wng_xt[ii]); it = ii; end
        end
        (ib * it == 0) && continue
        xta = xi - wng_xt[it]
        xba = xi - wng_xb[ib]
        xtb = wng_xt[it+1] - xi
        xbb = wng_xb[ib+1] - xi
        yt = (wng_yt[it] * xtb + wng_yt[it+1] * xta) / (xta + xtb)
        yb = (wng_yb[ib] * xbb + wng_yb[ib+1] * xba) / (xba + xbb)
        for j in 0:jg
            yyj = yy[j+4]
            if (yyj >= yb) && (yyj < yt)
                beta[i+1, j+1, 1] = 0.0
            end
            if (yyj <= yb) && (yy[j+1+4] > yb)
                gamma[i+1, j+1, 1] = (yb - yyj) / (yb - yy[j-1+4])
            end
            if (yyj >= yt) && (yy[j-1+4] < yt)
                gamma[i+1, j+1, 1] = (yyj - yt) / (yy[j+1+4] - yt)
            end
        end
    end

    for i in 0:ig
        ib = 0
        it = 0
        xxi = xx[i+4]
        for ii in 1:igs
            if (xxi < wng_xb[ii+1]) && (xxi >= wng_xb[ii]); ib = ii; end
            if (xxi < wng_xt[ii+1]) && (xxi >= wng_xt[ii]); it = ii; end
        end
        (ib * it == 0) && continue
        xta = xxi - wng_xt[it]
        xba = xxi - wng_xb[ib]
        xtb = wng_xt[it+1] - xxi
        xbb = wng_xb[ib+1] - xxi
        yt = (wng_yt[it] * xtb + wng_yt[it+1] * xta) / (xta + xtb)
        yb = (wng_yb[ib] * xbb + wng_yb[ib+1] * xba) / (xba + xbb)
        for j in 0:jg
            yj = y[j+4]
            if (yj >= yb) && (yj < yt)
                beta[i+1, j+1, 2] = 0.0
            end
        end
    end

    println(" GRID GENERATION OK!")

    for i in 1:ig, j in 1:(jg-1)
        idx[i+1, j+1] = 4 * Int(beta[i+1, j-1+1, 1]) + 2 * Int(beta[i+1, j+1, 1]) + 1 * Int(beta[i+1, j+1+1, 1])
        idy[i+1, j+1] = 4 * Int(beta[i-1+1, j+1, 2]) + 2 * Int(beta[i+1, j+1, 2]) + 1 * Int(beta[i+1+1, j+1, 2])
    end

    for i in 0:(ig+1)
        dx2[i+1] = dx[i+4] + dx[i+1+4]
        dxx2[i+1] = dxx[i+4] + dxx[i+1+4]
        pdx[i+1] = 1.0 / dx[i+4]
        pdxx[i+1] = 1.0 / dxx[i+4]
        pdx2[i+1] = 1.0 / (dx[i+4] + dx[i+1+4])
        pdx3[i+1] = 1.0 / (dx[i+4] * dx[i+1+4] * (dx[i+4] + dx[i+1+4]))
        pdxx3[i+1] = 1.0 / (dxx[i+4] * dxx[i+1+4] * (dxx[i+4] + dxx[i+1+4]))
    end
    for j in 0:(jg+1)
        dy2[j+1] = dy[j+4] + dy[j+1+4]
        dyy2[j+1] = dyy[j+4] + dyy[j+1+4]
        pdy[j+1] = 1.0 / dy[j+4]
        pdyy[j+1] = 1.0 / dyy[j+4]
        pdy2[j+1] = 1.0 / (dy[j+4] + dy[j+1+4])
        pdy3[j+1] = 1.0 / (dy[j+4] * dy[j+1+4] * (dy[j+4] + dy[j+1+4]))
        pdyy3[j+1] = 1.0 / (dyy[j+4] * dyy[j+1+4] * (dyy[j+4] + dyy[j+1+4]))
    end

    for i in 1:ig
        d_w[i+1] = pdx[i+1] * pdxx[i+1]
        d_e[i+1] = pdx[i+1] * pdxx[i+1+1]
    end
    for j in 1:jg
        d_s[j+1] = pdy[j+1] * pdyy[j+1]
        d_n[j+1] = pdy[j+1] * pdyy[j+1+1]
    end
    for i in 1:ig, j in 1:jg
        pdphi[i+1, j+1] = 1.0 / (d_w[i+1] + d_e[i+1] + d_s[j+1] + d_n[j+1])
    end

    dfn = 0.0
    for i in 1:ig, j in 1:jg
        dfnn = dt / re * (1.0 / (dx[i+4]^2) + 1.0 / (dy[j+4]^2))
        dfn = max(dfn, dfnn)
    end
    println(lpad("DIFFUSION NUMBER = ", 20), fE(dfn, 11, 4))

    return 0
end

# ----------------------------------------------------------------------
# 予測子ステップ (fs_predictor!) - @batch と @turbo の併用
# ----------------------------------------------------------------------
function fs_predictor!(s, first_step::Bool)
    ig, jg, re, dt = s.IG, s.JG, s.RE, s.DT
    u, ut, p, beta, gamma, idx = s.U, s.UT, s.P, s.BETA, s.GAMMA, s.IDX
    cnx_old, cny_old = s.CNX_OLD, s.CNY_OLD
    pdx, pdxx, dx, dx2, pdx2, pdx3, pdxx3, dxx, dxx2 = s.PDX, s.PDXX, s.DX, s.DX2, s.PDX2, s.PDX3, s.PDXX3, s.DXX, s.DXX2
    pdy, pdyy, dy, dy2, pdy2, pdy3, dyy, dyy2, pdyy3 = s.PDY, s.PDYY, s.DY, s.DY2, s.PDY2, s.PDY3, s.DYY, s.DYY2, s.PDYY3

    pre = 1.0 / re

    # u-運動量方程式
    @batch for j in 1:jg
        @turbo for i in 1:ig
            ududx = 0.5 * ((u[i+2, j+1, 1] + u[i+1, j+1, 1]) * (u[i+2, j+1, 1] - u[i+1, j+1, 1]) * pdx[i+2] * dx[i+4] +
                    (u[i+1, j+1, 1] + u[i, j+1, 1]) * (u[i+1, j+1, 1] - u[i, j+1, 1]) * pdx[i+1] * dx[i+5]) * pdx2[i+1]
            vdudy = 0.5 * ((u[i+2, j+1, 2] * dx[i+4] + u[i+1, j+1, 2] * dx[i+5]) * (u[i+1, j+2, 1] - u[i+1, j+1, 1]) * pdyy[j+2] +
                    (u[i+2, j, 2] * dx[i+4] + u[i+1, j, 2] * dx[i+5]) * (u[i+1, j+1, 1] - u[i+1, j, 1]) * pdyy[j+1]) * pdx2[i+1]
            u_con = ududx + vdudy

            uvisx = (dx[i+5] * u[i, j+1, 1] - dx2[i+1] * u[i+1, j+1, 1] + dx[i+4] * u[i+2, j+1, 1]) * 2.0 * pdx3[i+1] * pre
            uvisy = (dyy[j+5] * u[i+1, j, 1] - dyy2[j+1] * u[i+1, j+1, 1] + dyy[j+4] * u[i+1, j+2, 1]) * 2.0 * pdyy3[j+1] * pre
            u_vis = uvisx + uvisy

            dpdx = (p[i+6, j+3] - p[i+5, j+3]) * pdxx[i+2]

            val = u_con - u_vis + dpdx
            active = beta[i+1, j+1, 1] == 1.0
            rhs_val = first_step ? (u[i+1, j+1, 1] - dt * val) :
                      ((u[i+1, j+1, 1] - dt * (val * 1.5 - cnx_old[i+1, j+1] * 0.5)) * beta[i+1, j+1, 1])
            ut[i+1, j+1, 1] = ifelse(active, rhs_val, ut[i+1, j+1, 1])
            cnx_old[i+1, j+1] = ifelse(active, val, cnx_old[i+1, j+1])
        end
    end

    # v-運動量方程式
    @batch for j in 1:(jg-1)
        @turbo for i in 1:ig
            udvdx = 0.5 * ((u[i+1, j+2, 1] * dy[j+4] + u[i+1, j+1, 1] * dy[j+5]) * (u[i+2, j+1, 2] - u[i+1, j+1, 2]) * pdxx[i+2] +
                    (u[i, j+2, 1] * dy[j+4] + u[i, j+1, 1] * dy[j+5]) * (u[i+1, j+1, 2] - u[i, j+1, 2]) * pdxx[i+1]) * pdy2[j+1]
            vdvdy = 0.5 * ((u[i+1, j+2, 2] + u[i+1, j+1, 2]) * (u[i+1, j+2, 2] - u[i+1, j+1, 2]) * pdy[j+2] * dy[j+4] +
                    (u[i+1, j+1, 2] + u[i+1, j, 2]) * (u[i+1, j+1, 2] - u[i+1, j, 2]) * pdy[j+1] * dy[j+5]) * pdy2[j+1]
            v_con = udvdx + vdvdy

            vvisx = (dxx[i+5] * u[i, j+1, 2] - dxx2[i+1] * u[i+1, j+1, 2] + dxx[i+4] * u[i+2, j+1, 2]) * 2.0 * pdxx3[i+1] * pre
            vvisy = (dy[j+5] * u[i+1, j, 2] - dy2[j+1] * u[i+1, j+1, 2] + dy[j+4] * u[i+1, j+2, 2]) * 2.0 * pdy3[j+1] * pre
            v_vis = vvisx + vvisy

            dpdy = (p[i+5, j+4] - p[i+5, j+3]) * pdyy[j+2]

            val = v_con - v_vis + dpdy
            active = beta[i+1, j+1, 2] == 1.0
            rhs_val = first_step ? ((u[i+1, j+1, 2] - dt * val) * beta[i+1, j+1, 2]) :
                      ((u[i+1, j+1, 2] - dt * (val * 1.5 - cny_old[i+1, j+1] * 0.5)) * beta[i+1, j+1, 2])
            ut[i+1, j+1, 2] = ifelse(active, rhs_val, ut[i+1, j+1, 2])
            cny_old[i+1, j+1] = ifelse(active, val, cny_old[i+1, j+1])
        end
    end

    # 翼面近傍カットセル補正
    @batch for i in 1:ig
        @inbounds for j in 1:jg
            if idx[i+1, j+1] == 3
                ut[i+1, j+1, 1] = ut[i+1, j+2, 1] * gamma[i+1, j+1, 1]
            end
            if idx[i+1, j+1] == 5
                ut[i+1, j+1, 1] = ut[i+1, j, 1] * gamma[i+1, j+1, 1]
            end
        end
    end
    return nothing
end

function fs0!(s)
    ig, jg = s.IG, s.JG
    u, p, beta = s.U, s.P, s.BETA
    @batch for j in 0:(jg+1)
        @turbo for i in 0:(ig+1)
            u[i+1, j+1, 1] = 1.0 * beta[i+1, j+1, 1]
            u[i+1, j+1, 2] = 0.0 * beta[i+1, j+1, 2]
        end
    end
    fill!(p, 0.0)
    fs_predictor!(s, true)
end

fs1!(s) = fs_predictor!(s, false)

# ----------------------------------------------------------------------
# 圧力ポアソン方程式 (sor_row / fs2!)
# ----------------------------------------------------------------------
# sor_row はステップ幅 2 の更新と最大値探索を行うため @inbounds @fastmath を使用
function sor_row(phi, d_w, d_e, d_s, d_n, rhs, pdphi, omg, ig, j, color)
    j_start = 1 + ((j + color) & 1)
    m = 0.0
    pj = j + 3
    @inbounds @fastmath for i in j_start:2:ig
        pi = i + 5
        phi_calc = (d_w[i+1] * phi[pi-1, pj] + d_e[i+1] * phi[pi+1, pj] +
            d_s[j+1] * phi[pi, pj-1] + d_n[j+1] * phi[pi, pj+1] - rhs[i, j]) * pdphi[i+1, j+1]
        phidif = phi_calc - phi[pi, pj]
        phi[pi, pj] = phi[pi, pj] + omg * phidif
        a = abs(phidif)
        if a > m
            m = a
        end
    end
    return m
end

function fs2!(s)
    ig, jg, dt, xl = s.IG, s.JG, s.DT, s.XL
    u, ut, p, phi, x = s.U, s.UT, s.P, s.PHI, s.X
    pdx, pdy, pdxx, pdyy, pdphi = s.PDX, s.PDY, s.PDXX, s.PDYY, s.PDPHI
    d_w, d_e, d_s, d_n = s.D_W, s.D_E, s.D_S, s.D_N
    rowmax, rhs = s.ROWMAX, s.RHS

    igf = 12

    # 流出減衰領域
    @batch for j in 1:jg
        @inbounds for i in (div(ig*(igf-1), igf)+1):ig
            damp = -tanh(4.5 * (x[i+4] - xl))
            ut[i+1, j+1, 1] = (ut[i+1, j+1, 1] - 1.0) * damp + 1.0
            ut[i+1, j+1, 2] = ut[i+1, j+1, 2] * damp
        end
    end

    for i in 1:ig
        ut[i+1, 1, 1] = ut[i+1, 2, 1]
        ut[i+1, jg+2, 1] = ut[i+1, jg+1, 1]
        ut[i+1, 1, 2] = 0.0
        ut[i+1, jg+1, 2] = 0.0
    end
    for j in 1:jg
        ut[1, j+1, 1] = 1.0
        ut[ig+2, j+1, 1] = ut[ig+1, j+1, 1]
        ut[1, j+1, 2] = -ut[2, j+1, 2]
        ut[ig+2, j+1, 2] = -ut[ig+1, j+1, 2]
    end

    pdt = 1.0 / dt
    omg = 1.4
    eps = Float64(1.0f-4)
    ncount_final = 50001

    @batch for j in 1:jg
        @turbo for i in 1:ig
            rhs[i, j] = pdt * ((ut[i+1, j+1, 1] - ut[i, j+1, 1]) * pdx[i+1] +
                (ut[i+1, j+1, 2] - ut[i+1, j, 2]) * pdy[j+1])
        end
    end

    for nc_it in 1:50000
        ermax = 0.0
        for color in 0:1
            @batch for j in 1:jg
                @inbounds rowmax[j] = sor_row(phi, d_w, d_e, d_s, d_n, rhs, pdphi, omg, ig, j, color)
            end
            ermax_color = maximum(rowmax)
            if ermax_color > ermax
                ermax = ermax_color
            end

            for i in 1:ig
                phi[i+5, 3] = phi[i+5, 4]
                phi[i+5, jg+4] = phi[i+5, jg+3]
            end
            for j in 1:jg
                phi[5, j+3] = phi[6, j+3]
                phi[ig+6, j+3] = phi[ig+5, j+3]
            end
        end

        if ermax < eps
            ncount_final = nc_it
            break
        end
    end

    ncount = ncount_final

    p_stn = p[6, 4] + phi[6, 4]
    @batch for j in 1:jg
        @turbo for i in 1:ig
            u[i+1, j+1, 1] = ut[i+1, j+1, 1] - dt * (phi[i+6, j+3] - phi[i+5, j+3]) * pdxx[i+2]
            u[i+1, j+1, 2] = ut[i+1, j+1, 2] - dt * (phi[i+5, j+4] - phi[i+5, j+3]) * pdyy[j+2]
            p[i+5, j+3] = p[i+5, j+3] + phi[i+5, j+3] - p_stn
        end
    end

    return ncount
end

# ----------------------------------------------------------------------
# puv_field????.dat (8E13.5 per line)
# ----------------------------------------------------------------------
function write_puv(s, num)
    ig, jg = s.IG, s.JG
    u, p, beta, xx, yy, dx, dy = s.U, s.P, s.BETA, s.XX, s.YY, s.DX, s.DY
    nstep = @sprintf("%04d", num % 10000)
    open("TEXTDATA/puv_field" * nstep * ".dat", "w") do f20
        for i in 1:ig, j in 1:jg
            vort = -((u[i+1, j+1, 1] + u[i, j+1, 1]) * 0.5 - (u[i+1, j, 1] + u[i, j, 1]) * 0.5) / dy[j+4] +
                ((u[i+1, j+1, 2] + u[i+1, j, 2]) * 0.5 - (u[i, j+1, 2] + u[i, j, 2]) * 0.5) / dx[i+4]
            println(f20,
                fE(xx[i+4], 13, 5), fE(yy[j+4], 13, 5), fE(p[i+5, j+3], 13, 5),
                fE((u[i+1, j+1, 1] + u[i, j+1, 1]) * 0.5, 13, 5),
                fE((u[i+1, j+1, 2] + u[i+1, j, 2]) * 0.5, 13, 5),
                fE((beta[i+1, j+1, 2] + beta[i+1, j, 2]) * 0.5, 13, 5),
                fE((beta[i+1, j+1, 1] + beta[i, j+1, 1]) * 0.5, 13, 5),
                fE(vort, 13, 5))
        end
    end
end

# ----------------------------------------------------------------------
# UBC : 境界条件，DIV / CFL / 平均 U，C_L，C_D
# ----------------------------------------------------------------------
function ubc!(s, num)
    ig, jg, dt, yl = s.IG, s.JG, s.DT, s.YL
    u, p, beta, pdx, pdy, dx, dy = s.U, s.P, s.BETA, s.PDX, s.PDY, s.DX, s.DY

    for i in 1:ig
        u[i+1, 1, 1] = u[i+1, 2, 1]
        u[i+1, jg+2, 1] = u[i+1, jg+1, 1]
        u[i+1, 1, 2] = 0.0
        u[i+1, jg+1, 2] = 0.0
    end
    for j in 1:jg
        u[1, j+1, 1] = 1.0
        u[ig+2, j+1, 1] = u[ig+1, j+1, 1]
        u[1, j+1, 2] = -u[2, j+1, 2]
        u[ig+2, j+1, 2] = -u[ig+1, j+1, 2]
    end
    for i in -3:(ig+3)
        p[i+5, 3] = p[i+5, 4]
        p[i+5, jg+4] = p[i+5, jg+3]
    end
    for j in 1:jg
        p[5, j+3] = p[6, j+3]
        p[ig+6, j+3] = p[ig+5, j+3]
    end

    div = 0.0
    for j in 1:jg, i in 1:ig
        qqqq = (u[i+1, j+1, 1] - u[i, j+1, 1]) * pdx[i+1] + (u[i+1, j+1, 2] - u[i+1, j, 2]) * pdy[j+1]
        qqq = abs(qqqq)
        if qqq > div
            div = qqq
        end
    end

    cfl = 0.0
    for j in 1:jg, i in 1:ig
        cfll = (abs(u[i+1, j+1, 1] + u[i, j+1, 1]) * pdx[i+1] +
            abs(u[i+1, j+1, 2] + u[i+1, j, 2]) * pdy[j+1]) * 0.5 * dt
        if cfll > cfl
            cfl = cfll
        end
    end

    usum = 0.0
    for j in 1:jg, i in 1:ig
        usum = usum + u[i+1, j+1, 1] * dy[j+4]
    end
    usum = usum / Float64(ig) / yl

    if num % 5 == 0
        write_puv(s, num)
    end

    p_sum = 0.0
    for i in 1:ig
        p_top = 0.0
        p_btm = 0.0
        for j in 2:jg
            if beta[i+1, j+1, 2] - beta[i+1, j, 2] >= 0.50
                p_top = p[i+5, j+4]
            end
            if beta[i+1, j+1, 2] - beta[i+1, j, 2] <= -0.50
                p_btm = p[i+5, j+1]
            end
        end
        p_sum = p_sum + dx[i+4] * (p_btm - p_top)
    end
    cl = p_sum * 2.0

    p_sum = 0.0
    for j in 1:jg
        p_top = 0.0
        p_btm = 0.0
        for i in 2:ig
            if beta[i+1, j+1, 1] - beta[i, j+1, 1] >= 0.5
                p_top = p[i+8, j+3]
            end
            if beta[i+1, j+1, 1] - beta[i, j+1, 1] <= -0.5
                p_btm = p[i+1, j+3]
            end
        end
        p_sum = p_sum + dy[j+4] * (p_btm - p_top)
    end
    cd = p_sum * 2.0

    return div, cfl, usum, cl, cd
end

# ----------------------------------------------------------------------
# OUTPUT : TEXTDATA/output.dat
# ----------------------------------------------------------------------
function output(s, num, time)
    ig, jg, ly = s.IG, s.JG, s.LY
    p, beta = s.P, s.BETA
    x, dx, xx, dxx = s.X, s.DX, s.XX, s.DXX
    y, dy, yy, dyy = s.Y, s.DY, s.YY, s.DYY

    if num % 5 == 0
        write_puv(s, num)
    end

    cmtd1 = "#======================================================"
    cmtd2 = "# SIMULATION OF FLOW PAST AN 2D AIRFOIL : " * rpad(s.NWING, 64)
    cmtd3 = "# STATISTICS IN VELOCITY FIELD"
    cprm1 = "# CALC.NUM.     CALC.TIME"
    cprm2 = "# NX  NY  TIMESTEP  REYNOLDS"
    cprm3 = "# X-LENGTH  Y-LENGTH AOA(DEG)"
    cprm4 = "# I X DX XX DXX"
    cprm5 = "# J Y DY YY DYY"
    cprm6 = "# X DX P(UPPER-SURFACE)  P(BOTTOM SURFACE)"
    cprm7 = "# Y DY P(DOWNSTREAM-SURFACE)  P(UPSTREAM SURFACE)"

    cl = 0.0
    cd = 0.0
    open("TEXTDATA/output.dat", "w") do f14
        println(f14, a71(cmtd1))
        println(f14, a71(cmtd2))
        println(f14, a71(cmtd1))
        println(f14, a71(cprm1))
        println(f14, @sprintf("%10d", num), fE(time, 18, 10))
        println(f14, a71(cprm2))
        println(f14, @sprintf("%4d%4d%11.6f%9.1f", ig, jg, s.DT, s.RE))
        println(f14, a71(cprm3))
        println(f14, @sprintf("%10.6f%10.6f%7.2f", s.XL, s.YL, s.DEG))
        println(f14, a71(cmtd3))

        println(f14, a71(cprm4))
        for i in 0:(ig+1)
            println(f14, @sprintf("%3d", i), fE(x[i+4], 14, 6), fE(dx[i+4], 14, 6),
                fE(xx[i+4], 14, 6), fE(dxx[i+4], 14, 6))
        end

        println(f14, a71(cprm5))
        for j in 0:(ly+1)
            println(f14, @sprintf("%3d", j), fE(y[j+4], 14, 6), fE(dy[j+4], 14, 6),
                fE(yy[j+4], 14, 6), fE(dyy[j+4], 14, 6))
        end

        println(f14, a71(cprm6))
        p_sum = 0.0
        for i in 1:ig
            p_top = 0.0
            p_btm = 0.0
            for j in 2:jg
                if beta[i+1, j+1, 2] - beta[i+1, j, 2] >= 0.50
                    p_top = p[i+5, j+4]
                end
                if beta[i+1, j+1, 2] - beta[i+1, j, 2] <= -0.50
                    p_btm = p[i+5, j+1]
                end
            end
            if p_top != 0.0
                println(f14, fE(x[i+4], 14, 6), fE(dx[i+4], 14, 6), fE(p_top, 14, 6), fE(p_btm, 14, 6))
            end
            p_sum = p_sum + dx[i+4] * (p_btm - p_top)
        end
        cl = p_sum * 2.0

        println(f14, a71(cprm7))
        p_sum = 0.0
        for j in 1:jg
            p_top = 0.0
            p_btm = 0.0
            for i in 2:ig
                if beta[i+1, j+1, 1] - beta[i, j+1, 1] >= 0.5
                    p_top = p[i+8, j+3]
                end
                if beta[i+1, j+1, 1] - beta[i, j+1, 1] <= -0.5
                    p_btm = p[i+1, j+3]
                end
            end
            if p_top != 0.0
                println(f14, fE(y[j+4], 14, 6), fE(dy[j+4], 14, 6), fE(p_top, 14, 6), fE(p_btm, 14, 6))
            end
            p_sum = p_sum + dy[j+4] * (p_btm - p_top)
        end
        cd = p_sum * 2.0
    end

    @printf(" C_L =%6.3f , C_D =%6.3f\n", cl, cd)
end

# ----------------------------------------------------------------------
# SEED DATA : リスタートファイル
# ----------------------------------------------------------------------
function seedout(s, num, time)
    open("seed_u.dat", "w") do f31
        @printf(f31, "%d %.16e\n", num, time)
        @printf(f31, "%d %d %.16e %.16e %.16e %.16e\n", s.IG, s.JG, s.XL, s.YL, s.DT, s.RE)
        for v in s.U
            @printf(f31, "%.16e\n", v)
        end
    end
    open("seed_p.dat", "w") do f32
        for v in s.P
            @printf(f32, "%.16e\n", v)
        end
    end
    open("seed_old.dat", "w") do f33
        for v in s.CNX_OLD
            @printf(f33, "%.16e\n", v)
        end
        for v in s.CNY_OLD
            @printf(f33, "%.16e\n", v)
        end
    end
end

function seedin!(s)
    tokens(f) = split(read(f, String), r"[\s,]+"; keepempty=false)
    tu = tokens("seed_u.dat")
    num = parse(Int, tu[1])
    time = parse_f(tu[2])
    igd = parse(Int, tu[3])
    jgd = parse(Int, tu[4])
    xld = parse_f(tu[5])
    yld = parse_f(tu[6])
    dt = parse_f(tu[7])
    re = parse_f(tu[8])

    if igd != s.IG || jgd != s.JG
        println(" MESH DATA ERROR")
        return 1000, num, time, dt, re
    end
    if xld != s.XL || yld != s.YL
        println(" DOMAIN DATA ERROR")
        return 1000, num, time, dt, re
    end

    nu = length(s.U)
    copyto!(vec(s.U), parse_f.(tu[9:(8+nu)]))
    copyto!(vec(s.P), parse_f.(tokens("seed_p.dat")[1:length(s.P)]))
    to = tokens("seed_old.dat")
    nx = length(s.CNX_OLD)
    ny = length(s.CNY_OLD)
    copyto!(vec(s.CNX_OLD), parse_f.(to[1:nx]))
    copyto!(vec(s.CNY_OLD), parse_f.(to[(nx+1):(nx+ny)]))

    println(" READ VELOCITY-SEED DATA OK!")
    return 0, num, time, dt, re
end

# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
function run_simulation(setting_file::AbstractString="setting.dat")
    p = read_setting(setting_file)
    s = allocate_sim(p)
    dt, nc, ng = s.DT, s.NC, s.NG

    num = 0
    time = 0.0
    avetim = 0.0

    ierr = calbeta!(s)
    ierr != 0 && return

    mkpath("TEXTDATA")
    div, cfl, usum, cl, cd = ubc!(s, num)

    num += 1
    time += dt
    avetim += dt

    num == 1 ? fs0!(s) : fs1!(s)
    ncount = fs2!(s)
    div, cfl, usum, cl, cd = ubc!(s, num)

    open("TEXTDATA/time-series.dat", "w") do f1
        println(f1, TS_HEADER)
        println(TS_HEADER)
        line = ts_line(num, time, div, cfl, usum, cl, cd, ncount)
        println(f1, line)
        println(line)

        for nst in 2:nc
            num += 1
            time += dt
            avetim += dt

            fs1!(s)
            ncount = fs2!(s)
            div, cfl, usum, cl, cd = ubc!(s, num)
            if num % ng == 0
                line = ts_line(num, time, div, cfl, usum, cl, cd, ncount)
                println(f1, line)
                println(line)
            end
        end
    end

    seedout(s, num, time)
    output(s, num, time)
    return nothing
end

main() = run_simulation()

end # module AirfoilSimulation

if abspath(PROGRAM_FILE) == @__FILE__
    AirfoilSimulation.main()
end