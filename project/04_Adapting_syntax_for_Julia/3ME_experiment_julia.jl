# Two-dimensional airfoil flow simulation on a staggered finite-difference grid.
# Adapted to idiomatic Julia naming and local-variable conventions while
# preserving the numerical scheme, output formats, and public state fields.

using Printf
using OffsetArrays

const PI_F = 3.141592653589793          # REAL*8, PARAMETER :: PI

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
# zero-filled array with Fortran-style index ranges, e.g. oa(Float64, 0:IG+1, -2:JG+2)
oa(::Type{t}, rs::AbstractUnitRange...) where {t} =
    OffsetArray(zeros(t, map(length, rs)...), rs...)

# Fortran Ew.d edit descriptor  (0.dddddE+ee)
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

# FORMAT(A71): CHARACTER*71 variable, padded / truncated to 71 characters
a71(s::AbstractString) = first(rpad(s, 71), 71)

# list-directed style number parsing (accepts Fortran 'd' exponents)
parse_f(s) = parse(Float64, replace(strip(s), r"[dD]" => "e"))
first_token(line) = split(line, r"[\s,]+"; keepempty=false)[1]

const TS_HEADER = "    NN     NUM    TIME    DIV         CFL         AVE.U        C_L      C_D     NCOUNT"

# FORMAT(I6,I8,F10.4,3E12.4,2F9.4,I8)
ts_line(num, time, div, cfl, usum, cl, cd, ncount) =
    @sprintf("%6d%8d%10.4f", num, num, time) *
        fE(div, 12, 4) * fE(cfl, 12, 4) * fE(usum, 12, 4) *
        @sprintf("%9.4f%9.4f%8d", cl, cd, ncount)

# ----------------------------------------------------------------------
# setting.dat   (READ(50,*) value, comment  -> first token of each line)
# ----------------------------------------------------------------------
function read_setting(file::AbstractString)
    lines = readlines(file)
    ig = parse(Int, first_token(lines[1]))
    jg = parse(Int, first_token(lines[2]))
    igs = parse(Int, first_token(lines[3]))     # overwritten by auto-detection
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
# allocation (same bounds as the ALLOCATE statements in the Fortran code)
# ----------------------------------------------------------------------
function allocate_sim(p)
    ig, jg = p.IG, p.JG
    lx, ly = ig + 1, jg + 1

    # BETA/GAMMA are declared (0:IG,0:JG) in Fortran but FS0 touches
    # IG+1 / JG+1 as well; those extra cells are 1.0 (free stream) here.
    beta = oa(Float64, 0:(ig+1), 0:(jg+1), 1:2)
    fill!(beta, 1.0)
    gamma = oa(Float64, 0:(ig+1), 0:(jg+1), 1:2)

    return (;
        IG=ig, JG=jg, LX=lx, LY=ly,
        RE=p.RE, DEG=p.DEG, DT=p.DT, NC=p.NC, NG=p.NG,
        XL=p.XL, YL=p.YL, NWING=p.NWING,
        U=oa(Float64, 0:(ig+1), 0:(jg+1), 1:2),
        UT=oa(Float64, 0:(ig+1), 0:(jg+1), 1:2),
        CNX_OLD=oa(Float64, 0:ig, 0:jg),
        CNY_OLD=oa(Float64, 0:ig, 0:jg),
        P=oa(Float64, -4:(ig+4), -2:(jg+2)),
        PHI=oa(Float64, -4:(ig+4), -2:(jg+2)),      # PHI = 0 (warm start later)
        BETA=beta, GAMMA=gamma,
        IDX=oa(Int, 0:ig, 0:jg),
        IDY=oa(Int, 0:ig, 0:jg),
        X=oa(Float64, -3:(lx+5)), XX=oa(Float64, -3:(lx+5)),
        DX=oa(Float64, -3:(lx+5)), DXX=oa(Float64, -3:(lx+5)),
        Y=oa(Float64, -3:(ly+5)), YY=oa(Float64, -3:(ly+5)),
        DY=oa(Float64, -3:(ly+5)), DYY=oa(Float64, -3:(ly+5)),
        PDX=oa(Float64, 0:(lx+2)), PDXX=oa(Float64, 0:(lx+2)),
        DX2=oa(Float64, 0:(lx+2)), DXX2=oa(Float64, 0:(lx+2)),
        PDX2=oa(Float64, 0:(lx+2)), PDX3=oa(Float64, 0:(lx+2)),
        PDXX3=oa(Float64, 0:(lx+2)),
        PDY=oa(Float64, 0:(ly+2)), PDYY=oa(Float64, 0:(ly+2)),
        DY2=oa(Float64, 0:(ly+2)), DYY2=oa(Float64, 0:(ly+2)),
        PDY2=oa(Float64, 0:(ly+2)), PDY3=oa(Float64, 0:(ly+2)),
        PDYY3=oa(Float64, 0:(ly+2)),
        D_N=oa(Float64, 0:jg), D_E=oa(Float64, 0:ig),
        D_W=oa(Float64, 0:ig), D_S=oa(Float64, 0:jg),
        PDPHI=oa(Float64, 0:ig, 0:jg),
        ROWMAX=zeros(Float64, jg),
        RHS=oa(Float64, 1:ig, 1:jg),
    )
end

# ----------------------------------------------------------------------
# wing file : header line + (XT, YT, YB) in FORMAT(3F9.5), blank lines ignored
# ----------------------------------------------------------------------
field(line, a, b) = a > lastindex(line) ? "" : line[a:min(b, lastindex(line))]

function parse_f9(s0)
    s = replace(strip(s0), r"[dD]" => "e")
    isempty(s) && return 0.0
    (occursin('.', s) || occursin('e', s) || occursin('E', s)) && return parse(Float64, s)
    return parse(Int, s) / 1.0e5                 # implied decimal point (F9.5)
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
            v = parse_f.(split(l))               # fallback: whitespace separated
            (v[1], v[2], v[3])
        end
        xt[k], yt[k], yb[k] = vals
    end
    return xt, yt, yb
end

# ----------------------------------------------------------------------
# MESH GENERATION AND BETA CALCULATION
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
    igs = length(xt)                              # auto-detected
    @printf(" Detected wing data points (IGS) = %11d\n", igs)

    wng_xt = zeros(igs + 1)
    wng_xb = zeros(igs + 1)
    wng_yt = zeros(igs + 1)
    wng_yb = zeros(igs + 1)
    wng_xt[1:igs] .= xt
    wng_yt[1:igs] .= yt
    wng_yb[1:igs] .= yb

    x_temp = oa(Float64, -3:(lx+5))
    y_temp = oa(Float64, -3:(ly+5))

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
        x_temp[i] = xl * (((xalp1 * (x1^xalp2) + x1) * 0.5) / (xalp1 + 1.0) + 0.5)
    end
    yalp1 = 1.0
    yalp2 = 5.0
    for j in 1:jg
        y1 = 2.0 * Float64(j) / Float64(jg) - 1.0
        y_temp[j] = yl * (((yalp1 * (y1^yalp2) + y1) * 0.5) / (yalp1 + 1.0) + 0.5)
    end

    x_temp[0] = 0.0
    x_temp[-1] = -x_temp[1]
    x_temp[-2] = -x_temp[2]
    x_temp[-3] = -x_temp[3]
    x_temp[ig+1] = 2.0 * xl - x_temp[ig-1]
    x_temp[ig+2] = 2.0 * xl - x_temp[ig-2]
    x_temp[ig+3] = 2.0 * xl - x_temp[ig-3]
    x_temp[ig+4] = 2.0 * xl - x_temp[ig-4]
    x_temp[ig+5] = 2.0 * xl - x_temp[ig-5]

    y_temp[0] = 0.0
    y_temp[-1] = -y_temp[1]
    y_temp[-2] = -y_temp[2]
    y_temp[-3] = -y_temp[3]
    y_temp[jg+1] = 2.0 * yl - y_temp[jg-1]
    y_temp[jg+2] = 2.0 * yl - y_temp[jg-2]
    y_temp[jg+3] = 2.0 * yl - y_temp[jg-3]
    y_temp[jg+4] = 2.0 * yl - y_temp[jg-4]
    y_temp[jg+5] = 2.0 * yl - y_temp[jg-5]

    for i in -2:(lx+5)
        x[i] = x_temp[i]
        dx[i] = x[i] - x[i-1]
        xx[i] = (x[i] + x[i-1]) * 0.5
        dxx[i] = xx[i] - xx[i-1]
    end
    for j in -2:(ly+5)
        y[j] = y_temp[j]
        dy[j] = y[j] - y[j-1]
        yy[j] = (y[j] + y[j-1]) * 0.5
        dyy[j] = yy[j] - yy[j-1]
    end

    for i in 0:ig, j in 0:jg
        beta[i, j, 1] = 1.0
        beta[i, j, 2] = 1.0
    end

    # --- BETA(:,:,1)  (u-points) and GAMMA ---
    for i in 0:ig
        ib = 0
        it = 0
        for ii in 1:igs
            if (x[i] < wng_xb[ii+1]) && (x[i] >= wng_xb[ii]); ib = ii; end
            if (x[i] < wng_xt[ii+1]) && (x[i] >= wng_xt[ii]); it = ii; end
        end
        (ib * it == 0) && continue
        xta = x[i] - wng_xt[it]
        xba = x[i] - wng_xb[ib]
        xtb = wng_xt[it+1] - x[i]
        xbb = wng_xb[ib+1] - x[i]
        yt = (wng_yt[it] * xtb + wng_yt[it+1] * xta) / (xta + xtb)
        yb = (wng_yb[ib] * xbb + wng_yb[ib+1] * xba) / (xba + xbb)
        for j in 0:jg
            if (yy[j] >= yb) && (yy[j] < yt)
                beta[i, j, 1] = 0.0
            end
            if (yy[j] <= yb) && (yy[j+1] > yb)
                gamma[i, j, 1] = (yb - yy[j]) / (yb - yy[j-1])
            end
            if (yy[j] >= yt) && (yy[j-1] < yt)
                gamma[i, j, 1] = (yy[j] - yt) / (yy[j+1] - yt)
            end
        end
    end

    # --- BETA(:,:,2)  (v-points) ---
    for i in 0:ig
        ib = 0
        it = 0
        for ii in 1:igs
            if (xx[i] < wng_xb[ii+1]) && (xx[i] >= wng_xb[ii]); ib = ii; end
            if (xx[i] < wng_xt[ii+1]) && (xx[i] >= wng_xt[ii]); it = ii; end
        end
        (ib * it == 0) && continue
        xta = xx[i] - wng_xt[it]
        xba = xx[i] - wng_xb[ib]
        xtb = wng_xt[it+1] - xx[i]
        xbb = wng_xb[ib+1] - xx[i]
        yt = (wng_yt[it] * xtb + wng_yt[it+1] * xta) / (xta + xtb)
        yb = (wng_yb[ib] * xbb + wng_yb[ib+1] * xba) / (xba + xbb)
        for j in 0:jg
            if (y[j] >= yb) && (y[j] < yt)
                beta[i, j, 2] = 0.0
            end
        end
    end

    println(" GRID GENERATION OK!")

    for i in 1:ig, j in 1:(jg-1)
        idx[i, j] = 4 * Int(beta[i, j-1, 1]) + 2 * Int(beta[i, j, 1]) + 1 * Int(beta[i, j+1, 1])
        idy[i, j] = 4 * Int(beta[i-1, j, 2]) + 2 * Int(beta[i, j, 2]) + 1 * Int(beta[i+1, j, 2])
    end

    for i in 0:(ig+1)
        dx2[i] = dx[i] + dx[i+1]
        dxx2[i] = dxx[i] + dxx[i+1]
        pdx[i] = 1.0 / dx[i]
        pdxx[i] = 1.0 / dxx[i]
        pdx2[i] = 1.0 / (dx[i] + dx[i+1])
        pdx3[i] = 1.0 / (dx[i] * dx[i+1] * (dx[i] + dx[i+1]))
        pdxx3[i] = 1.0 / (dxx[i] * dxx[i+1] * (dxx[i] + dxx[i+1]))
    end
    for j in 0:(jg+1)
        dy2[j] = dy[j] + dy[j+1]
        dyy2[j] = dyy[j] + dyy[j+1]
        pdy[j] = 1.0 / dy[j]
        pdyy[j] = 1.0 / dyy[j]
        pdy2[j] = 1.0 / (dy[j] + dy[j+1])
        pdy3[j] = 1.0 / (dy[j] * dy[j+1] * (dy[j] + dy[j+1]))
        pdyy3[j] = 1.0 / (dyy[j] * dyy[j+1] * (dyy[j] + dyy[j+1]))
    end

    for i in 1:ig
        d_w[i] = pdx[i] * pdxx[i]
        d_e[i] = pdx[i] * pdxx[i+1]
    end
    for j in 1:jg
        d_s[j] = pdy[j] * pdyy[j]
        d_n[j] = pdy[j] * pdyy[j+1]
    end
    for i in 1:ig, j in 1:jg
        pdphi[i, j] = 1.0 / (d_w[i] + d_e[i] + d_s[j] + d_n[j])
    end

    dfn = 0.0
    for i in 1:ig, j in 1:jg
        dfnn = dt / re * (1.0 / (dx[i]^2) + 1.0 / (dy[j]^2))
        dfn = max(dfn, dfnn)
    end
    println(lpad("DIFFUSION NUMBER = ", 20), fE(dfn, 11, 4))

    return 0
end

# ----------------------------------------------------------------------
# FS0 / FS1 : temporary velocity
#   FS0 : Euler explicit (first step)       FS1 : Adams-Bashforth (2nd step on)
# The two Fortran routines differ only in the time integration, so they
# share one body here (identical arithmetic).
# ----------------------------------------------------------------------
function fs_predictor!(s, first_step::Bool)
    ig, jg, re, dt = s.IG, s.JG, s.RE, s.DT
    u, ut, p, beta, gamma, idx = s.U, s.UT, s.P, s.BETA, s.GAMMA, s.IDX
    cnx_old, cny_old = s.CNX_OLD, s.CNY_OLD
    pdx, pdxx, dx, dx2, pdx2, pdx3, pdxx3, dxx, dxx2 = s.PDX, s.PDXX, s.DX, s.DX2, s.PDX2, s.PDX3, s.PDXX3, s.DXX, s.DXX2
    pdy, pdyy, dy, dy2, pdy2, pdy3, dyy, dyy2, pdyy3 = s.PDY, s.PDYY, s.DY, s.DY2, s.PDY2, s.PDY3, s.DYY, s.DYY2, s.PDYY3

    pre = 1.0 / re

    # ---- u-momentum ----
    Threads.@threads for j in 1:jg
        for i in 1:ig
            if beta[i, j, 1] == 1.0
                ududx = 0.5 * ((u[i+1, j, 1] + u[i, j, 1]) * (u[i+1, j, 1] - u[i, j, 1]) * pdx[i+1] * dx[i] +
                        (u[i, j, 1] + u[i-1, j, 1]) * (u[i, j, 1] - u[i-1, j, 1]) * pdx[i] * dx[i+1]) * pdx2[i]
                vdudy = 0.5 * ((u[i+1, j, 2] * dx[i] + u[i, j, 2] * dx[i+1]) * (u[i, j+1, 1] - u[i, j, 1]) * pdyy[j+1] +
                        (u[i+1, j-1, 2] * dx[i] + u[i, j-1, 2] * dx[i+1]) * (u[i, j, 1] - u[i, j-1, 1]) * pdyy[j]) * pdx2[i]
                u_con = ududx + vdudy

                uvisx = (dx[i+1] * u[i-1, j, 1] - dx2[i] * u[i, j, 1] + dx[i] * u[i+1, j, 1]) * 2.0 * pdx3[i] * pre
                uvisy = (dyy[j+1] * u[i, j-1, 1] - dyy2[j] * u[i, j, 1] + dyy[j] * u[i, j+1, 1]) * 2.0 * pdyy3[j] * pre
                u_vis = uvisx + uvisy

                dpdx = (p[i+1, j] - p[i, j]) * pdxx[i+1]

                if first_step
                    ut[i, j, 1] = u[i, j, 1] - dt * (u_con - u_vis + dpdx)
                else
                    ut[i, j, 1] = (u[i, j, 1] - dt * ((u_con - u_vis + dpdx) * 1.5 - cnx_old[i, j] * 0.5)) * beta[i, j, 1]
                end
                cnx_old[i, j] = u_con - u_vis + dpdx
            end
        end
    end

    # ---- v-momentum ----
    Threads.@threads for j in 1:(jg-1)
        for i in 1:ig
            if beta[i, j, 2] == 1.0
                udvdx = 0.5 * ((u[i, j+1, 1] * dy[j] + u[i, j, 1] * dy[j+1]) * (u[i+1, j, 2] - u[i, j, 2]) * pdxx[i+1] +
                        (u[i-1, j+1, 1] * dy[j] + u[i-1, j, 1] * dy[j+1]) * (u[i, j, 2] - u[i-1, j, 2]) * pdxx[i]) * pdy2[j]
                vdvdy = 0.5 * ((u[i, j+1, 2] + u[i, j, 2]) * (u[i, j+1, 2] - u[i, j, 2]) * pdy[j+1] * dy[j] +
                        (u[i, j, 2] + u[i, j-1, 2]) * (u[i, j, 2] - u[i, j-1, 2]) * pdy[j] * dy[j+1]) * pdy2[j]
                v_con = udvdx + vdvdy

                vvisx = (dxx[i+1] * u[i-1, j, 2] - dxx2[i] * u[i, j, 2] + dxx[i] * u[i+1, j, 2]) * 2.0 * pdxx3[i] * pre
                vvisy = (dy[j+1] * u[i, j-1, 2] - dy2[j] * u[i, j, 2] + dy[j] * u[i, j+1, 2]) * 2.0 * pdy3[j] * pre
                v_vis = vvisx + vvisy

                dpdy = (p[i, j+1] - p[i, j]) * pdyy[j+1]

                if first_step
                    ut[i, j, 2] = (u[i, j, 2] - dt * (v_con - v_vis + dpdy)) * beta[i, j, 2]
                else
                    ut[i, j, 2] = (u[i, j, 2] - dt * ((v_con - v_vis + dpdy) * 1.5 - cny_old[i, j] * 0.5)) * beta[i, j, 2]
                end
                cny_old[i, j] = v_con - v_vis + dpdy
            end
        end
    end

    # ---- cut-cell correction near the wing surface ----
    Threads.@threads for i in 1:ig
        for j in 1:jg
            if idx[i, j] == 3
                ut[i, j, 1] = ut[i, j+1, 1] * gamma[i, j, 1]
            end
            if idx[i, j] == 5
                ut[i, j, 1] = ut[i, j-1, 1] * gamma[i, j, 1]
            end
        end
    end
    return nothing
end

function fs0!(s)
    ig, jg = s.IG, s.JG
    u, p, beta = s.U, s.P, s.BETA
    Threads.@threads for j in 0:(jg+1)
        for i in 0:(ig+1)
            u[i, j, 1] = 1.0 * beta[i, j, 1]
            u[i, j, 2] = 0.0 * beta[i, j, 2]
            p[i, j] = 0.0
        end
    end
    fs_predictor!(s, true)
end

fs1!(s) = fs_predictor!(s, false)

# ----------------------------------------------------------------------
# FS2 : Poisson equation for pressure correction (red-black SOR)
# ----------------------------------------------------------------------
function sor_row(phi, d_w, d_e, d_s, d_n, rhs, pdphi, omg, ig, j, color)
    j_start = 1 + ((j + color) & 1)
    m = 0.0
    @inbounds for i in j_start:2:ig
        phi_calc = (d_w[i] * phi[i-1, j] + d_e[i] * phi[i+1, j] +
            d_s[j] * phi[i, j-1] + d_n[j] * phi[i, j+1] - rhs[i, j]) * pdphi[i, j]
        phidif = phi_calc - phi[i, j]
        phi[i, j] = phi[i, j] + omg * phidif
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

    # outflow damping zone
    Threads.@threads for j in 1:jg
        for i in (div(ig*(igf-1), igf)+1):ig
            ut[i, j, 1] = (ut[i, j, 1] - 1.0) * (-tanh(4.5 * (x[i] - xl))) + 1.0
            ut[i, j, 2] = ut[i, j, 2] * (-tanh(4.5 * (x[i] - xl)))
        end
    end

    for i in 1:ig
        ut[i, 0, 1] = ut[i, 1, 1]
        ut[i, jg+1, 1] = ut[i, jg, 1]
        ut[i, 0, 2] = 0.0
        ut[i, jg, 2] = 0.0
    end
    for j in 1:jg
        ut[0, j, 1] = 1.0
        ut[ig+1, j, 1] = ut[ig, j, 1]
        ut[0, j, 2] = -ut[1, j, 2]
        ut[ig+1, j, 2] = -ut[ig, j, 2]
    end

    pdt = 1.0 / dt
    omg = 1.4
    eps = Float64(1.0f-4)          # Fortran "1.0E-4" is a single-precision literal

    ncount_final = 50001           # value left by a DO loop that never converges
    multi = Threads.nthreads() > 1

    Threads.@threads for j in 1:jg
        @inbounds @simd for i in 1:ig
            rhs[i, j] = pdt * ((ut[i, j, 1] - ut[i-1, j, 1]) * pdx[i] +
                (ut[i, j, 2] - ut[i, j-1, 2]) * pdy[j])
        end
    end

    for nc_it in 1:50000
        ermax = 0.0
        for color in 0:1
            if multi
                Threads.@threads for j in 1:jg
                    @inbounds rowmax[j] = sor_row(phi, d_w, d_e, d_s, d_n, rhs, pdphi, omg, ig, j, color)
                end
            else
                for j in 1:jg
                    @inbounds rowmax[j] = sor_row(phi, d_w, d_e, d_s, d_n, rhs, pdphi, omg, ig, j, color)
                end
            end
            ermax_color = maximum(rowmax)
            if ermax_color > ermax
                ermax = ermax_color
            end

            for i in 1:ig
                phi[i, 0] = phi[i, 1]
                phi[i, jg+1] = phi[i, jg]
            end
            for j in 1:jg
                phi[0, j] = phi[1, j]
                phi[ig+1, j] = phi[ig, j]
            end
        end

        if ermax < eps
            ncount_final = nc_it
            break
        end
    end

    ncount = ncount_final

    p_stn = p[1, 1] + phi[1, 1]
    Threads.@threads for j in 1:jg
        for i in 1:ig
            u[i, j, 1] = ut[i, j, 1] - dt * (phi[i+1, j] - phi[i, j]) * pdxx[i+1]
            u[i, j, 2] = ut[i, j, 2] - dt * (phi[i, j+1] - phi[i, j]) * pdyy[j+1]
        end
    end
    Threads.@threads for j in 1:jg
        for i in 1:ig
            p[i, j] = p[i, j] + phi[i, j] - p_stn
        end
    end

    return ncount
end

# ----------------------------------------------------------------------
# puv_field????.dat  (8E13.5 per line)
# ----------------------------------------------------------------------
function write_puv(s, num)
    ig, jg = s.IG, s.JG
    u, p, beta, xx, yy, dx, dy = s.U, s.P, s.BETA, s.XX, s.YY, s.DX, s.DY
    nstep = @sprintf("%04d", num % 10000)
    open("TEXTDATA/puv_field" * nstep * ".dat", "w") do f20
        for i in 1:ig, j in 1:jg
            vort = -((u[i, j, 1] + u[i-1, j, 1]) * 0.5 - (u[i, j-1, 1] + u[i-1, j-1, 1]) * 0.5) / dy[j] +
                ((u[i, j, 2] + u[i, j-1, 2]) * 0.5 - (u[i-1, j, 2] + u[i-1, j-1, 2]) * 0.5) / dx[i]
            println(f20,
                fE(xx[i], 13, 5), fE(yy[j], 13, 5), fE(p[i, j], 13, 5),
                fE((u[i, j, 1] + u[i-1, j, 1]) * 0.5, 13, 5),
                fE((u[i, j, 2] + u[i, j-1, 2]) * 0.5, 13, 5),
                fE((beta[i, j, 2] + beta[i, j-1, 2]) * 0.5, 13, 5),
                fE((beta[i, j, 1] + beta[i-1, j, 1]) * 0.5, 13, 5),
                fE(vort, 13, 5))
        end
    end
end

# ----------------------------------------------------------------------
# UBC : boundary conditions, DIV / CFL / mean U, C_L, C_D
# ----------------------------------------------------------------------
function ubc!(s, num)
    ig, jg, dt, yl = s.IG, s.JG, s.DT, s.YL
    u, p, beta, pdx, pdy, dx, dy = s.U, s.P, s.BETA, s.PDX, s.PDY, s.DX, s.DY

    for i in 1:ig
        u[i, 0, 1] = u[i, 1, 1]
        u[i, jg+1, 1] = u[i, jg, 1]
        u[i, 0, 2] = 0.0
        u[i, jg, 2] = 0.0
    end
    for j in 1:jg
        u[0, j, 1] = 1.0
        u[ig+1, j, 1] = u[ig, j, 1]
        u[0, j, 2] = -u[1, j, 2]
        u[ig+1, j, 2] = -u[ig, j, 2]
    end
    for i in -3:(ig+3)
        p[i, 0] = p[i, 1]
        p[i, jg+1] = p[i, jg]
    end
    for j in 1:jg
        p[0, j] = p[1, j]
        p[ig+1, j] = p[ig, j]
    end

    div = 0.0
    for j in 1:jg, i in 1:ig
        qqqq = (u[i, j, 1] - u[i-1, j, 1]) * pdx[i] + (u[i, j, 2] - u[i, j-1, 2]) * pdy[j]
        qqq = abs(qqqq)
        if qqq > div
            div = qqq
        end
    end

    cfl = 0.0
    for j in 1:jg, i in 1:ig
        cfll = (abs(u[i, j, 1] + u[i-1, j, 1]) * pdx[i] +
            abs(u[i, j, 2] + u[i, j-1, 2]) * pdy[j]) * 0.5 * dt
        if cfll > cfl
            cfl = cfll
        end
    end

    usum = 0.0
    for j in 1:jg, i in 1:ig
        usum = usum + u[i, j, 1] * dy[j]
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
            if beta[i, j, 2] - beta[i, j-1, 2] >= 0.50
                p_top = p[i, j+1]
            end
            if beta[i, j, 2] - beta[i, j-1, 2] <= -0.50
                p_btm = p[i, j-2]
            end
        end
        p_sum = p_sum + dx[i] * (p_btm - p_top)
    end
    cl = p_sum * 2.0

    p_sum = 0.0
    for j in 1:jg
        p_top = 0.0
        p_btm = 0.0
        for i in 2:ig
            if beta[i, j, 1] - beta[i-1, j, 1] >= 0.5
                p_top = p[i+3, j]
            end
            if beta[i, j, 1] - beta[i-1, j, 1] <= -0.5
                p_btm = p[i-4, j]
            end
        end
        p_sum = p_sum + dy[j] * (p_btm - p_top)
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
            println(f14, @sprintf("%3d", i), fE(x[i], 14, 6), fE(dx[i], 14, 6),
                fE(xx[i], 14, 6), fE(dxx[i], 14, 6))
        end

        println(f14, a71(cprm5))
        for j in 0:(ly+1)
            println(f14, @sprintf("%3d", j), fE(y[j], 14, 6), fE(dy[j], 14, 6),
                fE(yy[j], 14, 6), fE(dyy[j], 14, 6))
        end

        println(f14, a71(cprm6))
        p_sum = 0.0
        for i in 1:ig
            p_top = 0.0
            p_btm = 0.0
            for j in 2:jg
                if beta[i, j, 2] - beta[i, j-1, 2] >= 0.50
                    p_top = p[i, j+1]
                end
                if beta[i, j, 2] - beta[i, j-1, 2] <= -0.50
                    p_btm = p[i, j-2]
                end
            end
            if p_top != 0.0
                println(f14, fE(x[i], 14, 6), fE(dx[i], 14, 6), fE(p_top, 14, 6), fE(p_btm, 14, 6))
            end
            p_sum = p_sum + dx[i] * (p_btm - p_top)
        end
        cl = p_sum * 2.0

        println(f14, a71(cprm7))
        p_sum = 0.0
        for j in 1:jg
            p_top = 0.0
            p_btm = 0.0
            for i in 2:ig
                if beta[i, j, 1] - beta[i-1, j, 1] >= 0.5
                    p_top = p[i+3, j]
                end
                if beta[i, j, 1] - beta[i-1, j, 1] <= -0.5
                    p_btm = p[i-4, j]
                end
            end
            if p_top != 0.0
                println(f14, fE(y[j], 14, 6), fE(dy[j], 14, 6), fE(p_top, 14, 6), fE(p_btm, 14, 6))
            end
            p_sum = p_sum + dy[j] * (p_btm - p_top)
        end
        cd = p_sum * 2.0
    end

    @printf(" C_L =%6.3f , C_D =%6.3f\n", cl, cd)
end

# ----------------------------------------------------------------------
# SEED DATA : restart files (full precision so a restart is loss-free)
# ----------------------------------------------------------------------
function seedout(s, num, time)
    open("seed_u.dat", "w") do f31
        @printf(f31, "%d %.16e\n", num, time)
        @printf(f31, "%d %d %.16e %.16e %.16e %.16e\n", s.IG, s.JG, s.XL, s.YL, s.DT, s.RE)
        for v in parent(s.U)                 # column-major = Fortran memory order
            @printf(f31, "%.16e\n", v)
        end
    end
    open("seed_p.dat", "w") do f32
        for v in parent(s.P)
            @printf(f32, "%.16e\n", v)
        end
    end
    open("seed_old.dat", "w") do f33
        for v in parent(s.CNX_OLD)
            @printf(f33, "%.16e\n", v)
        end
        for v in parent(s.CNY_OLD)
            @printf(f33, "%.16e\n", v)
        end
    end
end

# Same as Fortran SEEDIN (not called by the main program, as in the Fortran code).
# Returns (IERR, NUM, TIME, DT, RE).
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
    copyto!(vec(parent(s.U)), parse_f.(tu[9:(8+nu)]))
    copyto!(vec(parent(s.P)), parse_f.(tokens("seed_p.dat")[1:length(s.P)]))
    to = tokens("seed_old.dat")
    nx = length(s.CNX_OLD)
    ny = length(s.CNY_OLD)
    copyto!(vec(parent(s.CNX_OLD)), parse_f.(to[1:nx]))
    copyto!(vec(parent(s.CNY_OLD)), parse_f.(to[(nx+1):(nx+ny)]))

    println(" READ VELOCITY-SEED DATA OK!")
    return 0, num, time, dt, re
end

# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
function main()
    p = read_setting("setting.dat")
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