# =======================================================================
# SIMULATION OF TWO DIMENSIONAL FLOW PAST AN AIRFOIL
# BY FINITE DIFFERENCE METHOD (Dynamic Grid & Setting File Version)
#
# Julia port of simulation_airfoil_2_3.f
#   * Arrays have the SAME index ranges as the Fortran code
#     (OffsetArrays), so each expression is a literal transcription.
#   * Setup : julia> using Pkg; Pkg.add("OffsetArrays")
#   * Run   : julia -t auto simulation_airfoil.jl
#     (results do not depend on the number of threads)
# =======================================================================

using Printf
using OffsetArrays

const PI_F = 3.141592653589793          # REAL*8, PARAMETER :: PI

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
# zero-filled array with Fortran-style index ranges, e.g. oa(Float64, 0:IG+1, -2:JG+2)
oa(::Type{T}, rs::AbstractUnitRange...) where {T} =
    OffsetArray(zeros(T, map(length, rs)...), rs...)

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
ts_line(NUM, TIME, DIV, CFL, USUM, CL, CD, NCOUNT) =
    @sprintf("%6d%8d%10.4f", NUM, NUM, TIME) *
        fE(DIV, 12, 4) * fE(CFL, 12, 4) * fE(USUM, 12, 4) *
        @sprintf("%9.4f%9.4f%8d", CL, CD, NCOUNT)

# ----------------------------------------------------------------------
# setting.dat   (READ(50,*) value, comment  -> first token of each line)
# ----------------------------------------------------------------------
function read_setting(file::AbstractString)
    lines = readlines(file)
    IG = parse(Int, first_token(lines[1]))
    JG = parse(Int, first_token(lines[2]))
    IGS = parse(Int, first_token(lines[3]))     # overwritten by auto-detection
    RE = parse_f(first_token(lines[4]))
    DEG = parse_f(first_token(lines[5]))
    DT = parse_f(first_token(lines[6]))
    NC = parse(Int, first_token(lines[7]))
    NG = parse(Int, first_token(lines[8]))
    XL = parse_f(first_token(lines[9]))
    YL = parse_f(first_token(lines[10]))
    NWING = String(strip(first_token(lines[11]), ['\'', '"']))
    return (; IG, JG, IGS, RE, DEG, DT, NC, NG, XL, YL, NWING)
end

# ----------------------------------------------------------------------
# allocation (same bounds as the ALLOCATE statements in the Fortran code)
# ----------------------------------------------------------------------
function allocate_sim(p)
    IG, JG = p.IG, p.JG
    LX, LY = IG + 1, JG + 1

    # BETA/GAMMA are declared (0:IG,0:JG) in Fortran but FS0 touches
    # IG+1 / JG+1 as well; those extra cells are 1.0 (free stream) here.
    BETA = oa(Float64, 0:(IG+1), 0:(JG+1), 1:2)
    fill!(BETA, 1.0)
    GAMMA = oa(Float64, 0:(IG+1), 0:(JG+1), 1:2)

    return (;
        IG, JG, LX, LY,
        RE=p.RE, DEG=p.DEG, DT=p.DT, NC=p.NC, NG=p.NG,
        XL=p.XL, YL=p.YL, NWING=p.NWING,
        U=oa(Float64, 0:(IG+1), 0:(JG+1), 1:2),
        UT=oa(Float64, 0:(IG+1), 0:(JG+1), 1:2),
        CNX_OLD=oa(Float64, 0:IG, 0:JG),
        CNY_OLD=oa(Float64, 0:IG, 0:JG),
        P=oa(Float64, -4:(IG+4), -2:(JG+2)),
        PHI=oa(Float64, -4:(IG+4), -2:(JG+2)),      # PHI = 0 (warm start later)
        BETA, GAMMA,
        IDX=oa(Int, 0:IG, 0:JG),
        IDY=oa(Int, 0:IG, 0:JG),
        X=oa(Float64, -3:(LX+5)), XX=oa(Float64, -3:(LX+5)),
        DX=oa(Float64, -3:(LX+5)), DXX=oa(Float64, -3:(LX+5)),
        Y=oa(Float64, -3:(LY+5)), YY=oa(Float64, -3:(LY+5)),
        DY=oa(Float64, -3:(LY+5)), DYY=oa(Float64, -3:(LY+5)),
        PDX=oa(Float64, 0:(LX+2)), PDXX=oa(Float64, 0:(LX+2)),
        DX2=oa(Float64, 0:(LX+2)), DXX2=oa(Float64, 0:(LX+2)),
        PDX2=oa(Float64, 0:(LX+2)), PDX3=oa(Float64, 0:(LX+2)),
        PDXX3=oa(Float64, 0:(LX+2)),
        PDY=oa(Float64, 0:(LY+2)), PDYY=oa(Float64, 0:(LY+2)),
        DY2=oa(Float64, 0:(LY+2)), DYY2=oa(Float64, 0:(LY+2)),
        PDY2=oa(Float64, 0:(LY+2)), PDY3=oa(Float64, 0:(LY+2)),
        PDYY3=oa(Float64, 0:(LY+2)),
        D_N=oa(Float64, 0:JG), D_E=oa(Float64, 0:IG),
        D_W=oa(Float64, 0:IG), D_S=oa(Float64, 0:JG),
        PDPHI=oa(Float64, 0:IG, 0:JG),
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
function calbeta!(S)
    IG, JG, LX, LY = S.IG, S.JG, S.LX, S.LY
    XL, YL, DEG, DT, RE = S.XL, S.YL, S.DEG, S.DT, S.RE
    X, XX, DX, DXX = S.X, S.XX, S.DX, S.DXX
    Y, YY, DY, DYY = S.Y, S.YY, S.DY, S.DYY
    BETA, GAMMA, IDX, IDY = S.BETA, S.GAMMA, S.IDX, S.IDY
    PDX, PDXX, DX2, DXX2, PDX2, PDX3, PDXX3 = S.PDX, S.PDXX, S.DX2, S.DXX2, S.PDX2, S.PDX3, S.PDXX3
    PDY, PDYY, DY2, DYY2, PDY2, PDY3, PDYY3 = S.PDY, S.PDYY, S.DY2, S.DYY2, S.PDY2, S.PDY3, S.PDYY3
    D_N, D_E, D_W, D_S, PDPHI = S.D_N, S.D_E, S.D_W, S.D_S, S.PDPHI

    NWING_FILE = S.NWING
    if !isfile(NWING_FILE)
        println(" ERROR: Cannot open wing file: ", NWING_FILE)
        return 1
    end

    xt, yt, yb = read_wing(NWING_FILE)
    IGS = length(xt)                              # auto-detected
    @printf(" Detected wing data points (IGS) = %11d\n", IGS)

    WNG_XT = zeros(IGS + 1)
    WNG_XB = zeros(IGS + 1)
    WNG_YT = zeros(IGS + 1)
    WNG_YB = zeros(IGS + 1)
    WNG_XT[1:IGS] .= xt
    WNG_YT[1:IGS] .= yt
    WNG_YB[1:IGS] .= yb

    X_TEMP = oa(Float64, -3:(LX+5))
    Y_TEMP = oa(Float64, -3:(LY+5))

    THETA = PI_F * DEG / 180.0

    for II in 1:IGS
        TWNG_XT = WNG_XT[II] * cos(-THETA) - WNG_YT[II] * sin(-THETA)
        TWNG_YT = WNG_XT[II] * sin(-THETA) + WNG_YT[II] * cos(-THETA)
        TWNG_XB = WNG_XT[II] * cos(-THETA) - WNG_YB[II] * sin(-THETA)
        TWNG_YB = WNG_XT[II] * sin(-THETA) + WNG_YB[II] * cos(-THETA)

        WNG_XT[II] = TWNG_XT + XL * 0.5 - cos(THETA) * 0.5
        WNG_XB[II] = TWNG_XB + XL * 0.5 - cos(THETA) * 0.5
        WNG_YT[II] = TWNG_YT + YL * 0.5 + sin(THETA) * 0.5 - 0.06
        WNG_YB[II] = TWNG_YB + YL * 0.5 + sin(THETA) * 0.5 - 0.06
    end

    XALP1 = 1.0
    XALP2 = 3.0
    for I in 1:IG
        X1 = 2.0 * Float64(I) / Float64(IG) - 1.0
        X_TEMP[I] = XL * (((XALP1 * (X1^XALP2) + X1) * 0.5) / (XALP1 + 1.0) + 0.5)
    end
    YALP1 = 1.0
    YALP2 = 5.0
    for J in 1:JG
        Y1 = 2.0 * Float64(J) / Float64(JG) - 1.0
        Y_TEMP[J] = YL * (((YALP1 * (Y1^YALP2) + Y1) * 0.5) / (YALP1 + 1.0) + 0.5)
    end

    X_TEMP[0] = 0.0
    X_TEMP[-1] = -X_TEMP[1]
    X_TEMP[-2] = -X_TEMP[2]
    X_TEMP[-3] = -X_TEMP[3]
    X_TEMP[IG+1] = 2.0 * XL - X_TEMP[IG-1]
    X_TEMP[IG+2] = 2.0 * XL - X_TEMP[IG-2]
    X_TEMP[IG+3] = 2.0 * XL - X_TEMP[IG-3]
    X_TEMP[IG+4] = 2.0 * XL - X_TEMP[IG-4]
    X_TEMP[IG+5] = 2.0 * XL - X_TEMP[IG-5]

    Y_TEMP[0] = 0.0
    Y_TEMP[-1] = -Y_TEMP[1]
    Y_TEMP[-2] = -Y_TEMP[2]
    Y_TEMP[-3] = -Y_TEMP[3]
    Y_TEMP[JG+1] = 2.0 * YL - Y_TEMP[JG-1]
    Y_TEMP[JG+2] = 2.0 * YL - Y_TEMP[JG-2]
    Y_TEMP[JG+3] = 2.0 * YL - Y_TEMP[JG-3]
    Y_TEMP[JG+4] = 2.0 * YL - Y_TEMP[JG-4]
    Y_TEMP[JG+5] = 2.0 * YL - Y_TEMP[JG-5]

    for I in -2:(LX+5)
        X[I] = X_TEMP[I]
        DX[I] = X[I] - X[I-1]
        XX[I] = (X[I] + X[I-1]) * 0.5
        DXX[I] = XX[I] - XX[I-1]
    end
    for J in -2:(LY+5)
        Y[J] = Y_TEMP[J]
        DY[J] = Y[J] - Y[J-1]
        YY[J] = (Y[J] + Y[J-1]) * 0.5
        DYY[J] = YY[J] - YY[J-1]
    end

    for I in 0:IG, J in 0:JG
        BETA[I, J, 1] = 1.0
        BETA[I, J, 2] = 1.0
    end

    # --- BETA(:,:,1)  (u-points) and GAMMA ---
    for I in 0:IG
        IB = 0
        IT = 0
        for II in 1:IGS
            if (X[I] < WNG_XB[II+1]) && (X[I] >= WNG_XB[II]); IB = II; end
            if (X[I] < WNG_XT[II+1]) && (X[I] >= WNG_XT[II]); IT = II; end
        end
        (IB * IT == 0) && continue
        XTA = X[I] - WNG_XT[IT]
        XBA = X[I] - WNG_XB[IB]
        XTB = WNG_XT[IT+1] - X[I]
        XBB = WNG_XB[IB+1] - X[I]
        YT = (WNG_YT[IT] * XTB + WNG_YT[IT+1] * XTA) / (XTA + XTB)
        YB = (WNG_YB[IB] * XBB + WNG_YB[IB+1] * XBA) / (XBA + XBB)
        for J in 0:JG
            if (YY[J] >= YB) && (YY[J] < YT)
                BETA[I, J, 1] = 0.0
            end
            if (YY[J] <= YB) && (YY[J+1] > YB)
                GAMMA[I, J, 1] = (YB - YY[J]) / (YB - YY[J-1])
            end
            if (YY[J] >= YT) && (YY[J-1] < YT)
                GAMMA[I, J, 1] = (YY[J] - YT) / (YY[J+1] - YT)
            end
        end
    end

    # --- BETA(:,:,2)  (v-points) ---
    for I in 0:IG
        IB = 0
        IT = 0
        for II in 1:IGS
            if (XX[I] < WNG_XB[II+1]) && (XX[I] >= WNG_XB[II]); IB = II; end
            if (XX[I] < WNG_XT[II+1]) && (XX[I] >= WNG_XT[II]); IT = II; end
        end
        (IB * IT == 0) && continue
        XTA = XX[I] - WNG_XT[IT]
        XBA = XX[I] - WNG_XB[IB]
        XTB = WNG_XT[IT+1] - XX[I]
        XBB = WNG_XB[IB+1] - XX[I]
        YT = (WNG_YT[IT] * XTB + WNG_YT[IT+1] * XTA) / (XTA + XTB)
        YB = (WNG_YB[IB] * XBB + WNG_YB[IB+1] * XBA) / (XBA + XBB)
        for J in 0:JG
            if (Y[J] >= YB) && (Y[J] < YT)
                BETA[I, J, 2] = 0.0
            end
        end
    end

    println(" GRID GENERATION OK!")

    for I in 1:IG, J in 1:(JG-1)
        IDX[I, J] = 4 * Int(BETA[I, J-1, 1]) + 2 * Int(BETA[I, J, 1]) + 1 * Int(BETA[I, J+1, 1])
        IDY[I, J] = 4 * Int(BETA[I-1, J, 2]) + 2 * Int(BETA[I, J, 2]) + 1 * Int(BETA[I+1, J, 2])
    end

    for I in 0:(IG+1)
        DX2[I] = DX[I] + DX[I+1]
        DXX2[I] = DXX[I] + DXX[I+1]
        PDX[I] = 1.0 / DX[I]
        PDXX[I] = 1.0 / DXX[I]
        PDX2[I] = 1.0 / (DX[I] + DX[I+1])
        PDX3[I] = 1.0 / (DX[I] * DX[I+1] * (DX[I] + DX[I+1]))
        PDXX3[I] = 1.0 / (DXX[I] * DXX[I+1] * (DXX[I] + DXX[I+1]))
    end
    for J in 0:(JG+1)
        DY2[J] = DY[J] + DY[J+1]
        DYY2[J] = DYY[J] + DYY[J+1]
        PDY[J] = 1.0 / DY[J]
        PDYY[J] = 1.0 / DYY[J]
        PDY2[J] = 1.0 / (DY[J] + DY[J+1])
        PDY3[J] = 1.0 / (DY[J] * DY[J+1] * (DY[J] + DY[J+1]))
        PDYY3[J] = 1.0 / (DYY[J] * DYY[J+1] * (DYY[J] + DYY[J+1]))
    end

    for I in 1:IG
        D_W[I] = PDX[I] * PDXX[I]
        D_E[I] = PDX[I] * PDXX[I+1]
    end
    for J in 1:JG
        D_S[J] = PDY[J] * PDYY[J]
        D_N[J] = PDY[J] * PDYY[J+1]
    end
    for I in 1:IG, J in 1:JG
        PDPHI[I, J] = 1.0 / (D_W[I] + D_E[I] + D_S[J] + D_N[J])
    end

    DFN = 0.0
    for I in 1:IG, J in 1:JG
        DFNN = DT / RE * (1.0 / (DX[I]^2) + 1.0 / (DY[J]^2))
        DFN = max(DFN, DFNN)
    end
    println(lpad("DIFFUSION NUMBER = ", 20), fE(DFN, 11, 4))

    return 0
end

# ----------------------------------------------------------------------
# FS0 / FS1 : temporary velocity
#   FS0 : Euler explicit (first step)       FS1 : Adams-Bashforth (2nd step on)
# The two Fortran routines differ only in the time integration, so they
# share one body here (identical arithmetic).
# ----------------------------------------------------------------------
function fs_predictor!(S, first_step::Bool)
    IG, JG, RE, DT = S.IG, S.JG, S.RE, S.DT
    U, UT, P, BETA, GAMMA, IDX = S.U, S.UT, S.P, S.BETA, S.GAMMA, S.IDX
    CNX_OLD, CNY_OLD = S.CNX_OLD, S.CNY_OLD
    PDX, PDXX, DX, DX2, PDX2, PDX3, PDXX3, DXX, DXX2 = S.PDX, S.PDXX, S.DX, S.DX2, S.PDX2, S.PDX3, S.PDXX3, S.DXX, S.DXX2
    PDY, PDYY, DY, DY2, PDY2, PDY3, DYY, DYY2, PDYY3 = S.PDY, S.PDYY, S.DY, S.DY2, S.PDY2, S.PDY3, S.DYY, S.DYY2, S.PDYY3

    PRE = 1.0 / RE

    # ---- u-momentum ----
    Threads.@threads for J in 1:JG
        for I in 1:IG
            if BETA[I, J, 1] == 1.0
                UDUDX = 0.5 * ((U[I+1, J, 1] + U[I, J, 1]) * (U[I+1, J, 1] - U[I, J, 1]) * PDX[I+1] * DX[I] +
                        (U[I, J, 1] + U[I-1, J, 1]) * (U[I, J, 1] - U[I-1, J, 1]) * PDX[I] * DX[I+1]) * PDX2[I]
                VDUDY = 0.5 * ((U[I+1, J, 2] * DX[I] + U[I, J, 2] * DX[I+1]) * (U[I, J+1, 1] - U[I, J, 1]) * PDYY[J+1] +
                        (U[I+1, J-1, 2] * DX[I] + U[I, J-1, 2] * DX[I+1]) * (U[I, J, 1] - U[I, J-1, 1]) * PDYY[J]) * PDX2[I]
                U_CON = UDUDX + VDUDY

                UVISX = (DX[I+1] * U[I-1, J, 1] - DX2[I] * U[I, J, 1] + DX[I] * U[I+1, J, 1]) * 2.0 * PDX3[I] * PRE
                UVISY = (DYY[J+1] * U[I, J-1, 1] - DYY2[J] * U[I, J, 1] + DYY[J] * U[I, J+1, 1]) * 2.0 * PDYY3[J] * PRE
                U_VIS = UVISX + UVISY

                DPDX = (P[I+1, J] - P[I, J]) * PDXX[I+1]

                if first_step
                    UT[I, J, 1] = U[I, J, 1] - DT * (U_CON - U_VIS + DPDX)
                else
                    UT[I, J, 1] = (U[I, J, 1] - DT * ((U_CON - U_VIS + DPDX) * 1.5 - CNX_OLD[I, J] * 0.5)) * BETA[I, J, 1]
                end
                CNX_OLD[I, J] = U_CON - U_VIS + DPDX
            end
        end
    end

    # ---- v-momentum ----
    Threads.@threads for J in 1:(JG-1)
        for I in 1:IG
            if BETA[I, J, 2] == 1.0
                UDVDX = 0.5 * ((U[I, J+1, 1] * DY[J] + U[I, J, 1] * DY[J+1]) * (U[I+1, J, 2] - U[I, J, 2]) * PDXX[I+1] +
                        (U[I-1, J+1, 1] * DY[J] + U[I-1, J, 1] * DY[J+1]) * (U[I, J, 2] - U[I-1, J, 2]) * PDXX[I]) * PDY2[J]
                VDVDY = 0.5 * ((U[I, J+1, 2] + U[I, J, 2]) * (U[I, J+1, 2] - U[I, J, 2]) * PDY[J+1] * DY[J] +
                        (U[I, J, 2] + U[I, J-1, 2]) * (U[I, J, 2] - U[I, J-1, 2]) * PDY[J] * DY[J+1]) * PDY2[J]
                V_CON = UDVDX + VDVDY

                VVISX = (DXX[I+1] * U[I-1, J, 2] - DXX2[I] * U[I, J, 2] + DXX[I] * U[I+1, J, 2]) * 2.0 * PDXX3[I] * PRE
                VVISY = (DY[J+1] * U[I, J-1, 2] - DY2[J] * U[I, J, 2] + DY[J] * U[I, J+1, 2]) * 2.0 * PDY3[J] * PRE
                V_VIS = VVISX + VVISY

                DPDY = (P[I, J+1] - P[I, J]) * PDYY[J+1]

                if first_step
                    UT[I, J, 2] = (U[I, J, 2] - DT * (V_CON - V_VIS + DPDY)) * BETA[I, J, 2]
                else
                    UT[I, J, 2] = (U[I, J, 2] - DT * ((V_CON - V_VIS + DPDY) * 1.5 - CNY_OLD[I, J] * 0.5)) * BETA[I, J, 2]
                end
                CNY_OLD[I, J] = V_CON - V_VIS + DPDY
            end
        end
    end

    # ---- cut-cell correction near the wing surface ----
    Threads.@threads for I in 1:IG
        for J in 1:JG
            if IDX[I, J] == 3
                UT[I, J, 1] = UT[I, J+1, 1] * GAMMA[I, J, 1]
            end
            if IDX[I, J] == 5
                UT[I, J, 1] = UT[I, J-1, 1] * GAMMA[I, J, 1]
            end
        end
    end
    return nothing
end

function fs0!(S)
    IG, JG = S.IG, S.JG
    U, P, BETA = S.U, S.P, S.BETA
    Threads.@threads for J in 0:(JG+1)
        for I in 0:(IG+1)
            U[I, J, 1] = 1.0 * BETA[I, J, 1]
            U[I, J, 2] = 0.0 * BETA[I, J, 2]
            P[I, J] = 0.0
        end
    end
    fs_predictor!(S, true)
end

fs1!(S) = fs_predictor!(S, false)

# ----------------------------------------------------------------------
# FS2 : Poisson equation for pressure correction (red-black SOR)
# ----------------------------------------------------------------------
function sor_row(PHI, UT, D_W, D_E, D_S, D_N, PDX, PDY, PDPHI, PDT, OMG, IG, J, COLOR)
    J_START = 1 + ((J + COLOR) & 1)
    m = 0.0
    for I in J_START:2:IG
        PHI_CALC = (D_W[I] * PHI[I-1, J] + D_E[I] * PHI[I+1, J] +
            D_S[J] * PHI[I, J-1] + D_N[J] * PHI[I, J+1] -
            PDT * ((UT[I, J, 1] - UT[I-1, J, 1]) * PDX[I] +
            (UT[I, J, 2] - UT[I, J-1, 2]) * PDY[J])) * PDPHI[I, J]
        PHIDIF = PHI_CALC - PHI[I, J]
        PHI[I, J] = PHI[I, J] + OMG * PHIDIF
        a = abs(PHIDIF)
        if a > m
            m = a
        end
    end
    return m
end

function fs2!(S)
    IG, JG, DT, XL = S.IG, S.JG, S.DT, S.XL
    U, UT, P, PHI, X = S.U, S.UT, S.P, S.PHI, S.X
    PDX, PDY, PDXX, PDYY, PDPHI = S.PDX, S.PDY, S.PDXX, S.PDYY, S.PDPHI
    D_W, D_E, D_S, D_N = S.D_W, S.D_E, S.D_S, S.D_N

    IGF = 12

    # outflow damping zone
    Threads.@threads for J in 1:JG
        for I in (div(IG*(IGF-1), IGF)+1):IG
            UT[I, J, 1] = (UT[I, J, 1] - 1.0) * (-tanh(4.5 * (X[I] - XL))) + 1.0
            UT[I, J, 2] = UT[I, J, 2] * (-tanh(4.5 * (X[I] - XL)))
        end
    end

    for I in 1:IG
        UT[I, 0, 1] = UT[I, 1, 1]
        UT[I, JG+1, 1] = UT[I, JG, 1]
        UT[I, 0, 2] = 0.0
        UT[I, JG, 2] = 0.0
    end
    for J in 1:JG
        UT[0, J, 1] = 1.0
        UT[IG+1, J, 1] = UT[IG, J, 1]
        UT[0, J, 2] = -UT[1, J, 2]
        UT[IG+1, J, 2] = -UT[IG, J, 2]
    end

    PDT = 1.0 / DT
    OMG = 1.4
    EPS = Float64(1.0f-4)          # Fortran "1.0E-4" is a single-precision literal

    NCOUNT_FINAL = 50001           # value left by a DO loop that never converges
    rowmax = zeros(Float64, JG)
    multi = Threads.nthreads() > 1

    for NC_IT in 1:50000
        ERMAX = 0.0
        for COLOR in 0:1
            if multi
                Threads.@threads for J in 1:JG
                    rowmax[J] = sor_row(PHI, UT, D_W, D_E, D_S, D_N, PDX, PDY, PDPHI, PDT, OMG, IG, J, COLOR)
                end
            else
                for J in 1:JG
                    rowmax[J] = sor_row(PHI, UT, D_W, D_E, D_S, D_N, PDX, PDY, PDPHI, PDT, OMG, IG, J, COLOR)
                end
            end
            ERMAX_COLOR = maximum(rowmax)
            if ERMAX_COLOR > ERMAX
                ERMAX = ERMAX_COLOR
            end

            for I in 1:IG
                PHI[I, 0] = PHI[I, 1]
                PHI[I, JG+1] = PHI[I, JG]
            end
            for J in 1:JG
                PHI[0, J] = PHI[1, J]
                PHI[IG+1, J] = PHI[IG, J]
            end
        end

        if ERMAX < EPS
            NCOUNT_FINAL = NC_IT
            break
        end
    end

    NCOUNT = NCOUNT_FINAL

    P_STN = P[1, 1] + PHI[1, 1]
    Threads.@threads for J in 1:JG
        for I in 1:IG
            U[I, J, 1] = UT[I, J, 1] - DT * (PHI[I+1, J] - PHI[I, J]) * PDXX[I+1]
            U[I, J, 2] = UT[I, J, 2] - DT * (PHI[I, J+1] - PHI[I, J]) * PDYY[J+1]
        end
    end
    Threads.@threads for J in 1:JG
        for I in 1:IG
            P[I, J] = P[I, J] + PHI[I, J] - P_STN
        end
    end

    return NCOUNT
end

# ----------------------------------------------------------------------
# puv_field????.dat  (8E13.5 per line)
# ----------------------------------------------------------------------
function write_puv(S, NUM)
    IG, JG = S.IG, S.JG
    U, P, BETA, XX, YY, DX, DY = S.U, S.P, S.BETA, S.XX, S.YY, S.DX, S.DY
    NSTEP = @sprintf("%04d", NUM % 10000)
    open("TEXTDATA/puv_field" * NSTEP * ".dat", "w") do f20
        for I in 1:IG, J in 1:JG
            vort = -((U[I, J, 1] + U[I-1, J, 1]) * 0.5 - (U[I, J-1, 1] + U[I-1, J-1, 1]) * 0.5) / DY[J] +
                ((U[I, J, 2] + U[I, J-1, 2]) * 0.5 - (U[I-1, J, 2] + U[I-1, J-1, 2]) * 0.5) / DX[I]
            println(f20,
                fE(XX[I], 13, 5), fE(YY[J], 13, 5), fE(P[I, J], 13, 5),
                fE((U[I, J, 1] + U[I-1, J, 1]) * 0.5, 13, 5),
                fE((U[I, J, 2] + U[I, J-1, 2]) * 0.5, 13, 5),
                fE((BETA[I, J, 2] + BETA[I, J-1, 2]) * 0.5, 13, 5),
                fE((BETA[I, J, 1] + BETA[I-1, J, 1]) * 0.5, 13, 5),
                fE(vort, 13, 5))
        end
    end
end

# ----------------------------------------------------------------------
# UBC : boundary conditions, DIV / CFL / mean U, C_L, C_D
# ----------------------------------------------------------------------
function ubc!(S, NUM)
    IG, JG, DT, YL = S.IG, S.JG, S.DT, S.YL
    U, P, BETA, PDX, PDY, DX, DY = S.U, S.P, S.BETA, S.PDX, S.PDY, S.DX, S.DY

    for I in 1:IG
        U[I, 0, 1] = U[I, 1, 1]
        U[I, JG+1, 1] = U[I, JG, 1]
        U[I, 0, 2] = 0.0
        U[I, JG, 2] = 0.0
    end
    for J in 1:JG
        U[0, J, 1] = 1.0
        U[IG+1, J, 1] = U[IG, J, 1]
        U[0, J, 2] = -U[1, J, 2]
        U[IG+1, J, 2] = -U[IG, J, 2]
    end
    for I in -3:(IG+3)
        P[I, 0] = P[I, 1]
        P[I, JG+1] = P[I, JG]
    end
    for J in 1:JG
        P[0, J] = P[1, J]
        P[IG+1, J] = P[IG, J]
    end

    DIV = 0.0
    for J in 1:JG, I in 1:IG
        QQQQ = (U[I, J, 1] - U[I-1, J, 1]) * PDX[I] + (U[I, J, 2] - U[I, J-1, 2]) * PDY[J]
        QQQ = abs(QQQQ)
        if QQQ > DIV
            DIV = QQQ
        end
    end

    CFL = 0.0
    for J in 1:JG, I in 1:IG
        CFLL = (abs(U[I, J, 1] + U[I-1, J, 1]) * PDX[I] +
            abs(U[I, J, 2] + U[I, J-1, 2]) * PDY[J]) * 0.5 * DT
        if CFLL > CFL
            CFL = CFLL
        end
    end

    USUM = 0.0
    for J in 1:JG, I in 1:IG
        USUM = USUM + U[I, J, 1] * DY[J]
    end
    USUM = USUM / Float64(IG) / YL

    if NUM % 5 == 0
        write_puv(S, NUM)
    end

    P_SUM = 0.0
    for I in 1:IG
        P_TOP = 0.0
        P_BTM = 0.0
        for J in 2:JG
            if BETA[I, J, 2] - BETA[I, J-1, 2] >= 0.50
                P_TOP = P[I, J+1]
            end
            if BETA[I, J, 2] - BETA[I, J-1, 2] <= -0.50
                P_BTM = P[I, J-2]
            end
        end
        P_SUM = P_SUM + DX[I] * (P_BTM - P_TOP)
    end
    CL = P_SUM * 2.0

    P_SUM = 0.0
    for J in 1:JG
        P_TOP = 0.0
        P_BTM = 0.0
        for I in 2:IG
            if BETA[I, J, 1] - BETA[I-1, J, 1] >= 0.5
                P_TOP = P[I+3, J]
            end
            if BETA[I, J, 1] - BETA[I-1, J, 1] <= -0.5
                P_BTM = P[I-4, J]
            end
        end
        P_SUM = P_SUM + DY[J] * (P_BTM - P_TOP)
    end
    CD = P_SUM * 2.0

    return DIV, CFL, USUM, CL, CD
end

# ----------------------------------------------------------------------
# OUTPUT : TEXTDATA/output.dat
# ----------------------------------------------------------------------
function output(S, NUM, TIME)
    IG, JG, LY = S.IG, S.JG, S.LY
    P, BETA = S.P, S.BETA
    X, DX, XX, DXX = S.X, S.DX, S.XX, S.DXX
    Y, DY, YY, DYY = S.Y, S.DY, S.YY, S.DYY

    if NUM % 5 == 0
        write_puv(S, NUM)
    end

    CMTD1 = "#======================================================"
    CMTD2 = "# SIMULATION OF FLOW PAST AN 2D AIRFOIL : " * rpad(S.NWING, 64)
    CMTD3 = "# STATISTICS IN VELOCITY FIELD"
    CPRM1 = "# CALC.NUM.     CALC.TIME"
    CPRM2 = "# NX  NY  TIMESTEP  REYNOLDS"
    CPRM3 = "# X-LENGTH  Y-LENGTH AOA(DEG)"
    CPRM4 = "# I X DX XX DXX"
    CPRM5 = "# J Y DY YY DYY"
    CPRM6 = "# X DX P(UPPER-SURFACE)  P(BOTTOM SURFACE)"
    CPRM7 = "# Y DY P(DOWNSTREAM-SURFACE)  P(UPSTREAM SURFACE)"

    CL = 0.0
    CD = 0.0
    open("TEXTDATA/output.dat", "w") do f14
        println(f14, a71(CMTD1))
        println(f14, a71(CMTD2))
        println(f14, a71(CMTD1))
        println(f14, a71(CPRM1))
        println(f14, @sprintf("%10d", NUM), fE(TIME, 18, 10))
        println(f14, a71(CPRM2))
        println(f14, @sprintf("%4d%4d%11.6f%9.1f", IG, JG, S.DT, S.RE))
        println(f14, a71(CPRM3))
        println(f14, @sprintf("%10.6f%10.6f%7.2f", S.XL, S.YL, S.DEG))
        println(f14, a71(CMTD3))

        println(f14, a71(CPRM4))
        for I in 0:(IG+1)
            println(f14, @sprintf("%3d", I), fE(X[I], 14, 6), fE(DX[I], 14, 6),
                fE(XX[I], 14, 6), fE(DXX[I], 14, 6))
        end

        println(f14, a71(CPRM5))
        for J in 0:(LY+1)
            println(f14, @sprintf("%3d", J), fE(Y[J], 14, 6), fE(DY[J], 14, 6),
                fE(YY[J], 14, 6), fE(DYY[J], 14, 6))
        end

        println(f14, a71(CPRM6))
        P_SUM = 0.0
        for I in 1:IG
            P_TOP = 0.0
            P_BTM = 0.0
            for J in 2:JG
                if BETA[I, J, 2] - BETA[I, J-1, 2] >= 0.50
                    P_TOP = P[I, J+1]
                end
                if BETA[I, J, 2] - BETA[I, J-1, 2] <= -0.50
                    P_BTM = P[I, J-2]
                end
            end
            if P_TOP != 0.0
                println(f14, fE(X[I], 14, 6), fE(DX[I], 14, 6), fE(P_TOP, 14, 6), fE(P_BTM, 14, 6))
            end
            P_SUM = P_SUM + DX[I] * (P_BTM - P_TOP)
        end
        CL = P_SUM * 2.0

        println(f14, a71(CPRM7))
        P_SUM = 0.0
        for J in 1:JG
            P_TOP = 0.0
            P_BTM = 0.0
            for I in 2:IG
                if BETA[I, J, 1] - BETA[I-1, J, 1] >= 0.5
                    P_TOP = P[I+3, J]
                end
                if BETA[I, J, 1] - BETA[I-1, J, 1] <= -0.5
                    P_BTM = P[I-4, J]
                end
            end
            if P_TOP != 0.0
                println(f14, fE(Y[J], 14, 6), fE(DY[J], 14, 6), fE(P_TOP, 14, 6), fE(P_BTM, 14, 6))
            end
            P_SUM = P_SUM + DY[J] * (P_BTM - P_TOP)
        end
        CD = P_SUM * 2.0
    end

    @printf(" C_L =%6.3f , C_D =%6.3f\n", CL, CD)
end

# ----------------------------------------------------------------------
# SEED DATA : restart files (full precision so a restart is loss-free)
# ----------------------------------------------------------------------
function seedout(S, NUM, TIME)
    open("seed_u.dat", "w") do f31
        @printf(f31, "%d %.16e\n", NUM, TIME)
        @printf(f31, "%d %d %.16e %.16e %.16e %.16e\n", S.IG, S.JG, S.XL, S.YL, S.DT, S.RE)
        for v in parent(S.U)                 # column-major = Fortran memory order
            @printf(f31, "%.16e\n", v)
        end
    end
    open("seed_p.dat", "w") do f32
        for v in parent(S.P)
            @printf(f32, "%.16e\n", v)
        end
    end
    open("seed_old.dat", "w") do f33
        for v in parent(S.CNX_OLD)
            @printf(f33, "%.16e\n", v)
        end
        for v in parent(S.CNY_OLD)
            @printf(f33, "%.16e\n", v)
        end
    end
end

# Same as Fortran SEEDIN (not called by the main program, as in the Fortran code).
# Returns (IERR, NUM, TIME, DT, RE).
function seedin!(S)
    tokens(f) = split(read(f, String), r"[\s,]+"; keepempty=false)
    tu = tokens("seed_u.dat")
    NUM = parse(Int, tu[1])
    TIME = parse_f(tu[2])
    IGD = parse(Int, tu[3])
    JGD = parse(Int, tu[4])
    XLD = parse_f(tu[5])
    YLD = parse_f(tu[6])
    DT = parse_f(tu[7])
    RE = parse_f(tu[8])

    if IGD != S.IG || JGD != S.JG
        println(" MESH DATA ERROR")
        return 1000, NUM, TIME, DT, RE
    end
    if XLD != S.XL || YLD != S.YL
        println(" DOMAIN DATA ERROR")
        return 1000, NUM, TIME, DT, RE
    end

    nU = length(S.U)
    copyto!(vec(parent(S.U)), parse_f.(tu[9:(8+nU)]))
    copyto!(vec(parent(S.P)), parse_f.(tokens("seed_p.dat")[1:length(S.P)]))
    to = tokens("seed_old.dat")
    nx = length(S.CNX_OLD)
    ny = length(S.CNY_OLD)
    copyto!(vec(parent(S.CNX_OLD)), parse_f.(to[1:nx]))
    copyto!(vec(parent(S.CNY_OLD)), parse_f.(to[(nx+1):(nx+ny)]))

    println(" READ VELOCITY-SEED DATA OK!")
    return 0, NUM, TIME, DT, RE
end

# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
function main()
    p = read_setting("setting.dat")
    S = allocate_sim(p)
    DT, NC, NG = S.DT, S.NC, S.NG

    NUM = 0
    TIME = 0.0
    AVETIM = 0.0

    IERR = calbeta!(S)
    IERR != 0 && return

    mkpath("TEXTDATA")
    DIV, CFL, USUM, CL, CD = ubc!(S, NUM)

    NUM += 1
    TIME += DT
    AVETIM += DT

    NUM == 1 ? fs0!(S) : fs1!(S)
    NCOUNT = fs2!(S)
    DIV, CFL, USUM, CL, CD = ubc!(S, NUM)

    open("TEXTDATA/time-series.dat", "w") do f1
        println(f1, TS_HEADER)
        println(TS_HEADER)
        line = ts_line(NUM, TIME, DIV, CFL, USUM, CL, CD, NCOUNT)
        println(f1, line)
        println(line)

        for NST in 2:NC
            NUM += 1
            TIME += DT
            AVETIM += DT

            fs1!(S)
            NCOUNT = fs2!(S)
            DIV, CFL, USUM, CL, CD = ubc!(S, NUM)
            if NUM % NG == 0
                line = ts_line(NUM, TIME, DIV, CFL, USUM, CL, CD, NCOUNT)
                println(f1, line)
                println(line)
            end
        end
    end

    seedout(S, NUM, TIME)
    output(S, NUM, TIME)
    return nothing
end