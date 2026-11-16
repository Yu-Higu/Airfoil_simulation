function _fortran_tokens(path::AbstractString)
    vals = Float64[]
    for token in split(read(path, String))
        push!(vals, parse(Float64, replace(token, r"[dD]" => "E")))
    end
    return vals
end

function _take_values!(vals::Vector{Float64}, pos::Base.RefValue{Int}, n::Int)
    pos[] + n - 1 <= length(vals) || error("seed file のデータ数が不足しています")
    out = @view vals[pos[]:pos[]+n-1]
    pos[] += n
    return out
end

function SEEDIN!(d::CFDData; directory::AbstractString=".")
    upath = joinpath(directory, "seed_u.dat")
    ppath = joinpath(directory, "seed_p.dat")
    opath = joinpath(directory, "seed_old.dat")
    # Fortran list-directed input: whitespace/newlines are interchangeable.
    uvals = _fortran_tokens(upath)
    pvals = _fortran_tokens(ppath)
    ovals = _fortran_tokens(opath)
    pos = Ref(1)
    length(uvals) >= 8 || error("seed_u.dat header is incomplete")
    d.NUM = Int(round(uvals[pos[]])); pos[] += 1
    d.TIME = uvals[pos[]]; pos[] += 1
    igd = Int(round(uvals[pos[]])); pos[] += 1
    jgd = Int(round(uvals[pos[]])); pos[] += 1
    xld = uvals[pos[]]; pos[] += 1
    yld = uvals[pos[]]; pos[] += 1
    dtseed = uvals[pos[]]; pos[] += 1
    re_seed = uvals[pos[]]; pos[] += 1
    if igd != IG || jgd != JG
        println(" MESH DATA ERROR")
        return 1000
    end
    if xld != d.XL || yld != d.YL
        println(" DOMAIN DATA ERROR")
        return 1000
    end
    d.DT = dtseed; d.RE = re_seed
    uv = _take_values!(uvals, pos, length(d.U))
    copyto!(vec(parent(d.U)), uv)
    length(pvals) >= length(d.P) || error("seed_p.dat のデータ数が不足しています")
    copyto!(vec(parent(d.P)), view(pvals, 1:length(d.P)))
    opos = Ref(1)
    copyto!(vec(parent(d.CNX_OLD)), _take_values!(ovals, opos, length(d.CNX_OLD)))
    copyto!(vec(parent(d.CNY_OLD)), _take_values!(ovals, opos, length(d.CNY_OLD)))
    println("READ VELOCITY-SEED DATA OK!")
    return 0
end

function SEEDOUT!(d::CFDData; directory::AbstractString=".")
    open(joinpath(directory, "seed_u.dat"), "w") do io
        @printf(io, "%d %24.16E\n", d.NUM, d.TIME)
        @printf(io, "%d %d %24.16E %24.16E %24.16E %24.16E\n", IG, JG, d.XL, d.YL, d.DT, d.RE)
        for v in parent(d.U); @printf(io, " %.16E\n", v); end
    end
    open(joinpath(directory, "seed_p.dat"), "w") do io
        for v in parent(d.P); @printf(io, " %.16E\n", v); end
    end
    open(joinpath(directory, "seed_old.dat"), "w") do io
        for v in parent(d.CNX_OLD); @printf(io, " %.16E\n", v); end
        for v in parent(d.CNY_OLD); @printf(io, " %.16E\n", v); end
    end
    return nothing
end

function OUTPUT!(d::CFDData, wingname::AbstractString; outdir::AbstractString="TEXTDATA")
    # OUTPUT also writes the final instantaneous velocity field at NUM divisible by 5.
    if mod(d.NUM, 5) == 0
        _write_velocity_field(d, outdir)
    end
    open(joinpath(outdir, "output.dat"), "w") do io
        header(s) = rpad(first(s, min(length(s),71)), 71)
        @printf(io, "%s\n", header("#======================================================"))
        @printf(io, "%s\n", header("# SIMULATION OF FLOW PAST AN 2D AIRFOIL : "*wingname))
        @printf(io, "%s\n", header("#======================================================"))
        @printf(io, "%s\n", header("# CALC.NUM.     CALC.TIME"))
        @printf(io, "%10d%18.10E\n", d.NUM, d.TIME)
        @printf(io, "%s\n", header("# NX  NY  TIMESTEP  REYNOLDS"))
        @printf(io, "%4d%4d%11.6f%9.1f\n", IG, JG, d.DT, d.RE)
        @printf(io, "%s\n", header("# X-LENGTH  Y-LENGTH AOA(DEG)"))
        @printf(io, "%10.6f%10.6f%7.2f\n", d.XL, d.YL, d.DEG)
        @printf(io, "%s\n", header("# STATISTICS IN VELOCITY FIELD"))
        @printf(io, "%s\n", header("# I X DX XX DXX"))
        for i in 0:IG+1
            @printf(io, "%3d%14.6E%14.6E%14.6E%14.6E\n", i, d.X[i], d.DX[i], d.XX[i], d.DXX[i])
        end
        @printf(io, "%s\n", header("# J Y DY YY DYY"))
        for j in 0:JG+1
            @printf(io, "%3d%14.6E%14.6E%14.6E%14.6E\n", j, d.Y[j], d.DY[j], d.YY[j], d.DYY[j])
        end
        @printf(io, "%s\n", header("# X DX P(UPPER-SURFACE)  P(BOTTOM SURFACE)"))
        psum = 0.0
        for i in 1:IG
            ptop = 0.0; pbtm = 0.0
            for j in 2:JG
                if d.BETA[i,j,2]-d.BETA[i,j-1,2] >= 0.50; ptop = d.P[i,j]; end
                if d.BETA[i,j,2]-d.BETA[i,j-1,2] <= -0.50; pbtm = d.P[i,j-1]; end
            end
            if ptop != 0.0
                @printf(io, "%14.6E%14.6E%14.6E%14.6E\n", d.X[i], d.DX[i], ptop, pbtm)
            end
            psum += d.DX[i]*(pbtm-ptop)
        end
        cl = psum*2.0
        @printf(io, "%s\n", header("# Y DY P(DOWNSTREAM-SURFACE)  P(UPSTREAM SURFACE)"))
        psum = 0.0
        for j in 1:JG
            ptop = 0.0; pbtm = 0.0
            for i in 2:IG
                if d.BETA[i,j,1]-d.BETA[i-1,j,1] >= 0.5; ptop = d.P[i+1,j]; end
                if d.BETA[i,j,1]-d.BETA[i-1,j,1] <= -0.5; pbtm = d.P[i-2,j]; end
            end
            if ptop != 0.0
                @printf(io, "%14.6E%14.6E%14.6E%14.6E\n", d.Y[j], d.DY[j], ptop, pbtm)
            end
            psum += d.DY[j]*(pbtm-ptop)
        end
        cd = psum*2.0
        println("C_L =$(round(cl,digits=3)), C_D =$(round(cd,digits=3))")
    end
    return nothing
end
