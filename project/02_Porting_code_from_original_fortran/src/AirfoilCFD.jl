module AirfoilCFD

using OffsetArrays
using Printf

include("types.jl")
include("mesh.jl")
include("solver.jl")
include("boundary.jl")
include("io.jl")

export CFDData, CALBETA!, FS0!, FS1!, FS2!, UBC!, OUTPUT!, SEEDIN!, SEEDOUT!, run_simulation

function run_simulation(re::Real, deg::Real; nc::Int=50, ng::Int=5,
                        airfoil_file::AbstractString="data/clark_y.dat",
                        output_dir::AbstractString="TEXTDATA", write_seed::Bool=true)
    nc >= 1 || throw(ArgumentError("NC must be at least 1"))
    ng >= 1 || throw(ArgumentError("NG must be at least 1"))
    isfile(airfoil_file) || error("翼座標ファイルが見つかりません: $airfoil_file")
    mkpath(output_dir)

    d = CFDData()
    d.RE = Float64(re)
    d.DEG = Float64(deg)
    d.DT = 0.002500
    d.XL = 8.0
    d.YL = 3.0
    d.NUM = 0
    d.TIME = 0.0

    println("REYNOLDS NUMBER, RE = $(d.RE)")
    println("ANGLE OF ATTACK, Degree = $(d.DEG)")
    wingname = CALBETA!(d, airfoil_file)

    div, cfl, usum, cl, cd = UBC!(d; outdir=output_dir)
    nn = 0
    # FIRST STEP
    nn += 1; d.NUM += 1; d.TIME += d.DT
    FS0!(d)
    ncount = FS2!(d)
    div, cfl, usum, cl, cd = UBC!(d; outdir=output_dir)

    ts_path = joinpath(output_dir, "time-series.dat")
    open(ts_path, "w") do ts
        println(ts, "    NN     NUM    TIME    DIV         CFL         AVE.U     C_L      C_D     NCOUNT")
        @printf(ts, "%6d%8d%10.4f%12.4E%12.4E%12.4E%9.4f%9.4f%8d\n",
            nn, d.NUM, d.TIME, div, cfl, usum, cl, cd, ncount)
        @printf("%6d%8d%10.4f%12.4E%12.4E%12.4E%9.4f%9.4f%8d\n",
            nn, d.NUM, d.TIME, div, cfl, usum, cl, cd, ncount)

        # SECOND STEP onward; preserve the original NST=2:NC loop.
        for nst in 2:nc
            nn += 1; d.NUM += 1; d.TIME += d.DT
            FS1!(d)
            ncount = FS2!(d)
            div, cfl, usum, cl, cd = UBC!(d; outdir=output_dir)
            if mod(nn, ng) == 0
                @printf(ts, "%6d%8d%10.4f%12.4E%12.4E%12.4E%9.4f%9.4f%8d\n",
                    nn, d.NUM, d.TIME, div, cfl, usum, cl, cd, ncount)
                @printf("%6d%8d%10.4f%12.4E%12.4E%12.4E%9.4f%9.4f%8d\n",
                    nn, d.NUM, d.TIME, div, cfl, usum, cl, cd, ncount)
            end
        end
    end
    if write_seed
        SEEDOUT!(d; directory=".")
    end
    OUTPUT!(d, wingname; outdir=output_dir)
    return d
end

end # module
