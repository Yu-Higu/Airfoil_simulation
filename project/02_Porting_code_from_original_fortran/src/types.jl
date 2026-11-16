# Fortran incinit.f / incdata.f / incmesh.f の共通データ領域を移植
# 配列の添字範囲をFortran版に合わせるため OffsetArrays を使用する。

using OffsetArrays

const IG  = 384
const JG  = 144
const LX  = IG + 1
const LY  = JG + 1
const IGS = 52
const PI  = 4.0 * atan(1.0)

mutable struct CFDData
    # COMMON /DAT1/, /DAT2/, /STEP/
    RE::Float64
    DT::Float64
    TIME::Float64
    XL::Float64
    YL::Float64
    DEG::Float64
    NDRAW::Int
    NFIN::Int
    NUM::Int

    # COMMON /VELO/, /TEMP/, /OLDC/, /SCAL/, /PHIP/, /BIBM/
    U::OffsetArray{Float64,3}
    UT::OffsetArray{Float64,3}
    CNX_OLD::OffsetArray{Float64,2}
    CNY_OLD::OffsetArray{Float64,2}
    P::OffsetArray{Float64,2}
    PHI::OffsetArray{Float64,2}
    BETA::OffsetArray{Float64,3}

    # COMMON /MESH/
    X::OffsetArray{Float64,1}
    XX::OffsetArray{Float64,1}
    DX::OffsetArray{Float64,1}
    DXX::OffsetArray{Float64,1}
    Y::OffsetArray{Float64,1}
    YY::OffsetArray{Float64,1}
    DY::OffsetArray{Float64,1}
    DYY::OffsetArray{Float64,1}

    # COMMON /COE2/
    PDX::OffsetArray{Float64,1}
    PDXX::OffsetArray{Float64,1}
    DX2::OffsetArray{Float64,1}
    DXX2::OffsetArray{Float64,1}
    PDX2::OffsetArray{Float64,1}
    PDX3::OffsetArray{Float64,1}
    PDXX3::OffsetArray{Float64,1}
    PDY::OffsetArray{Float64,1}
    PDYY::OffsetArray{Float64,1}
    DY2::OffsetArray{Float64,1}
    DYY2::OffsetArray{Float64,1}
    PDY2::OffsetArray{Float64,1}
    PDY3::OffsetArray{Float64,1}
    PDYY3::OffsetArray{Float64,1}

    # COMMON /CONS/, /IDXY/, /GAMM/
    D_N::OffsetArray{Float64,1}
    D_E::OffsetArray{Float64,1}
    D_W::OffsetArray{Float64,1}
    D_S::OffsetArray{Float64,1}
    PDPHI::OffsetArray{Float64,2}
    IDX::OffsetArray{Int,2}
    IDY::OffsetArray{Int,2}
    GAMMA::OffsetArray{Float64,3}
end

# Fortran の DIMENSION 宣言と同じ添字範囲で確保する。
z1(lo, hi) = OffsetArray(zeros(Float64, hi - lo + 1), lo:hi)
z2(i1, i2, j1, j2) =
    OffsetArray(zeros(Float64, i2-i1+1, j2-j1+1), i1:i2, j1:j2)
z3(i1, i2, j1, j2, k1, k2) =
    OffsetArray(zeros(Float64, i2-i1+1, j2-j1+1, k2-k1+1),
               i1:i2, j1:j2, k1:k2)
i2(i1, i2, j1, j2) =
    OffsetArray(zeros(Int, i2-i1+1, j2-j1+1), i1:i2, j1:j2)

function CFDData()
    return CFDData(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0,
        z3(0, IG+1, 0, JG+1, 1, 2),
        z3(0, IG+1, 0, JG+1, 1, 2),
        z2(0, IG, 0, JG), z2(0, IG, 0, JG),
        z2(0, IG+1, 0, JG+1), z2(0, IG+1, 0, JG+1),
        z3(0, IG+1, 0, JG+1, 1, 2),
        z1(-3, LX+5), z1(-3, LX+5),
        z1(-3, LX+5), z1(-3, LX+5),
        z1(-3, LY+5), z1(-3, LY+5),
        z1(-3, LY+5), z1(-3, LY+5),
        z1(0, LX+2), z1(0, LX+2),
        z1(0, LX+2), z1(0, LX+2),
        z1(0, LX+2), z1(0, LX+2), z1(0, LX+2),
        z1(0, LY+2), z1(0, LY+2),
        z1(0, LY+2), z1(0, LY+2),
        z1(0, LY+2), z1(0, LY+2), z1(0, LY+2),
        z1(0, JG), z1(0, IG), z1(0, IG), z1(0, JG),
        z2(0, IG, 0, JG),
        i2(0, IG, 0, JG), i2(0, IG, 0, JG),
        z3(0, IG, 0, JG, 1, 2)
    )
end
