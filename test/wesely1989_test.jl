@testsnippet WeselySetup begin
    using AtmosphericDeposition
    using StaticArrays

    """
    Results from Wesely (1989) table 3; updated to values from Walmsley (1996) table 1.
    """
    const SO2 = SA_F32[
        130.0 140 160 380 1000 100 1200
        1400 1400 1400 1400 1500 100 1300
        1100 1100 1100 1100 1200 90 1000
        1000 1000 1000 1000 1100 1100 1100
        270 290 330 620 1100 90 1000
    ]

    const O3 = SA_F32[
        100.0 110 130 320 960 960 580
        430 470 520 710 1300 950 580
        390 420 460 610 960 770 510
        560 620 710 1100 3200 3200 3200
        180 200 230 440 950 820 530
    ]

    const NO2 = SA_F32[
        120.0 130 160 480 2900 2700 2300
        1900 1900 1900 2000 2700 2500 2200
        1700 1700 1800 1900 2400 2300 2000
        3900 4000 4100 4500 9999 9999 9999
        270 290 350 850 2500 2300 2000
    ]

    const H2O2 = SA_F32[
        90.0 90 110 250 640 90 80
        400 430 480 650 1100 90 90
        370 390 430 550 840 90 80
        400 430 470 620 1000 1000 1000
        160 170 200 370 750 90 80
    ]

    const ALD = SA_F32[
        330.0 340 370 800 9999 9999 9999
        9999 9999 9999 9999 9999 9999 9999
        9999 9999 9999 9999 9999 9999 9999
        9999 9999 9999 9999 9999 9999 9999
        520 550 630 1700 9999 9999 9999
    ]

    const HCHO = SA_F32[
        100.0 110 140 450 6700 1400 1400
        8700 8700 8700 8700 8700 1400 1400
        8300 8300 8300 8300 8400 1400 1400
        2900 2900 2900 2900 2900 2900 2900
        250 270 340 1000 7500 1400 1400
    ]

    const OP = SA_F32[
        120.0 130 160 480 2800 2500 2200
        1900 1900 1900 2000 2700 2400 2000
        1700 1700 1800 1800 2400 2100 1900
        3700 3700 3800 4200 8600 8600 8600
        270 290 350 850 2500 2200 1900
    ]

    const PAA = SA_F32[
        150.0 160 200 580 2800 2400 2000
        1900 1900 1900 2000 2700 2200 1900
        1700 1700 1700 1800 2400 2000 1800
        3400 3400 3500 3800 7200 7200 7200
        330 350 420 960 2400 2100 1800
    ]

    const ORA = SA_F32[
        30.0 30 30 40 50 10 10
        140 140 150 170 190 10 10
        130 140 140 160 180 10 10
        310 340 390 550 910 910 910
        60 60 70 80 90 10 10
    ]

    const NH3 = SA_F32[
        80.0 80 100 320 2700 430 430
        3400 3400 3400 3400 3400 440 440
        3000 3000 3000 3000 3100 430 430
        1500 1500 1500 1500 1500 1500 1500
        180 200 240 680 2800 430 430
    ]
    const PAN = SA_F32[
        190.0 210 250 700 2900 2700 2300
        1900 1900 1900 2000 2700 2500 2200
        1700 1700 1800 1900 2400 2300 2000
        3900 4000 4100 4500 9999 9999 9999
        410 430 510 1100 2500 2300 2000
    ]
    const HNO2 = SA_F32[
        110.0 120 140 330 950 90 90
        1000 1000 1000 1100 1400 90 90
        860 860 870 910 1100 90 90
        820 830 830 870 1000 1000 1000
        220 240 280 530 1000 90 90
    ]

    function different(a::Float64, b::Float32) #where T<:AbstractFloat
        c = abs(a - b)
        return c / b > 0.1 && c >= 11.0
    end

    function TestWesely()
        iLandUse = 4                       # deciduous forest
        Ts = [25.0, 10, 2, 0, 10]        # Surface Temperature [C]
        Garr = [800.0, 500, 300, 100, 0] # Solar radiation [W m-2]
        θ = 0.0                                  # Slope [radians]

        polNames = [
            "SO2",
            "O3",
            "NO2",
            "H2O2",
            "ALD",
            "HCHO",
            "OP",
            "PAA",
            "ORA",
            "NH3",
            "PAN",
            "HNO2",
        ]

        testData = [SO2, O3, NO2, H2O2, ALD, HCHO, OP, PAA, ORA, NH3, PAN, HNO2]
        gasData = [
            AtmosphericDeposition.So2Data,
            AtmosphericDeposition.O3Data,
            AtmosphericDeposition.No2Data,
            AtmosphericDeposition.H2o2Data,
            AtmosphericDeposition.AldData,
            AtmosphericDeposition.HchoData,
            AtmosphericDeposition.OpData,
            AtmosphericDeposition.PaaData,
            AtmosphericDeposition.OraData,
            AtmosphericDeposition.Nh3Data,
            AtmosphericDeposition.PanData,
            AtmosphericDeposition.Hno2Data,
        ]

        for i in 1:12
            pol = polNames[i]
            polData = testData[i]
            isSO2, isO3 = false, false
            if pol == "SO2"
                isSO2 = true
            end
            if pol == "O3"
                isO3 = true
            end
            for iSeason in 1:5
                for ig in 1:5
                    G = Garr[ig]
                    r_c = WesleySurfaceResistance(
                        gasData[i],
                        G,
                        Ts[iSeason],
                        θ,
                        iSeason,
                        iLandUse,
                        false,
                        false,
                        isSO2,
                        isO3
                    )
                    if different(r_c, polData[iSeason, ig])
                        println(pol, iSeason, G, r_c, polData[iSeason, ig])
                        return false
                    end
                end
                r_c = WesleySurfaceResistance(
                    gasData[i],
                    0.0,
                    Ts[iSeason],
                    θ,
                    iSeason,
                    iLandUse,
                    false,
                    true,
                    isSO2,
                    isO3
                ) # dew
                if different(r_c, polData[iSeason, 6])
                    println(pol, iSeason, "dew", r_c, polData[iSeason, 6])
                    return false
                end
                r_c = WesleySurfaceResistance(
                    gasData[i],
                    0.0,
                    Ts[iSeason],
                    θ,
                    iSeason,
                    iLandUse,
                    true,
                    false,
                    isSO2,
                    isO3
                ) # rain
                if different(r_c, polData[iSeason, 7])
                    println(pol, iSeason, "rain", r_c, polData[iSeason, 7])
                    return false
                end
            end
        end
        return true
    end
end

@testitem "wesley1989.jl" setup = [WeselySetup] begin
    @test TestWesely() == true
    @test AtmosphericDeposition.r_s(1.0, 1.0, 1, 1, true) ≈ 3.0772306344553016e30
    @test AtmosphericDeposition.r_dc(1.0, 1.0) ≈ 9.181727363545544
    @test AtmosphericDeposition.r_mx(1.0, 1.0) ≈ 0.009999966666777778
    @test AtmosphericDeposition.r_smx(1.0, 1.0, 1.0) ≈ 2.0
    @test AtmosphericDeposition.r_lux(1.0, 1.0, 1, 1, true, false, true, false) ≈ 50
    @test AtmosphericDeposition.r_clx(1.0, 1.0, 1, 1) ≈ 9.999899563027895e24
    @test AtmosphericDeposition.r_gsx(1.0, 1.0, 1, 1) ≈ 299.9977500168749
    @test WesleySurfaceResistance(
        AtmosphericDeposition.So2Data,
        1.0,
        1.0,
        1.0,
        1,
        1,
        true,
        true,
        true,
        false
    ) ≈ 45.45454545454546
end

@testitem "r_s closes stomata outside 0.1-39.9 degC" setup = [WeselySetup] begin
    # Wesely (1989) eq. 3 is only defined for 0 < Ts < 40 degC.
    for Ts in (-5.0, 0.0, 40.5, 45.0, 55.0)
        @test AtmosphericDeposition.r_s(300.0, Ts, 1, 6, false) == AtmosphericDeposition.inf
    end
    @test AtmosphericDeposition.r_s(300.0, 20.0, 1, 6, false) < 1.0e3

    # Rc must stay in a physical range rather than hitting either clamp.
    for Ts in (45.0, 55.0)
        rc = WesleySurfaceResistance(
            AtmosphericDeposition.O3Data, 300.0, Ts, 0.0, 1, 6, false, false, false, true)
        @test 100.0 < rc < 9999.0
    end
end

@testitem "obtain_value returns exact table entries" setup = [WeselySetup] begin
    # Regression guard. `obtain_value` must be an exact categorical lookup:
    # (iSeason, iLandUse) are class indices, not continuous coordinates, so
    # every query must reproduce the table entry bit-for-bit.
    #
    # This has regressed twice. From 2024-09-06 (46379366) until 2026-05-22 the
    # lookup flattened the 5×11 table and ran `LinearInterpolation` over the flat
    # index. That is exact at interior nodes but wrong at the LAST one: for
    # `r_i[5,11]` the flat index is 55, so the closing segment
    # `u[54] + (u[55]-u[54])` is evaluated, and `u[54]` is `r_i[4,11] = inf =
    # 1.0e25`. Since `ulp(1e25) ≈ 2.1e9 ≫ 300`, the 300 is annihilated and the
    # result was `0.0` instead of `300.0` (a magnitude problem, not a Float32
    # one — it happens in Float64 too). It was fixed on 2026-05-22 and reverted
    # a day later by 9f5cebf0.
    #
    # Consequence if it regresses: `r_i = 0` → `r_s = 0` → `rsmx = r_mx`, and in
    # the parallel-conductance sum `1/rsmx` swamps every other pathway, so `Rc`
    # clamps to its 10 s/m floor and Vd is maximal over rocky shrubs in
    # Transitional season (Mar–May). Measured Vd overestimate: O3 4.8×, NO2 10.5×.
    tables = [
        "r_i" => AtmosphericDeposition.r_i,
        "r_lu" => AtmosphericDeposition.r_lu,
        "r_ac" => AtmosphericDeposition.r_ac,
        "r_gsS" => AtmosphericDeposition.r_gsS,
        "r_gsO" => AtmosphericDeposition.r_gsO,
        "r_clS" => AtmosphericDeposition.r_clS,
        "r_clO" => AtmosphericDeposition.r_clO,
    ]
    for (name, m) in tables, iSeason in 1:5, iLandUse in 1:11
        got = AtmosphericDeposition.obtain_value(iSeason, iLandUse, m)
        want = Float64(m[iSeason, iLandUse])
        @test got ≈ want rtol = 1.0e-9
    end

    # The specific entry that regressed, called out so a failure is unmistakable.
    @test AtmosphericDeposition.obtain_value(5, 11, AtmosphericDeposition.r_i) ≈ 300.0

    # ...and its physical consequence: Rc must NOT sit on the 10 s/m floor here.
    rc_o3 = WesleySurfaceResistance(
        AtmosphericDeposition.O3Data, 800.0, 20.0, 0.0, 5, 11,
        false, false, false, true
    )
    @test rc_o3 > 100.0
end

@testitem "zero r_ac and r_gsS entries do not zero the ground resistance" setup = [WeselySetup] begin
    # Water has r_ac = r_gsS = 0 in every season. Insoluble gases must not reach the
    # 10 s/m floor there, while highly soluble gases still do.
    Rc(g, iLandUse; isSO2 = false, isO3 = false) = WesleySurfaceResistance(
        g, 300.0, 25.0, 0.0, 1, iLandUse, false, false, isSO2, isO3)
    water, wetland = 7, 9
    @test Rc(AtmosphericDeposition.No2Data, water) == 9999.0
    @test Rc(AtmosphericDeposition.NoData, water) == 9999.0
    @test Rc(AtmosphericDeposition.PANData, water) ≈ 1866.67 rtol = 1.0e-4
    @test Rc(AtmosphericDeposition.No2Data, wetland) ≈ 188.27 rtol = 1.0e-4
    @test Rc(AtmosphericDeposition.Hno3Data, water) == 10.0
    @test Rc(AtmosphericDeposition.H2o2Data, water) == 10.0
    @test Rc(AtmosphericDeposition.So2Data, water; isSO2 = true) == 10.0
    @test Rc(AtmosphericDeposition.O3Data, water; isO3 = true) ≈ 2001.0
end
