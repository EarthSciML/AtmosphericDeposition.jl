using DepositionMTK
using Test
using StaticArrays

include("../src/DepositionMTK.jl")

# Results from Wesely (1989) table 3; updated to values from Walmsley (1996) table 1.

begin
	# r_i represents the minimum bulk canopy resistances for water vapor.
	const r_i = SA_F32[
	inf 60 120 70 130 100 inf inf 80 100 150
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 400 800 inf inf inf inf inf
	inf 120 240 140 250 190 inf inf 160 200 300]

	# r_lu signifies leaf cuticles in healthy vegetation and otherwise the outer surfaces in the upper canopy.
	const r_lu = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf inf inf 6000 9000 inf inf 9000 9000 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

	# r_ac signifies transfer that depends only on canopy height and density.
	const r_ac = SA_F32[
	100.0 200 100 2000 2000 2000 0 0 300 150 200
	100 150 100 1500 2000 1700 0 0 200 120 140
	100 10 100 1000 2000 1500 0 0 100 50 120
	100 10 10 1000 2000 1500 0 0 50 10 50
	100 50 80 1200 2000 1500 0 0 200 60 120]

	# r_gs signifies uptake at the "ground" by soil, leaf litter, snow, water etc. 'S' and 'O' stand for SO2 and O3 respectively.
	const r_gsS = SA_F32[
	400.0 150 350 500 500 100 0 1000 0 220 40
	400 200 350 500 500 100 0 1000 0 300 400
	400 150 350 500 500 200 0 1000 0 200 400
	100 100 100 100 100 100 0 1000 100 100 50
	500 150 350 500 500 200 0 1000 0 250 40]

	const r_gsO = SA_F32[
	300.0 150 200 200 200 300 2000 400 1000 180 200
	300 150 200 200 200 300 2000 400 800 180 200
	300 150 200 200 200 300 2000 400 1000 180 20
	600 3500 3500 3500 3500 3500 2000 400 3500 3500 3500
	300 150 200 200 200 300 2000 400 1000 180 200]

	# r_cl is meant to account for uptake pathways at the leaves, bark, etc. 'S' and 'O' stand for SO2 and O3 respectively.
	const r_clS = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 2000 4000 inf inf 9000 9000 9000
	inf inf 9000 9000 3000 6000 inf inf 9000 9000 9000
	inf inf inf 9000 200 400 inf inf 9000 inf 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

	const r_clO = SA_F32[
	inf 1000 1000 1000 1000 1000 inf inf 1000 1000 1000
	inf 400 400 400 1000 600 inf inf 400 400 400
	inf 1000 400 400 1000 600 inf inf 800 600 600
	inf 1000 1000 400 1500 600 inf inf 800 1000 800
	inf 1000 500 500 1500 700 inf inf 600 800 800]

	# Holder for gas properties from Wesely (1989) Table 2.'
	struct GasData
		Dh2oPerDx::AbstractFloat
		Hstar::AbstractFloat
		Fo::AbstractFloat
	end

	const So2Data = GasData(1.9, 1.e5, 0)
	const O3Data  = GasData(1.6, 0.01, 1)
	const No2Data = GasData(1.6, 0.01, 0.1) # Wesely (1989) suggests that,
	# in general, the sum of NO and NO2 should be considered rather
	# than NO2 alone because rapid in-air chemical reactions can cause
	# a significant change of NO and NO2 vertical fluxes between the
	# surface and the point at which deposition velocities are applied,
	# but the sum of NO and NO2 fluxes should be practically unchanged.
	const NoData   = GasData(1.3, 3.e-3, 0) # Changed according to Walmsley (1996)
	const Hno3Data = GasData(1.9, 1.e14, 0)
	const H2o2Data = GasData(1.4, 1.e5, 1)
	const AldData  = GasData(1.6, 15, 0)     # Acetaldehyde (aldehyde class)
	const HchoData = GasData(1.3, 6.e3, 0)   # Formaldehyde
	const OpData   = GasData(1.6, 240, 0.1)  
	# Methyl hydroperoxide (organic peroxide class)
	const PaaData  = GasData(2.0, 540, 0.1)  # Peroxyacetyl nitrate
	const OraData  = GasData(1.6, 4.e6, 0)   # Formic acid (organic acid class)
	const Nh3Data  = GasData(0.97, 2.e4, 0)  # Changed according to Walmsley (1996)
	const PanData  = GasData(2.6, 3.6, 0.1)  # Peroxyacetyl nitrate
	const Hno2Data = GasData(1.6, 1.e5, 0.1) # Nitrous acid
end

begin
	const SO2 = SA_F32[
		130.0 140 160 380 1000 100 1200
		1400 1400 1400 1400 1500 100 1300
		1100 1100 1100 1100 1200 90 1000
		1000 1000 1000 1000 1100 1100 1100
		270 290 330 620 1100 90 1000]

	const O3 = SA_F32[
		100.0 110 130 320 960 960 580
		430 470 520 710 1300 950 580
		390 420 460 610 960 770 510
		560 620 710 1100 3200 3200 3200
		180 200 230 440 950 820 530]

	const NO2 = SA_F32[
		120.0 130 160 480 2900 2700 2300
		1900 1900 1900 2000 2700 2500 2200
		1700 1700 1800 1900 2400 2300 2000
		3900 4000 4100 4500 9999 9999 9999
		270 290 350 850 2500 2300 2000]

	const H2O2 = SA_F32[
		90.0 90 110 250 640 90 80
		400 430 480 650 1100 90 90
		370 390 430 550 840 90 80
		400 430 470 620 1000 1000 1000
		160 170 200 370 750 90 80]

	const ALD = SA_F32[
		330.0 340 370 800 9999 9999 9999
		9999 9999 9999 9999 9999 9999 9999
		9999 9999 9999 9999 9999 9999 9999
		9999 9999 9999 9999 9999 9999 9999
		520 550 630 1700 9999 9999 9999]

	const HCHO = SA_F32[
		100.0 110 140 450 6700 1400 1400
		8700 8700 8700 8700 8700 1400 1400
		8300 8300 8300 8300 8400 1400 1400
		2900 2900 2900 2900 2900 2900 2900
		250 270 340 1000 7500 1400 1400]

	const OP = SA_F32[
		120.0 130 160 480 2800 2500 2200
		1900 1900 1900 2000 2700 2400 2000
		1700 1700 1800 1800 2400 2100 1900
		3700 3700 3800 4200 8600 8600 8600
		270 290 350 850 2500 2200 1900]

	const PAA = SA_F32[
		150.0 160 200 580 2800 2400 2000
		1900 1900 1900 2000 2700 2200 1900
		1700 1700 1700 1800 2400 2000 1800
		3400 3400 3500 3800 7200 7200 7200
		330 350 420 960 2400 2100 1800]

	const ORA = SA_F32[
		30.0 30 30 40 50 10 10
		140 140 150 170 190 10 10
		130 140 140 160 180 10 10
		310 340 390 550 910 910 910
		60 60 70 80 90 10 10]

	const NH3 = SA_F32[
		80.0 80 100 320 2700 430 430
		3400 3400 3400 3400 3400 440 440
		3000 3000 3000 3000 3100 430 430
		1500 1500 1500 1500 1500 1500 1500
		180 200 240 680 2800 430 430]
	const PAN = SA_F32[
		190.0 210 250 700 2900 2700 2300
		1900 1900 1900 2000 2700 2500 2200
		1700 1700 1800 1900 2400 2300 2000
		3900 4000 4100 4500 9999 9999 9999
		410 430 510 1100 2500 2300 2000]
	const HNO2 = SA_F32[
		110.0 120 140 330 950 90 90
		1000 1000 1000 1100 1400 90 90
		860 860 870 910 1100 90 90
		820 830 830 870 1000 1000 1000
		220 240 280 530 1000 90 90]
end

function different(a::Float64, b::Float32) where T<:AbstractFloat
	c = abs(a - b)
	return c/b > 0.1 && c >= 11.0
end

function TestWesely()
	iLandUse = 4                       			# deciduous forest
	Ts = [25.0, 10, 2, 0, 10]        			# Surface Temperature [C]
	Garr = [800.0, 500, 300, 100, 0] 			# Solar radiation [W m-2]
	θ = 0.0                                  	# Slope [radians]

	polNames = [
		"SO2", "O3", "NO2", "H2O2", 
		"ALD", "HCHO", "OP", "PAA", 
		"ORA", "NH3", "PAN", "HNO2"]
	
	testData = [SO2, O3, NO2, H2O2, ALD, HCHO, OP, PAA, ORA, NH3, PAN, HNO2]
	gasData = [
		So2Data, O3Data, No2Data, 
		H2o2Data, AldData, HchoData, 
		OpData, PaaData, OraData, 
		Nh3Data, PanData, Hno2Data]
	
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
				r_c = SurfaceResistance(gasData[i], G, Ts[iSeason], θ,
					iSeason, iLandUse, false, false, isSO2, isO3)
				if different(r_c, polData[iSeason, ig])
					# println(pol, iSeason, G, r_c, polData[iSeason, ig])
					return false
				end
			end
			r_c = SurfaceResistance(gasData[i], 0.0, Ts[iSeason], θ,
				iSeason, iLandUse, false, true, isSO2, isO3) # dew
			if different(r_c, polData[iSeason, 6])
				# println(pol, iSeason, "dew", r_c, polData[iSeason, 6])
				return false;
			end
			r_c = SurfaceResistance(gasData[i], 0.0, Ts[iSeason], θ,
				iSeason, iLandUse, true, false, isSO2, isO3) # rain
			if different(r_c, polData[iSeason, 7])
				# println(pol, iSeason, "rain", r_c, polData[iSeason, 7])
				return false;
			end
		end
	end
	return true
end

@testset "DepositionMTK.jl" begin
    @test TestWesely() == true
    @test r_s(1.0, 1.0, 1, 1, true) ≈ 3.0772306344553016e30
    @test r_dc(1.0, 1.0) ≈ 9.181727363545544
    @test r_mx(1.0, 1.0) ≈ 0.009999966666777778
    @test r_smx(1.0, 1.0, 1.0) ≈ 2.0
    @test r_lux(1.0, 1.0, 1, 1, true, false, true, false) ≈ 50
    @test r_clx(1.0, 1.0, 1, 1) ≈ 9.999899563027895e24
    @test r_gsx(1.0, 1.0, 1, 1) ≈ 299.9977500168749
    @test SurfaceResistance(So2Data, 1.0, 1.0, 1.0, 1, 1, true, true, true, false) ≈ 45.45454545454546
end
