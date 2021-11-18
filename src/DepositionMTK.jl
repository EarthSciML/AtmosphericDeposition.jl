module DepositionMTK

# Write your package code here.

### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 91119eae-3023-4097-bd50-8079e3ae69c7
using StaticArrays

# ╔═╡ 7363fcd7-b115-4c06-ab0f-429daaf56a83
md"""
# Data file
"""

# ╔═╡ 96a4a0c0-321b-11ec-09d7-3bce044720f0
const inf = 1.e25

# ╔═╡ e655d9ac-8049-4950-831d-4ef4c421cdbb
md"""
r_i represents the minimum bulk canopy resistances for water vapor.
"""

# ╔═╡ ef83b8f7-53fa-46b7-86de-2300532004b1
const r_i = SA_F32[
	inf 60 120 70 130 100 inf inf 80 100 150
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 400 800 inf inf inf inf inf
	inf 120 240 140 250 190 inf inf 160 200 300]

# ╔═╡ 9d86d682-92bd-4b26-ae23-c4a385374a7c
md"""
r_lu signifies leaf cuticles in healthy vegetation and otherwise the outer surfaces in the upper canopy.
"""

# ╔═╡ 301eaff9-06c8-4726-9067-27a5e813f86c
const r_lu = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf inf inf 6000 9000 inf inf 9000 9000 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

# ╔═╡ 9a7cb29b-85ba-4722-8bd4-329684837b83
md"""
r_ac signifies transfer that depends only on canopy height and density.
"""

# ╔═╡ 2b446731-9bbe-4f64-8c2c-e4f1dd7e5f41
const r_ac = SA_F32[
	100.0 200 100 2000 2000 2000 0 0 300 150 200
	100 150 100 1500 2000 1700 0 0 200 120 140
	100 10 100 1000 2000 1500 0 0 100 50 120
	100 10 10 1000 2000 1500 0 0 50 10 50
	100 50 80 1200 2000 1500 0 0 200 60 120]

# ╔═╡ ad0cd7a5-5e67-4f0d-b1cb-9d9661e3548e
md"""
r_gs signifies uptake at the "ground" by soil, leaf litter, snow, water etc. 'S' and 'O' stand for SO2 and O3 respectively.
"""

# ╔═╡ c36b7c6e-648a-4bf8-8797-d4826c83ce2a
const r_gsS = SA_F32[
	400.0 150 350 500 500 100 0 1000 0 220 40
	400 200 350 500 500 100 0 1000 0 300 400
	400 150 350 500 500 200 0 1000 0 200 400
	100 100 100 100 100 100 0 1000 100 100 50
	500 150 350 500 500 200 0 1000 0 250 40]

# ╔═╡ 1d2853eb-3a6d-4313-b6a1-4858a8708c1d
const r_gsO = SA_F32[
	300.0 150 200 200 200 300 2000 400 1000 180 200
	300 150 200 200 200 300 2000 400 800 180 200
	300 150 200 200 200 300 2000 400 1000 180 20
	600 3500 3500 3500 3500 3500 2000 400 3500 3500 3500
	300 150 200 200 200 300 2000 400 1000 180 200]

# ╔═╡ 15daae3e-8acd-4ee2-9b8f-f36a3e111c26
md"""
r_cl is meant to account for uptake pathways at the leaves, bark, etc. 'S' and 'O' stand for SO2 and O3 respectively.
"""

# ╔═╡ e791cac6-4890-42d9-b770-99b9ff8d97c7
const r_clS = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 2000 4000 inf inf 9000 9000 9000
	inf inf 9000 9000 3000 6000 inf inf 9000 9000 9000
	inf inf inf 9000 200 400 inf inf 9000 inf 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

# ╔═╡ 8a78d619-8632-4fdd-800b-74f1a0b8a9d7
const r_clO = SA_F32[
	inf 1000 1000 1000 1000 1000 inf inf 1000 1000 1000
	inf 400 400 400 1000 600 inf inf 400 400 400
	inf 1000 400 400 1000 600 inf inf 800 600 600
	inf 1000 1000 400 1500 600 inf inf 800 1000 800
	inf 1000 500 500 1500 700 inf inf 600 800 800]

# ╔═╡ 431c0a99-9a71-4284-bc6a-2f95a879b3e1
md"""
Holder for gas properties from Wesely (1989) Table 2.
"""

# ╔═╡ f160601e-e9e5-4181-8c4a-425c35838374
struct GasData
	Dh2oPerDx::AbstractFloat
	Hstar::AbstractFloat
	Fo::AbstractFloat
end

# ╔═╡ adfb7cbf-3dc0-4fa7-9adc-66e203accffc
md"""
Properties of various gases from Wesely (1989) Table 2.
"""

# ╔═╡ 34f744cb-d9fc-44b4-a61e-3262ee2b1e65
begin
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

# ╔═╡ 31510e91-8f94-44ee-8163-2905249ba2fd
md"""
# Surface resistance file
"""

# ╔═╡ 5b0071b0-aab8-4800-a9ed-9dfe8947cb90
begin
	const Midsummer =  1 		# 0: Midsummer with lush vegetation
	const Autumn = 2        	# 1: Autumn with unharvested cropland
	const LateAutumn = 3    	# 2: Late autumn after frost, no snow
	const Winter = 4        	# 3: Winter, snow on ground and subfreezing
	const Transitional = 5  	# 4: Transitional spring with partially green short 								     annuals
end

# ╔═╡ fa59f3f5-1115-4329-9b1d-a6f943ec70aa
begin
	const Urban = 1         	# 0: Urban land
	const Agricultural = 2      # 1: Agricultural land
	const Range = 3             # 2: Range land
	const Deciduous = 4         # 3: Deciduous forest
	const Coniferous = 5        # 4: Coniferous forest
	const MixedForest = 6       # 5: Mixed forest including wetland
	const Water = 7             # 6: Water, both salt and fresh
	const Barren = 8            # 7: Barren land, mostly desert
	const Wetland = 9           # 8: Nonforested wetland
	const RangeAg = 10          # 9: Mixed agricultural and range land
	const RockyShrubs = 11      # 10: Rocky open areas with low-growing shrubs
end

# ╔═╡ 344d4e1a-8c6e-4050-bdea-945b60d2aa6e
function max(a::T, b::T) where T<:AbstractFloat
	if a > b
		return a
	else
		return b
	end
end

# ╔═╡ 4f47597b-b383-4ebd-9d1c-662c23ed8396
function min(a::T, b::T) where T<:AbstractFloat
	if a < b
		return a
	else
		return b
	end
end

# ╔═╡ c7c1aade-4eba-46e3-9d4a-e89a15a8d21b
function r_s(G::T, Ts::T, iSeason::Int, iLandUse::Int, rainOrDew::Bool) where T<:AbstractFloat
	rs = 0.0
	if Ts >= 39.9 || Ts <= 0.1
		rs = inf
	else
		rs = r_i[iSeason, iLandUse] * (1 + (200.0 * 1.0 / (G + 1.0))^2) *
			(400.0 * 1.0 / (Ts * (40.0 - Ts)))
	end
	
	if rainOrDew
		rs *= 3
	end
	return rs
end

# ╔═╡ f0481175-a676-4df0-ad1e-aabcfefb0995
function r_dc(G::T, θ::T) where T<:AbstractFloat
	return 100.0 * (1.0 + 1000.0/(G+10.0)) / (1.0 + 1000.0*θ)
end

# ╔═╡ 09352c4a-0a0c-4640-b36d-c69b805ee24b
function r_mx(Hstar::T, fo::T) where T<:AbstractFloat
	return 1.0 / (Hstar/3000.0 + 100.0*fo)
end

# ╔═╡ 0d24bfd4-bc76-4f0a-8814-a3e7cf1c9db1
function r_smx(r_s::T, Dh2oPerDx::T, r_mx::T) where T<:AbstractFloat
	return r_s * Dh2oPerDx + r_mx
end

# ╔═╡ 2e36c0e5-8037-4421-8370-2fbb49bdfe4a
1.e-3*3

# ╔═╡ b2af8534-556c-4278-81f3-e6d9e81604ec
function r_lux(Hstar::T, fo::T, iSeason::Int, iLandUse::Int, rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool) where T<:AbstractFloat
	rlux = 0.0
	if dew && (iSeason != 4)
		if isSO2
			if iLandUse == 1
				rlux = 50.0
			else
				rlux = 100
			end
		elseif isO3
			rlux = 1.0 / (1.0/3000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
		else
			rluO = 1.0 / (1.0/3000.0 + 1.0/(3*r_lu[iSeason, iLandUse])) 
			rlux = 1.0 / (1.0/(3*r_lu[iSeason, iLandUse]/(1.0e-5*Hstar+fo)) 
				+ 1.0e-7*Hstar + fo/rluO)
		end
	elseif rain && (iSeason != 4)
		if isSO2
			if iLandUse == 1
				rlux = 50
			else
				rlux = 1.0 / (1.0/5000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
			end
		elseif isO3
			rlux = 1.0 / (1.0/1000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
		else 
			rluO = 1.0 / (1.0/1000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
			rlux = 1.0 / (1.0/(3*r_lu[iSeason, iLandUse]/(1.e-5*Hstar+fo)) + 1.0e-7*Hstar +fo/rluO)
		end				
	else
		rlux = r_lu[iSeason, iLandUse] / (1.0e-5*Hstar + fo)
	end
	return rlux
end

# ╔═╡ fbfe0521-8abe-43e3-8453-4f2453dd635e
function r_clx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
	return 1.0 / (Hstar/(1.0e5*r_clS[iSeason, iLandUse]) +
		fo/r_clO[iSeason, iLandUse])
end

# ╔═╡ 31d79b1d-dc72-467d-9af0-4a7482124f47
function r_gsx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
	return 1.0 / (Hstar/(1.0e5*r_gsS[iSeason, iLandUse]) +
		fo/r_gsO[iSeason, iLandUse])
end

# ╔═╡ 907b6856-0702-4b9d-82e4-6057799dadae
function SurfaceResistance(gasData::GasData, G::T, Ts::T, θ::T,
	iSeason::Int, iLandUse::Int, rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool) where T<:AbstractFloat
	rs = r_s(G, Ts, iSeason, iLandUse, rain || dew)
	rmx = r_mx(gasData.Hstar, gasData.Fo)
	rsmx = r_smx(rs, gasData.Dh2oPerDx, rmx)
	rdc = r_dc(G, θ)
	rlux = r_lux(gasData.Hstar, gasData.Fo, iSeason, iLandUse,
		rain, dew, isSO2, isO3)
	rclx = 0.0
	rgsx = 0.0
	if isSO2
		rclx = r_clS[iSeason, iLandUse]
		rgsx = r_gsS[iSeason, iLandUse]
	elseif isO3
		rclx = r_clO[iSeason, iLandUse]
		rgsx = r_gsO[iSeason, iLandUse]
	else
		rclx = r_clx(gasData.Hstar, gasData.Fo, iSeason, iLandUse)
		rgsx = r_gsx(gasData.Hstar, gasData.Fo, iSeason, iLandUse)
	end
	
	rac = r_ac[iSeason, iLandUse]

	# Correction for cold temperatures from page 4 column 1.
	if Ts < 0.0
		correction = 1000.0 * exp(-Ts-4) # [s m-1] #mark
		rlux += correction
		rclx += correction
		rgsx += correction
	end
	
	r_c = 1.0 / (1.0/(rsmx) + 1.0/rlux + 1.0/(rdc+rclx) + 1.0/(rac+rgsx))
	r_c = max(r_c, 10.0) # From "Results and conclusions" section
	r_c = min(r_c, 9999.0)
	return r_c
end

# ╔═╡ dc996b49-c386-4048-a87a-0590dc37cd94
md"""
# Test case
"""

# ╔═╡ 859a71ad-21d6-4926-aca8-2b441e3aa5f5
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

# ╔═╡ caf60740-02e4-4a25-b3be-e1ebd4ffccab
function different(a::Float64, b::Float32) where T<:AbstractFloat
	c = abs(a - b)
	return c/b > 0.1 && c >= 11.0
end

# ╔═╡ 01701c34-d937-42eb-9dca-3513d7beda5c
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

# ╔═╡ 15600208-4263-4535-94d2-29f39de5b11f
TestWesely()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
StaticArrays = "~1.2.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
"""

# ╔═╡ Cell order:
# ╟─7363fcd7-b115-4c06-ab0f-429daaf56a83
# ╠═91119eae-3023-4097-bd50-8079e3ae69c7
# ╠═96a4a0c0-321b-11ec-09d7-3bce044720f0
# ╟─e655d9ac-8049-4950-831d-4ef4c421cdbb
# ╠═ef83b8f7-53fa-46b7-86de-2300532004b1
# ╟─9d86d682-92bd-4b26-ae23-c4a385374a7c
# ╠═301eaff9-06c8-4726-9067-27a5e813f86c
# ╟─9a7cb29b-85ba-4722-8bd4-329684837b83
# ╠═2b446731-9bbe-4f64-8c2c-e4f1dd7e5f41
# ╟─ad0cd7a5-5e67-4f0d-b1cb-9d9661e3548e
# ╠═c36b7c6e-648a-4bf8-8797-d4826c83ce2a
# ╠═1d2853eb-3a6d-4313-b6a1-4858a8708c1d
# ╟─15daae3e-8acd-4ee2-9b8f-f36a3e111c26
# ╠═e791cac6-4890-42d9-b770-99b9ff8d97c7
# ╠═8a78d619-8632-4fdd-800b-74f1a0b8a9d7
# ╟─431c0a99-9a71-4284-bc6a-2f95a879b3e1
# ╠═f160601e-e9e5-4181-8c4a-425c35838374
# ╟─adfb7cbf-3dc0-4fa7-9adc-66e203accffc
# ╠═34f744cb-d9fc-44b4-a61e-3262ee2b1e65
# ╟─31510e91-8f94-44ee-8163-2905249ba2fd
# ╠═5b0071b0-aab8-4800-a9ed-9dfe8947cb90
# ╠═fa59f3f5-1115-4329-9b1d-a6f943ec70aa
# ╠═344d4e1a-8c6e-4050-bdea-945b60d2aa6e
# ╠═4f47597b-b383-4ebd-9d1c-662c23ed8396
# ╠═c7c1aade-4eba-46e3-9d4a-e89a15a8d21b
# ╠═f0481175-a676-4df0-ad1e-aabcfefb0995
# ╠═09352c4a-0a0c-4640-b36d-c69b805ee24b
# ╠═0d24bfd4-bc76-4f0a-8814-a3e7cf1c9db1
# ╠═2e36c0e5-8037-4421-8370-2fbb49bdfe4a
# ╠═b2af8534-556c-4278-81f3-e6d9e81604ec
# ╠═fbfe0521-8abe-43e3-8453-4f2453dd635e
# ╠═31d79b1d-dc72-467d-9af0-4a7482124f47
# ╠═907b6856-0702-4b9d-82e4-6057799dadae
# ╟─dc996b49-c386-4048-a87a-0590dc37cd94
# ╠═859a71ad-21d6-4926-aca8-2b441e3aa5f5
# ╠═caf60740-02e4-4a25-b3be-e1ebd4ffccab
# ╠═01701c34-d937-42eb-9dca-3513d7beda5c
# ╠═15600208-4263-4535-94d2-29f39de5b11f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

end
