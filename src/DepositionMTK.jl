module DepositionMTK
using Markdown
using InteractiveUtils

# ╔═╡ 2e6e8182-5253-11ec-1249-e16119693727
begin
	using StaticArrays
	using Test
end

# ╔═╡ 2fc50067-d03f-412b-951f-9619e9cc66fb
const inf = 1.e25

# ╔═╡ 0583dbd7-95e2-4c17-a6ba-2a0b66db7195
md"""
r_i represents the minimum bulk canopy resistances for water vapor.
"""

# ╔═╡ 9fdeb749-1774-432b-90da-8410443de455
const r_i = SA_F32[
	inf 60 120 70 130 100 inf inf 80 100 150
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 250 500 inf inf inf inf inf
	inf inf inf inf 400 800 inf inf inf inf inf
	inf 120 240 140 250 190 inf inf 160 200 300]

# ╔═╡ 84a56271-9589-4626-9c86-eef331b72692
md"""
r_lu signifies leaf cuticles in healthy vegetation and otherwise the outer surfaces in the upper canopy.
"""

# ╔═╡ bd011600-1956-4b8c-ad93-222c2ce36c33
const r_lu = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf 9000 9000 4000 8000 inf inf 9000 9000 9000
	inf inf inf inf 6000 9000 inf inf 9000 9000 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

# ╔═╡ 3861e155-f082-429f-bc5a-9bcc5ed0a2fc
md"""
r_ac signifies transfer that depends only on canopy height and density.
"""

# ╔═╡ 9d49357e-3a94-4065-a58a-716c58b38898
const r_ac = SA_F32[
	100.0 200 100 2000 2000 2000 0 0 300 150 200
	100 150 100 1500 2000 1700 0 0 200 120 140
	100 10 100 1000 2000 1500 0 0 100 50 120
	100 10 10 1000 2000 1500 0 0 50 10 50
	100 50 80 1200 2000 1500 0 0 200 60 120]

# ╔═╡ d3fbe7d7-b432-42a8-81e0-e7d2b13066d5
md"""
r_gs signifies uptake at the "ground" by soil, leaf litter, snow, water etc. 'S' and 'O' stand for SO2 and O3 respectively.
"""

# ╔═╡ ddd27845-4a38-44e8-8261-ceca29c19af6
const r_gsS = SA_F32[
	400.0 150 350 500 500 100 0 1000 0 220 40
	400 200 350 500 500 100 0 1000 0 300 400
	400 150 350 500 500 200 0 1000 0 200 400
	100 100 100 100 100 100 0 1000 100 100 50
	500 150 350 500 500 200 0 1000 0 250 40]

# ╔═╡ 4a000b1a-0aa5-46ee-a5b7-25158a73d35f
const r_gsO = SA_F32[
	300.0 150 200 200 200 300 2000 400 1000 180 200
	300 150 200 200 200 300 2000 400 800 180 200
	300 150 200 200 200 300 2000 400 1000 180 20
	600 3500 3500 3500 3500 3500 2000 400 3500 3500 3500
	300 150 200 200 200 300 2000 400 1000 180 200]

# ╔═╡ 669cab8e-a55f-40b3-aa14-7b67c9941ac4
md"""
r_cl is meant to account for uptake pathways at the leaves, bark, etc. 'S' and 'O' stand for SO2 and O3 respectively.
"""

# ╔═╡ 59048b8f-ef68-4c04-8e1c-ee0da4c5e892
const r_clS = SA_F32[
	inf 2000 2000 2000 2000 2000 inf inf 2500 2000 4000
	inf 9000 9000 9000 2000 4000 inf inf 9000 9000 9000
	inf inf 9000 9000 3000 6000 inf inf 9000 9000 9000
	inf inf inf 9000 200 400 inf inf 9000 inf 9000
	inf 4000 4000 4000 2000 3000 inf inf 4000 4000 8000]

# ╔═╡ dfc0f57a-03ca-4cac-872a-e78fe6a7d4c2
const r_clO = SA_F32[
	inf 1000 1000 1000 1000 1000 inf inf 1000 1000 1000
	inf 400 400 400 1000 600 inf inf 400 400 400
	inf 1000 400 400 1000 600 inf inf 800 600 600
	inf 1000 1000 400 1500 600 inf inf 800 1000 800
	inf 1000 500 500 1500 700 inf inf 600 800 800]

# ╔═╡ 52cc1697-fce3-4e31-b530-80121aa8ef00
md"""
Holder for gas properties from Wesely (1989) Table 2.
"""

# ╔═╡ 35c53bc4-7816-40ae-976a-2f774593802e
struct GasData
	Dh2oPerDx::AbstractFloat
	Hstar::AbstractFloat
	Fo::AbstractFloat
end

# ╔═╡ 0146649e-bd35-414d-bcb0-34f9d0aceae7
md"""
Properties of various gases from Wesely (1989) Table 2.
"""

# ╔═╡ 6b385d2e-420b-4250-98b7-a89be83f73db
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

# ╔═╡ 505b1901-06c4-4e35-8ef1-844627725e6a
begin
	const Midsummer =  1 		# 0: Midsummer with lush vegetation
	const Autumn = 2        	# 1: Autumn with unharvested cropland
	const LateAutumn = 3    	# 2: Late autumn after frost, no snow
	const Winter = 4        	# 3: Winter, snow on ground and subfreezing
	const Transitional = 5  	# 4: Transitional spring with partially green short 								     annuals
end

# ╔═╡ 744f73f4-36c8-4540-8fcb-61da1f236635
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

# ╔═╡ baa3c332-6c4c-4e25-b4c9-3679defe4793
md"""
Calculate bulk canopy stomatal resistance [s m-1] based on Wesely (1989) equation 3 when given the solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]), the season index (iSeason), the land use index (iLandUse), and whether there is currently rain or dew.
"""

# ╔═╡ 9c002cf1-fc6e-4750-af52-dd83a13cfe0c
function r_s(G::T, Ts::T, iSeason::Int, iLandUse::Int, rainOrDew::Bool) where T<:AbstractFloat
	rs = 0.0
	if Ts >= 39.9 || Ts <= 0.1
		rs = inf
	else
		rs = r_i[iSeason, iLandUse] * (1 + (200.0 * 1.0 / (G + 1.0))^2) *
			(400.0 * 1.0 / (Ts * (40.0 - Ts)))
	end
	# Adjust for dew and rain (from "Effects of dew and rain" section).
	if rainOrDew
		rs *= 3
	end
	return rs
end

# ╔═╡ 3731fea0-3b03-476c-a4dd-ac0c8e6cc81e
md"""
Calculate the resistance from the effects of mixing forced by buoyant convection when sunlight heats the ground or lower canopy and by penetration of wind into canopies on the sides of hills [s m-1] when given the solar irradiation (G [W m-2]) and the slope of the local terrain (θ [radians]). From Wesely (1989) equation 5.
"""

# ╔═╡ 988007aa-a0d2-431c-b220-e7ae983b5cfd
function r_dc(G::T, θ::T) where T<:AbstractFloat
	return 100.0 * (1.0 + 1000.0/(G+10.0)) / (1.0 + 1000.0*θ)
end

# ╔═╡ 0988dc3a-5687-4552-8a51-524e74ae4a5a
md"""
Calculate mesophyll resistance [s m-1] based on Wesely (1989) equation 6 when given the effective Henry's law coefficient (Hstar [M atm-1]) and the reactivity factor (fo [-]), both available in Wesely (1989) table 2.
"""

# ╔═╡ 87f2e68f-6c3e-49f1-9283-aec8384994b5
function r_mx(Hstar::T, fo::T) where T<:AbstractFloat
	return 1.0 / (Hstar/3000.0 + 100.0*fo)
end

# ╔═╡ c5f2979b-29fc-4c19-a1e9-bf9a070393af
md"""
Calculate combined minimum stomatal and mesophyll resistance [s m-1] based on Wesely (1989) equation 4 when given stomatal resistance (r_s [s m-1]), ratio of water to chemical-of-interest diffusivities (Dh2oPerDx [-]), and mesophyll resistance (r_mx [s m-1]).
"""

# ╔═╡ d4d1128c-9200-488c-b055-31122b54b9a3
function r_smx(r_s::T, Dh2oPerDx::T, r_mx::T) where T<:AbstractFloat
	return r_s * Dh2oPerDx + r_mx
end

# ╔═╡ 2d676758-2fb6-4295-9997-b2f8aab98184
md"""
Calculate the resistance of the outer surfaces in the upper canopy (leaf cuticular resistance in healthy vegetation) based on Wesely (1989) equations 7 and 10-14 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), the land use index (iLandUse), whether there is currently rain or dew, and whether the chemical of interest is either SO2 or O3.
"""

# ╔═╡ 00303208-de33-4afe-9fad-520a51a6e5df
function r_lux(Hstar::T, fo::T, iSeason::Int, iLandUse::Int, rain::Bool, dew::Bool, isSO2::Bool, isO3::Bool) where T<:AbstractFloat
	rlux = 0.0
	if dew && (iSeason != 4) # Dew doesn't have any effect in the winter
		if isSO2
			if iLandUse == 1
				rlux = 50.0 # equation 13 and a half
			else
				rlux = 100 # equation 10.
			end
		elseif isO3
			# equation 11
			rlux = 1.0 / (1.0/3000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
		else
			rluO = 1.0 / (1.0/3000.0 + 1.0/(3*r_lu[iSeason, iLandUse])) #equation 11
			rlux = 1.0 / (1.0/(3*r_lu[iSeason, iLandUse]/(1.0e-5*Hstar+fo)) 
				+ 1.0e-7*Hstar + fo/rluO) # equation 14, modified to match Walmsley eq. 5g
		end
	elseif rain && (iSeason != 4)
		if isSO2
			if iLandUse == 1
				rlux = 50 #equation 13 and a half
			else
				# equation 12
				rlux = 1.0 / (1.0/5000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
			end
		elseif isO3
			# equation 13
			rlux = 1.0 / (1.0/1000.0 + 1.0/(3*r_lu[iSeason, iLandUse]))
		else 
			rluO = 1.0 / (1.0/1000.0 + 1.0/(3*r_lu[iSeason, iLandUse])) # equation 13
			rlux = 1.0 / (1.0/(3*r_lu[iSeason, iLandUse]/(1.e-5*Hstar+fo)) + 1.0e-7*Hstar +fo/rluO) # equation 14, modified to match Walmsley eq. 5g
		end				
	else
		rlux = r_lu[iSeason, iLandUse] / (1.0e-5*Hstar + fo)
	end
	return rlux
end

# ╔═╡ 1f4b159b-fe58-4b86-aa69-aa6130977bbb
md"""
Calculate the resistance of the exposed surfaces in the lower portions of structures (canopies, buildings) above the ground based on Wesely (1989) equation 8 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
"""

# ╔═╡ 8abc6c93-0496-4ec2-9a1d-8fd4693bfb0d
function r_clx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
	return 1.0 / (Hstar/(1.0e5*r_clS[iSeason, iLandUse]) +
		fo/r_clO[iSeason, iLandUse])
end

# ╔═╡ 6e05e60f-7a9c-4fdf-9e8a-5054289fd355
md"""
Calculate the resistance to uptake at the 'ground' surface based on Wesely (1989) equation 9 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
"""

# ╔═╡ a72dc26e-45d1-460d-bf43-cf7150b15bd8
function r_gsx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
	return 1.0 / (Hstar/(1.0e5*r_gsS[iSeason, iLandUse]) +
		fo/r_gsO[iSeason, iLandUse])
end

# ╔═╡ 9922c06a-abe0-4f91-82e6-92d2001212f9
md"""
Calculates surface resistance to dry depostion [s m-1] based on Wesely (1989)
equation 2 when given information on the chemical of interest (gasData),
solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]),
the slope of the local terrain (Θ [radians]),
the season index (iSeason), the land use index (iLandUse), whether there is
currently rain or dew, and whether the chemical of interest is either SO2
(isSO2) or O3 (isO3).

From Wesely (1989) regarding rain and dew inputs:
	"A direct computation of the surface wetness would be most desirable, e.g.
	by estimating the amount of free surface water accumulated and then
	evaporated. Alternatively, surface relative humidity might be a useful
	index. After dewfall and rainfall events are completed, surface wetness
	often disappears as a result of evaporation after approximately 2
	hours of good atmospheric mixing, the period of time recommended earlier
	(Sheih et al., 1986)".
"""

# ╔═╡ 421a3505-e1dd-43a3-8708-589c0b1018c5
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
	r_c = max(r_c, 10.0) # From "Results and conclusions" section to avoid extremely high deposition velocities over extremely rough surfaces.
	r_c = min(r_c, 9999.0)
	return r_c
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
StaticArrays = "~1.2.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

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

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ Cell order:
# ╠═2e6e8182-5253-11ec-1249-e16119693727
# ╠═2fc50067-d03f-412b-951f-9619e9cc66fb
# ╟─0583dbd7-95e2-4c17-a6ba-2a0b66db7195
# ╠═9fdeb749-1774-432b-90da-8410443de455
# ╟─84a56271-9589-4626-9c86-eef331b72692
# ╠═bd011600-1956-4b8c-ad93-222c2ce36c33
# ╟─3861e155-f082-429f-bc5a-9bcc5ed0a2fc
# ╠═9d49357e-3a94-4065-a58a-716c58b38898
# ╟─d3fbe7d7-b432-42a8-81e0-e7d2b13066d5
# ╠═ddd27845-4a38-44e8-8261-ceca29c19af6
# ╠═4a000b1a-0aa5-46ee-a5b7-25158a73d35f
# ╟─669cab8e-a55f-40b3-aa14-7b67c9941ac4
# ╠═59048b8f-ef68-4c04-8e1c-ee0da4c5e892
# ╠═dfc0f57a-03ca-4cac-872a-e78fe6a7d4c2
# ╟─52cc1697-fce3-4e31-b530-80121aa8ef00
# ╠═35c53bc4-7816-40ae-976a-2f774593802e
# ╠═0146649e-bd35-414d-bcb0-34f9d0aceae7
# ╠═6b385d2e-420b-4250-98b7-a89be83f73db
# ╠═505b1901-06c4-4e35-8ef1-844627725e6a
# ╠═744f73f4-36c8-4540-8fcb-61da1f236635
# ╟─baa3c332-6c4c-4e25-b4c9-3679defe4793
# ╠═9c002cf1-fc6e-4750-af52-dd83a13cfe0c
# ╟─3731fea0-3b03-476c-a4dd-ac0c8e6cc81e
# ╠═988007aa-a0d2-431c-b220-e7ae983b5cfd
# ╟─0988dc3a-5687-4552-8a51-524e74ae4a5a
# ╠═87f2e68f-6c3e-49f1-9283-aec8384994b5
# ╟─c5f2979b-29fc-4c19-a1e9-bf9a070393af
# ╠═d4d1128c-9200-488c-b055-31122b54b9a3
# ╟─2d676758-2fb6-4295-9997-b2f8aab98184
# ╠═00303208-de33-4afe-9fad-520a51a6e5df
# ╟─1f4b159b-fe58-4b86-aa69-aa6130977bbb
# ╠═8abc6c93-0496-4ec2-9a1d-8fd4693bfb0d
# ╟─6e05e60f-7a9c-4fdf-9e8a-5054289fd355
# ╠═a72dc26e-45d1-460d-bf43-cf7150b15bd8
# ╟─9922c06a-abe0-4f91-82e6-92d2001212f9
# ╠═421a3505-e1dd-43a3-8708-589c0b1018c5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002


end
