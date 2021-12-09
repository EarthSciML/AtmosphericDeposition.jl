module DepositionMTK
	using StaticArrays
	using Test
	
	# Calculate bulk canopy stomatal resistance [s m-1] based on Wesely (1989) equation 3 when given the solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]), the season index (iSeason), the land use index (iLandUse), and whether there is currently rain or dew.
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

	# Calculate the resistance from the effects of mixing forced by buoyant convection when sunlight heats the ground or lower canopy and by penetration of wind into canopies on the sides of hills [s m-1] when given the solar irradiation (G [W m-2]) and the slope of the local terrain (θ [radians]). From Wesely (1989) equation 5.
	function r_dc(G::T, θ::T) where T<:AbstractFloat
		return 100.0 * (1.0 + 1000.0/(G+10.0)) / (1.0 + 1000.0*θ)
	end

	#Calculate mesophyll resistance [s m-1] based on Wesely (1989) equation 6 when given the effective Henry's law coefficient (Hstar [M atm-1]) and the reactivity factor (fo [-]), both available in Wesely (1989) table 2.
	function r_mx(Hstar::T, fo::T) where T<:AbstractFloat
		return 1.0 / (Hstar/3000.0 + 100.0*fo)
	end

	#Calculate combined minimum stomatal and mesophyll resistance [s m-1] based on Wesely (1989) equation 4 when given stomatal resistance (r_s [s m-1]), ratio of water to chemical-of-interest diffusivities (Dh2oPerDx [-]), and mesophyll resistance (r_mx [s m-1]).
	function r_smx(r_s::T, Dh2oPerDx::T, r_mx::T) where T<:AbstractFloat
		return r_s * Dh2oPerDx + r_mx
	end

	# Calculate the resistance of the outer surfaces in the upper canopy (leaf cuticular resistance in healthy vegetation) based on Wesely (1989) equations 7 and 10-14 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), the land use index (iLandUse), whether there is currently rain or dew, and whether the chemical of interest is either SO2 or O3.
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

	# Calculate the resistance of the exposed surfaces in the lower portions of structures (canopies, buildings) above the ground based on Wesely (1989) equation 8 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
	function r_clx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
		return 1.0 / (Hstar/(1.0e5*r_clS[iSeason, iLandUse]) +
			fo/r_clO[iSeason, iLandUse])
	end

	# Calculate the resistance to uptake at the 'ground' surface based on Wesely (1989) equation 9 when given the effective Henry's law coefficient (Hstar [M atm-1]), the reactivity factor (fo [-]) (both available in Wesely (1989) table 2), the season index (iSeason), and the land use index (iLandUse).
	function r_gsx(Hstar::T, fo::T, iSeason::Int, iLandUse::Int) where T<:AbstractFloat
		return 1.0 / (Hstar/(1.0e5*r_gsS[iSeason, iLandUse]) +
			fo/r_gsO[iSeason, iLandUse])
	end

	# Calculates surface resistance to dry depostion [s m-1] based on Wesely (1989) equation 2 when given information on the chemical of interest (gasData), solar irradiation (G [W m-2]), the surface air temperature (Ts [°C]), the slope of the local terrain (Θ [radians]), the season index (iSeason), the land use index (iLandUse), whether there is currently rain or dew, and whether the chemical of interest is either SO2 (isSO2) or O3 (isO3).
	# From Wesely (1989) regarding rain and dew inputs: "A direct computation of the surface wetness would be most desirable, e.g. by estimating the amount of free surface water accumulated and then evaporated. Alternatively, surface relative humidity might be a useful	index. After dewfall and rainfall events are completed, surface wetness	often disappears as a result of evaporation after approximately 2	hours of good atmospheric mixing, the period of time recommended earlier (Sheih et al., 1986)".
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
end
