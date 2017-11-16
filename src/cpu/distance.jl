function distance(H::Array{Int64, 2}, H_sim::Array{Int64, 2})

	dH::Array{Float64, 2} = convert(Array{Float64, 2}, H - H_sim)
	cum_dH::Array{Float64, 2} = cumsum(cumsum(dH, 1), 2)
	d::Float64 = sum( cum_dH.^2 )

	return d
end
