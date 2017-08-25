function distance(H::Array{Int64, 2}, H_sim::Array{Int64, 2})
	
	cum_H::Array{Float64, 2} = cumsum(cumsum(H, 1), 2)
	cum_H_sim::Array{Float64, 2} = cumsum(cumsum(H_sim, 1), 2)
	d::Float64 = sum( (cum_H - cum_H_sim).^2 )
	
	return d
end
