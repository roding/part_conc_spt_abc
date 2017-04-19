function distance(n_K_real::Array{Int64, 1}, n_DE_real::Array{Int64, 1}, n_K_sim::Array{Int64, 1}, n_DE_sim::Array{Int64, 1})
	
	d::Int64 = sum( (cumsum(n_K_real) - cumsum(n_K_sim)).^2 ) + sum( (cumsum(n_DE_real) - cumsum(n_DE_sim)).^2 )

	return convert(Float64, d)
end
