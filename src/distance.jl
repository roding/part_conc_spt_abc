function distance(n_K_real::Array{Int64, 1}, n_DE_real::Array{Float64, 1}, n_K_sim::Array{Int64, 1}, n_DE_sim::Array{Float64, 1})
	
	d::Float64 = sum( (n_K_real - n_K_sim).^2 ) + sum( (n_DE_real - n_DE_sim).^2 )

	return d
end