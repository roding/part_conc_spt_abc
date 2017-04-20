function distance(n_K_real::Array{Int64, 1}, n_DE_real::Array{Int64, 1}, n_K_sim::Array{Int64, 1}, n_DE_sim::Array{Int64, 1}, d_de::Float64)
	
	d::Float64 = maximum( (cumsum(n_K_real) - cumsum(n_K_sim)).^2 ) + maximum( (cumsum(n_DE_real) - cumsum(n_DE_sim)).^2)
	#d::Float64 = sum( (cumsum(n_K_real) - cumsum(n_K_sim)).^2 ) + sum( (cumsum(n_DE_real) - cumsum(n_DE_sim)).^2 * d_de )
	#println( ( sum( (cumsum(n_K_real) - cumsum(n_K_sim)).^2 ), sum( (cumsum(n_DE_real) - cumsum(n_DE_sim)).^2 ) * d_de ) )
	return d
end
