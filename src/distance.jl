function distance(n_K_real::Array{Int64, 1}, n_DE_real::Array{Int64, 1}, n_K_sim::Array{Int64, 1}, n_DE_sim::Array{Int64, 1})
	
#	d::Int64 = 0
#	
#	for i = 1:length(n_K_real)
#		d = d + (n_K_real[i] - n_K_sim[i])^2
#	end
#	
#	for i = 1:length(n_DE_real)
#		d = d + (n_DE_real[i] - n_DE_sim[i])^2
#	end
	
	d::Int64 = sum( (cumsum(n_K_real) - cumsum(n_K_sim)).^2 ) + sum( (cumsum(n_DE_real) - cumsum(n_DE_sim)).^2 )

	return convert(Float64, d)
end
