function rand_poisson(lambda)
	x::Int64 = 0	
	t::Float64 = 0.0
	
	while t < lambda
		x = x + 1
		t = t + randexp()
	end
	
	return x - 1
end