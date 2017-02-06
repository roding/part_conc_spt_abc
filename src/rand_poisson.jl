# Algorithm for generating Poisson distributed random numbers.
# Attributed to Donald Knuth, see Wikipedia.

function rand_poisson(lambda)
	x::Int64 = 0	
#	if lambda <= 15.0	
#		L::Float64 = exp(-lambda)
#		p::Float64 = 1
#   
#		while p > L
#		     x = x + 1
#		     p = p * rand()
#		end
#		return x - 1
#	else
#		
#	end

	t::Float64 = 0.0
	while t < lambda
		x = x + 1
		t = t + randexp()
	end
	
	return x - 1
	
end
