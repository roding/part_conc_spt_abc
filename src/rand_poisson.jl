# Algorithm for generating Poisson distributed random numbers.
# Attributed to Donald Knuth, see Wikipedia.

function rand_poisson(lambda)
	x::Int64 = 0	
	if lambda <= 15.0	
		L::Float64 = exp(-lambda)
		p::Float64 = 1
    
		while p > L
		     x = x + 1
		     p = p * rand()
		end
		return x - 1
	else
		
end

x = zeros(1000000)
lambda = 5000.0
for i = 1:length(x)
	x[i] = rand_poisson(lambda)
end

println(mean(x))
