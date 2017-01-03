include("rand_poisson.jl")

function test_rand_poisson()
	lambda::Float64 = 10.0
	
	x = Array(Int64, 1000)
	
	for i = 1:length(x)
		x[i] = rand_poisson(lambda)
	end
	
	println(mean(x))
	
	nothing
end

test_rand_poisson()