workspace()

include("../src/rand_poisson.jl")

function test_rand_poisson()
	lambda::Float64 = 100000.0
	
	x::Int64 = 0
	n::Int64 = 10000
	for i = 1:n
		x = x + rand_poisson(lambda)
	end
	
	println(x / convert(Float64, n))
	
	nothing
end

test_rand_poisson()