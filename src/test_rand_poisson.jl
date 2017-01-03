workspace()

include("rand_poisson.jl")

function test_rand_poisson()
	lambda::Float64 = 2.5
	
	x::Int64 = 0
	n::Int64 = 100000000
	for i = 1:n
		x = x + rand_poisson(lambda)
	end
	
	println(x / convert(Float64, n))
	
	nothing
end

test_rand_poisson()