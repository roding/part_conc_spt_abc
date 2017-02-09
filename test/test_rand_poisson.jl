workspace()

include("../src/rand_poisson.jl")

function test_rand_poisson()
	lambda::Float64 = 10000.0
	
	n::Int64 = 1000
	x::Array{Int64, 1} = zeros(n)
	
	for i = 1:n
		x[i] = rand_poisson(lambda)
	end
	
	println(mean(convert(Array{Float64, 1}, x)))
	println(var(convert(Array{Float64, 1}, x)))
	
	nothing
end

test_rand_poisson()
