include("simulate_system.jl")

function test_simulate_system()
	mu::Float64
	sigma::Float64
	c::Float64
	ax::Float64
	ay::Float64
	az::Float64
	L::Float64
	number_of_frames::Array{Int64, 1}
	deltat::Float64
	
	(K::Array{Int64, 1}, R2::Array{Float64, 1}) = simulate_system(mu, sigma, c, ax, ay, az, L, number_of_frames, deltat)
	
	println((K,R2))
	
end

test_simulate_system()