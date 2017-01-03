include("distance.jl")
include("simulate_system.jl")

function test_distance()
	mu::Float64 = 1.0 # µm^2/s.
	sigma::Float64 = 1.0 # µm^2/s.
	c::Float64 = 1e8 # part/ml.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	az::Float64 = 5.0 # µm.
	L::Float64 = 100.0 # µm.
	number_of_frames::Array{Int64, 1} = [1000]
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 3
	
	(K1, RSQ1) = simulate_system(mu, sigma, c, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	(K2, RSQ2) = simulate_system(mu, sigma, c, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	println(distance(K1, RSQ1, K2, RSQ2))
	
end

test_distance()
