include("simulate_system.jl")

function test_simulate_system()
	mu::Float64 = 1.0 # µm^2/s.
	sigma::Float64 = 1.0 # µm^2/s.
	c::Float64 = 1e8 # part/ml.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	az::Float64 = 5.0 # µm.
	L::Float64 = 100.0 # µm.
	number_of_frames::Array{Int64, 1} = [1000]
	deltat::Float64 = 0.05 # seconds
	
	(K, R2) = simulate_system(mu, sigma, c, ax, ay, az, L, number_of_frames, deltat)
	
	println(K)
	println(R2)
	
end

test_simulate_system()