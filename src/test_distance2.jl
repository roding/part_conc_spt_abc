workspace()

include("distance2.jl")
include("simulate_system.jl")
include("rand_poisson.jl")

function test_distance2()
	
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	mu1::Float64 = 0.0 # µm^2/s.
	sigma1::Float64 = 1.0 # µm^2/s.
	mu2::Float64 = 0.0 # µm^2/s.
	sigma2::Float64 = 1.0 # µm^2/s.
	c1::Float64 = 1e8 # part/ml.
	c2::Float64 = 1e8 # part/ml.
	
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	az::Float64 = 5.0 # µm.
	L::Float64 = 100.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 3
	
	(K1, RSQ1) = simulate_system(mu1, sigma1, c1, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	(K2, RSQ2) = simulate_system(mu2, sigma2, c2, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	return distance2(K1, RSQ1, K2, RSQ2)
	
end

n = 100
x = 0.0
for i = 1:n
	x = x + test_distance2()
end
println(x/n)
