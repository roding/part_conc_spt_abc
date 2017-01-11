workspace()

include("../src/simulate_system.jl")

function test_simulate_system()
	# Diffusion coefficient distribution.
#	distribution_class::String = "monodisperse"
#	D::Float64 = 2.5 # µm^2/s.
#	distribution_parameters::Array{Float64, 1} = [D]
	
	distribution_class::String = "lognormal"
	m::Float64 = 2.5 # µm^2/s.
	s::Float64 = 0.5 # µm^2/s.
	mu::Float64 = log(m) - 0.5 * log(1 + s^2/m^2)
	sigma::Float64 = sqrt(log(1 + s^2/m^2))
	distribution_parameters::Array{Float64, 1} = [mu, sigma]
	
	# Concentration.
	c::Float64 = 1e10 # part/ml.
	
	# Acquisition parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	az::Float64 = 2.0 # µm.
	L::Float64 = 100.0 # µm.
	number_of_frames::Array{Int64, 1} = 100 * ones(20000)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	
	# Simulate system.
	(K, DE) = simulate_system(distribution_class, distribution_parameters, c, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	# Save results.
	file_name_output = "simulated_system.dat"
	writedlm(file_name_output, [K DE], ',')
	
	nothing
end

test_simulate_system()