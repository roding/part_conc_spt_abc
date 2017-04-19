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
	number_of_frames::Array{Int64, 1} = 100 * ones(2000)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	
	de_number_of_bins::Int64 = 1250
	de_max::Float64 = 12.5
	
	# Simulate system.
	(n_K::Array{Int64, 1}, n_DE::Array{Int64, 1}) = simulate_system(
		distribution_class, 
		distribution_parameters, 
		c, 
		ax, 
		ay, 
		az, 
		Lx, 
		Ly, 
		Lz, 
		number_of_frames, 
		deltat,
		kmin, 
		de_number_of_bins, 
		de_max)
	
	# Save results.
	file_name_output = "simulated_system.dat"
	writedlm(file_name_output, [n_K n_DE], ',')
	
	nothing
end

test_simulate_system()