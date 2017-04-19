workspace()

include("../src/simulate_system.jl")

function test_simulate_system()
	# Diffusion coefficient distribution.
#	distribution_class::String = "monodisperse"
#	D::Float64 = 2.5 # µm^2/s.
#	distribution_parameters::Array{Float64, 1} = [D]
	
	distribution_class::String = "lognormal"
	m::Float64 = 0.5 # µm^2/s.
	s::Float64 = 0.05 # µm^2/s.
	distribution_parameters::Array{Float64, 1} = [m, s]
	
	# Concentration.
	c::Float64 = 5e8 # part/ml.
	
	# Acquisition parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	az::Float64 = 2.0 # µm.
	Lx::Float64 = 60.0 # µm.
	Ly::Float64 = 60.0 # µm.
	Lz::Float64 = 10.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 1#2
	kmax::Int64 = maximum(number_of_frames)
	
	de_number_of_bins::Int64 = 1000#1250
	de_max::Float64 = 10.0#12.5
	
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
	
	#println(n_K')
	#println(n_DE')
	
	println(sum(n_K))
	mean_K = sum(n_K .* (1:kmax)) / sum(n_K)
	println(mean_K)
	println(size(n_DE))
	println(size( de_max/de_number_of_bins:de_max/de_number_of_bins:de_max ))
	mean_DE = sum(n_DE .* ( de_max/de_number_of_bins:de_max/de_number_of_bins:de_max )) / sum(n_DE)	
	println(mean_DE)
	
	npart = sum(n_K) / sum(number_of_frames) * (Lx*Ly*Lz) / (ax*ay*az) 
	println(npart)
	# Save results.
	#file_name_output::String = "simulated_system.dat"
	#file_stream_output::IOStream = open(file_name_output, "w")
	#write(file_stream_output, n_K)
	#write(file_stream_output, n_DE)
	#close(file_stream_output)
	
	nothing
end

test_simulate_system()
