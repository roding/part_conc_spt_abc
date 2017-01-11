workspace()

include("../src/simulate_system.jl")

function test_simulate_system()
	D::Float64 = 1.0 # µm^2/s.
	c::Float64 = 1e10 # part/ml.
	ax::Float64 = 80.0 # µm.
	ay::Float64 = 80.0 # µm.
	az::Float64 = 0.1 # µm.
	L::Float64 = 100.0 # µm.
	number_of_frames::Array{Int64, 1} = 100 * ones(2000)
	deltat::Float64 = 0.01 # seconds
	kmin::Int64 = 2
	
	(K, DE) = simulate_system_monodisperse(D, c, ax, ay, az, L, number_of_frames, deltat, kmin)
	
	#println(DE)
	file_name_output = "simulated_system.dat"
	writedlm(file_name_output, [K DE], ',')
	#file_stream_output = open(file_name_output, "w")
	#write(file_stream_output, K, DE)
	#close(file_stream_output)
	
	nothing
end

test_simulate_system()