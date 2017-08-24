workspace()

include("../src/file_io/read_key.jl")
include("../src/file_io/read_input.jl")
include("../src/file_io/read_data.jl")
include("../src/file_io/write_key.jl")
include("../src/file_io/write_input.jl")
include("../src/file_io/write_data.jl")
include("../src/rand_component.jl")
include("../src/generate_experiment.jl")
include("../src/rand_poisson.jl")
include("../src/position_periodic.jl")

function test_part_conc_spt_abc()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Experimental parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	Lx::Float64 = 60.0 # µm.
	Ly::Float64 = 60.0 # µm.
	Lz::Float64 = 10.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(50)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	distribution_class::String = "discrete"
	D_real::Array{Float64, 1} = [1.5] # µm^2/s.
	w_real::Array{Float64, 1} = [1.0]
	c_real::Float64 = 1e8 # part/ml.
	az_real::Float64 = 2.0 # µm.
	
	# Simulate experiment.
	distribution_parameters::Array{Float64, 1} = vcat(D_real, w_real)
	(K::Array{Int64, 1}, DE::Array{Float64, 1}) = generate_experiment(	distribution_class, 
																distribution_parameters, 
																c_real, 
																ax, 
																ay, 
																az_real, 
																Lx, 
																Ly, 
																Lz, 
																number_of_frames, 
																deltat, 
																kmin)
	
	# Write data to file.
	data_file_path::String = "M:/part_conc_spt_abc/dev/part_conc_spt_abc/test/data.xml"
	write_data(	data_file_path,
				ax,
				ay,
				number_of_frames,
				deltat,
				K,
				DE)
				
	nothing
end

test_part_conc_spt_abc()