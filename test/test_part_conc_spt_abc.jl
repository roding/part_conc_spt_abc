workspace()

include("../src/file_io/read_key.jl")
include("../src/file_io/read_input.jl")
include("../src/file_io/read_data.jl")
include("../src/file_io/write_key.jl")
include("../src/file_io/write_input.jl")
include("../src/file_io/write_data.jl")
include("../src/rand_weighted_index.jl")
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
	m_real::Array{Float64, 1} = [0.66666, 2.0] # µm^2/s.
	s_real::Array{Float64, 1} = [0.0, 0.0] # µm^2/s. Just put to zero for discrete model.
	c_real::Array{Float64, 1} = [1e8, 1e8] # part/ml.
	az_real::Float64 = 2.0 # µm.
	
	# Simulate experiment.
	(K::Array{Int64, 1}, DE::Array{Float64, 1}) = generate_experiment(	distribution_class, 
																m_real,
																s_real,
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
	data_file_path::String = abspath("data.xml")
	write_data(	data_file_path,
				ax,
				ay,
				number_of_frames,
				deltat,
				K,
				DE)
				
	# Write input to file.
	input_file_path::String = abspath("input.xml")
	number_of_components::Int64 = 2
	number_of_de_bins::Int64 = 2000
	ub_de::Float64 = 4.0 * maximum(m_real)
	lb_m::Float64 = 0.25 * minimum(m_real)
	ub_m::Float64 = 4.0 * maximum(m_real)
	lb_s::Float64 = 0.0 # Just put to zero for discrete model.
	ub_s::Float64 = 0.0 # Just put to zero for discrete model.
	lb_c::Float64 = 0.25 * minimum(c_real)
	ub_c::Float64 = 4.0 * maximum(c_real)
	lb_az::Float64 = 0.25 * az_real
	ub_az::Float64 = 4.0 * az_real
	number_of_abc_samples::Int64 = 512
	gamma_initial::Float64 = 9.0
	delta_gamma::Float64 = 0.01
	output_file_path::String = abspath("output.xml")
	write_input(	input_file_path,
				data_file_path,
				distribution_class,
				number_of_components,
				Lx,
				Ly,
				Lz,
				kmin,
				number_of_de_bins,
				ub_de,
				lb_m,
				ub_m,
				lb_s,
				ub_s,
				lb_c,
				ub_c,
				lb_az,
				ub_az,
				number_of_abc_samples,
				gamma_initial,
				delta_gamma,
				output_file_path)
	
	# Run inference.
	program_path::String = abspath("../src/run_part_conc_spt_abc.jl")
	#cmd::Cmd = `julia $program_path $input_file_path`
	cmd::Cmd = `julia -p 32 $program_path $input_file_path`
	run(cmd)
				
	nothing
end

test_part_conc_spt_abc()
