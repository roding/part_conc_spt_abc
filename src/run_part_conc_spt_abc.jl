workspace()

include("file_io/read_key.jl")
include("file_io/read_input.jl")
include("file_io/read_data.jl")
include("rand_poisson.jl")
include("position_periodic.jl")

function run_part_conc_spt_abc()
	# Start time.
	t_start_ns::Int64 = convert(Int64, time_ns())

	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Assert that input is file and store path.
	input_file_path::String = ""
	if isfile(ARGS[1])
		input_file_path = ARGS[1]
	else
		println("No input file specified or specified input file does not exist. Aborting.")
		return -1
	end

	# Read input.
	(	data_file_path::String,
		distribution_class::String,
		number_of_components::Int64,
		Lx::Float64,
		Ly::Float64,
		Lz::Float64,
		kmin::Int64,
		number_of_de_bins::Int64,
		ub_de::Float64,
		lb_m::Float64,
		ub_m::Float64,
		lb_s::Float64,
		ub_s::Float64,
		lb_c::Float64,
		ub_c::Float64,
		lb_az::Float64,
		ub_az::Float64,
		number_of_abc_samples::Int64,
		gamma_initial::Float64,
		delta_gamma::Float64,
		output_file_path::String) = read_input(input_file_path)
		
	# Read data.
	(	ax::Float64,
		ay::Float64,
		number_of_frames::Array{Int64, 1},
		deltat::Float64,
		K::Array{Int64, 1},
		DE::Array{Float64, 1}) = read_data(data_file_path)
		
	


	return 0
end

run_part_conc_spt_abc()
