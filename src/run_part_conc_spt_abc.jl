workspace()

include("file_io/read_key.jl")
include("file_io/read_input.jl")
include("file_io/read_data.jl")
include("estimate.jl")

foo = @__FILE__
@eval @everywhere f = $foo
@everywhere (program_file_dir, program_file_name) = splitdir(f)
@everywhere include(joinpath(program_file_dir, "simulate.jl"))
@everywhere include(joinpath(program_file_dir, "distance.jl"))
@everywhere include(joinpath(program_file_dir, "rand_weighted_index.jl"))
@everywhere include(joinpath(program_file_dir, "position_periodic.jl"))
@everywhere include(joinpath(program_file_dir, "displace.jl"))
@everywhere include(joinpath(program_file_dir, "rand_poisson.jl"))

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
		
	println(distribution_class)
	println(number_of_components)
	println(Lx)
	println(Ly)
	println(Lz)
	println(kmin)
	println(number_of_de_bins)
	println((ub_de,lb_m,ub_m,lb_s,ub_s,lb_c,ub_c,lb_az,ub_az,number_of_abc_samples,gamma_initial,delta_gamma))
	println(output_file_path)
	println(ax)
	println(ay)
	println(number_of_frames)
	println(deltat)
	println(length(K))
	println(length(DE))
	
	# Run ABC.
	(	m::Array{Float64, 2}, 
		s::Array{Float64, 2}, 
		c::Array{Float64, 2}, 
		az::Array{Float64, 1}, 
		dist::Array{Float64, 1}, 
		w::Array{Float64, 1}, 
		epsilon::Float64) = estimate(	distribution_class,
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
									ax,
									ay,
									number_of_frames,
									deltat,
									K,
									DE)

	return 0
end

run_part_conc_spt_abc()
