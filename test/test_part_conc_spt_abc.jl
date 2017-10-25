workspace()

include("../src/cpu/file_io/read_key.jl")
include("../src/cpu/file_io/read_input.jl")
include("../src/cpu/file_io/read_data.jl")
include("../src/cpu/file_io/write_key.jl")
include("../src/cpu/file_io/write_input.jl")
include("../src/cpu/file_io/write_data.jl")
include("../src/cpu/rand_weighted_index.jl")
include("../src/cpu/generate_experiment.jl")
include("../src/cpu/rand_poisson.jl")
include("../src/cpu/position_periodic.jl")

function test_part_conc_spt_abc()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Timing.
	time_start::Int64 = convert(Int64, time_ns())

	# Experimental parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	Lx::Float64 = 60.0 # µm.
	Ly::Float64 = 60.0 # µm.
	Lz::Float64 = 10.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(10)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	model::String = "discrete-fixed-depth"
	m_real::Array{Float64, 1} = [1.0, 3.0] # µm^2/s.
	c_real::Array{Float64, 1} = [1e8, 1e8] # part/ml.
	az_real::Array{Float64, 1} = [2.0, 2.0] # µm.

	# Simulate experiment.
	(K::Array{Int64, 1}, DE::Array{Float64, 1}) = generate_experiment(	model,
																		m_real,
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
	number_of_de_bins::Int64 = 200
	ub_de::Float64 = 4.0 * maximum(m_real)
	lb_m::Float64 = 0.25 * minimum(m_real)
	ub_m::Float64 = 4.0 * maximum(m_real)
	lb_c::Float64 = 0.25 * minimum(c_real)
	ub_c::Float64 = 4.0 * maximum(c_real)
	lb_az::Float64 = 0.25 * minimum(az_real)
	ub_az::Float64 = 4.0 * maximum(az_real)
	number_of_abc_samples::Int64 = 128
	gamma_initial::Float64 = 15.0
	gamma_adaptive::Bool = true#false
	delta_gamma::Float64 = 0.01
	weighting_scheme::String = "pmc-standard"
	ub_average_number_of_trials::Int64 = 50#500
	output_file_path::String = abspath("output.xml")
	write_input(input_file_path,
				data_file_path,
				model,
				number_of_components,
				Lx,
				Ly,
				Lz,
				kmin,
				number_of_de_bins,
				ub_de,
				lb_m,
				ub_m,
				lb_c,
				ub_c,
				lb_az,
				ub_az,
				number_of_abc_samples,
				gamma_initial,
				gamma_adaptive,
				delta_gamma,
				weighting_scheme,
				ub_average_number_of_trials,
				output_file_path)

	# Run inference.
	program_path::String = abspath("../src/cpu/run_part_conc_spt_abc.jl")
	#cmd::Cmd = `julia $program_path $input_file_path`
	str_number_of_cores::String = string(Sys.CPU_CORES)
	cmd::Cmd = `julia -p $str_number_of_cores $program_path $input_file_path`
	#cmd::Cmd = `julia -p 1 $program_path $input_file_path`
	#cmd::Cmd = `julia --math-mode=fast --check-bounds=no --optimize=3 -p $str_number_of_cores $program_path $input_file_path`
	run(cmd)

	time_end::Int64 = convert(Int64, time_ns())
	time_exec::Float64 = convert(Float64, time_end - time_start) / 1e9
	println(time_exec)
	nothing
end

test_part_conc_spt_abc()
