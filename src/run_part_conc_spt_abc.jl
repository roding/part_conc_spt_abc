workspace()

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
		return nothing
	end

	# Read input.
	


	nothing
end

run_part_conc_spt_abc()
