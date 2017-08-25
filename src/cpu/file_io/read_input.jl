function read_input(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	data_file_path::String = read_key(file_string, "data_file_path", String)
	distribution_class::String = read_key(file_string, "distribution_class", String)
	number_of_components::Int64 = read_key(file_string, "number_of_components", Int64)
	Lx::Float64 = read_key(file_string, "Lx", Float64)
	Ly::Float64 = read_key(file_string, "Ly", Float64)
	Lz::Float64 = read_key(file_string, "Lz", Float64)
	kmin::Int64 = read_key(file_string, "kmin", Int64)
	number_of_de_bins::Int64 = read_key(file_string, "number_of_de_bins", Int64)
	ub_de::Float64 = read_key(file_string, "ub_de", Float64)
	lb_m::Float64 = read_key(file_string, "lb_m", Float64)
	ub_m::Float64 = read_key(file_string, "ub_m", Float64)
	lb_s::Float64 = read_key(file_string, "lb_s", Float64)
	ub_s::Float64 = read_key(file_string, "ub_s", Float64)
	lb_c::Float64 = read_key(file_string, "lb_c", Float64)
	ub_c::Float64 = read_key(file_string, "ub_c", Float64)
	lb_az::Float64 = read_key(file_string, "lb_az", Float64)
	ub_az::Float64 = read_key(file_string, "ub_az", Float64)
	number_of_abc_samples::Int64 = read_key(file_string, "number_of_abc_samples", Int64)
	gamma_initial::Float64 = read_key(file_string, "gamma_initial", Float64)
	delta_gamma::Float64 = read_key(file_string, "delta_gamma", Float64)
	output_file_path::String = read_key(file_string, "output_file_path", String)
	
	return (
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
end
