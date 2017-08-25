function write_input(	file_path::String,
					data_file_path::String,
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
					ub_average_number_of_trials::Int64,
					output_file_path::String)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<input>\n")

	write_key(file_stream, "data_file_path", data_file_path)
	write_key(file_stream, "distribution_class", distribution_class)
	write_key(file_stream, "number_of_components", number_of_components)
	write_key(file_stream, "Lx", Lx)
	write_key(file_stream, "Ly", Ly)
	write_key(file_stream, "Lz", Lz)
	write_key(file_stream, "kmin", kmin)
	write_key(file_stream, "number_of_de_bins", number_of_de_bins)
	write_key(file_stream, "ub_de", ub_de)
	write_key(file_stream, "lb_m", lb_m)
	write_key(file_stream, "ub_m", ub_m)
	write_key(file_stream, "lb_s", lb_s)
	write_key(file_stream, "ub_s", ub_s)
	write_key(file_stream, "lb_c", lb_c)
	write_key(file_stream, "ub_c", ub_c)
	write_key(file_stream, "lb_az", lb_az)
	write_key(file_stream, "ub_az", ub_az)
	write_key(file_stream, "number_of_abc_samples", number_of_abc_samples)
	write_key(file_stream, "gamma_initial", gamma_initial)
	write_key(file_stream, "delta_gamma", delta_gamma)
	write_key(file_stream, "ub_average_number_of_trials", ub_average_number_of_trials)
	write_key(file_stream, "output_file_path", output_file_path)

	@printf(file_stream, "%s", "</input>")

	close(file_stream)

	nothing
end
