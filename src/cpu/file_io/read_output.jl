function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	model::String = read_key(file_string, "model", String)
	number_of_components::Int64 = read_key(file_string, "number_of_components", Int64)
	number_of_abc_samples::Int64 = read_key(file_string, "number_of_abc_samples", Int64)
	m::Array{Float64, 2} = reshape(read_key(file_string, "m", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	c::Array{Float64, 2} = reshape(read_key(file_string, "c", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	az::Array{Float64, 2} = reshape(read_key(file_string, "az", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	dist::Array{Float64, 1} = read_key(file_string, "dist", Array{Float64, 1})
	w::Array{Float64, 1} = read_key(file_string, "w", Array{Float64, 1})
	weighting_scheme::String = read_key(file_string, "weighting_scheme", String)
	gamma::Float64 = read_key(file_string, "gamma", Float64)
	t_exec::Float64 = read_key(file_string, "t_exec", Float64)
	number_of_iterations::Int64 = read_key(file_string, "number_of_iterations", Int64)
	number_of_simulations::Int64 = read_key(file_string, "number_of_simulations", Int64)

	return (
		model,
		number_of_components,
		number_of_abc_samples,
		m,
		c,
		az,
		dist,
		w,
		weighting_scheme,
		gamma,
		t_exec,
		number_of_iterations,
		number_of_simulations)
end
