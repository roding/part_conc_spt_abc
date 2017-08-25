function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	distribution_class::String = read_key(file_string, "distribution_class", String)
	number_of_components::Int64 = read_key(file_string, "number_of_components", Int64)
	number_of_abc_samples::Int64 = read_key(file_string, "number_of_abc_samples", Int64)
	m::Array{Float64, 2} = reshape(read_key(file_string, "m", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	s::Array{Float64, 2} = reshape(read_key(file_string, "s", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	c::Array{Float64, 2} = reshape(read_key(file_string, "c", Array{Float64, 1}), (number_of_components, number_of_abc_samples))
	az::Array{Float64, 1} = read_key(file_string, "az", Array{Float64, 1})
	dist::Array{Float64, 1} = read_key(file_string, "dist", Array{Float64, 1})
	w::Array{Float64, 1} = read_key(file_string, "w", Array{Float64, 1})
	epsilon::Float64 = read_key(file_string, "epsilon", Float64)

	return (
		distribution_class,
		number_of_components,
		number_of_abc_samples,
		m,
		s,
		c,
		az,
		dist,
		w,
		epsilon)
end
