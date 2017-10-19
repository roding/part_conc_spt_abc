function write_output(	file_path::String,
					model::String,
					number_of_components::Int64,
					number_of_abc_samples::Int64,
					m::Array{Float64, 2},
					s::Array{Float64, 2},
					c::Array{Float64, 2},
					az::Array{Float64, 1},
					dist::Array{Float64, 1},
					w::Array{Float64, 1},
					epsilon::Float64)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<output>\n")

	write_key(file_stream, "model", model)
	write_key(file_stream, "number_of_components", number_of_components)
	write_key(file_stream, "number_of_abc_samples", number_of_abc_samples)
	write_key(file_stream, "m", m[:])
	write_key(file_stream, "s", s[:])
	write_key(file_stream, "c", c[:])
	write_key(file_stream, "az", az)
	write_key(file_stream, "dist", dist)
	write_key(file_stream, "w", w)
	write_key(file_stream, "epsilon", epsilon)

	@printf(file_stream, "%s", "</output>")

	close(file_stream)

	nothing
end
