function write_data(	file_path::String,
					ax::Float64,
					ay::Float64,
					number_of_frames::Array{Int64, 1},
					deltat::Float64,
					K::Array{Int64, 1},
					DE::Array{Float64, 1})

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<data>\n")

	write_key(file_stream, "ax", ax)
	write_key(file_stream, "ay", ay)
	write_key(file_stream, "number_of_frames", number_of_frames)
	write_key(file_stream, "deltat", deltat)
	write_key(file_stream, "K", K)
	write_key(file_stream, "DE", DE)


	@printf(file_stream, "%s", "</data>")

	close(file_stream)

	nothing
end
