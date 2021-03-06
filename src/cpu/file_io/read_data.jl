function read_data(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	ax::Float64 = read_key(file_string, "ax", Float64)
	ay::Float64 = read_key(file_string, "ay", Float64)
	number_of_frames::Array{Int64, 1} = read_key(file_string, "number_of_frames", Array{Int64, 1})
	deltat::Float64 = read_key(file_string, "deltat", Float64)
	K::Array{Int64, 1} = read_key(file_string, "K", Array{Int64, 1})
	DE::Array{Float64, 1} = read_key(file_string, "DE", Array{Float64, 1})
	
	return (
		ax,
		ay,
		number_of_frames,
		deltat,
		K,
		DE)
end
