function histogram_de(DE::Array{Float64, 1}, de_bin_edges::FloatRange{Float64})
	n_DE::Array{Int64, 1} = zeros(length(de_bin_edges) - 1)
	
	
	
	#x = x - x_bin_edges[1]
	#x_bin_edges = x_bin_edges - x_bin_edges[1]
	#println(x)
	x = x / dx
	
	x = convert(Array{Int64, 1}, ceil(x))
	#println(x)
	for current_x = 1:length(x)
		if x[current_x] == 0
			x[current_x] = 1
		elseif x[current_x] == length(n) + 1
			x[current_x] = length(n)
		end
		
		n[x[current_x]] = n[x[current_x]] + 1
	end
	
	return n
	
end