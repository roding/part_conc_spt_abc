function distance(K1::Array{Int64, 1}, RSQ1::Array{Float64, 1}, K2::Array{Int64, 1}, RSQ2::Array{Float64, 1})
	n1::Int64 = length(K1)
	n2::Int64 = length(K2)
	
	d::Float64 = 0.0
	p1::Float64 = 0.0
	p2::Float64 = 0.0
	
	# Evaluate 2-D ECDF for both data sets, at points in data set 1.
	for current_point_1 = 1:n1
		p1 = sum( (K1 .<= K1[current_point_1]) & (RSQ1 .<= RSQ1[current_point_1]) )
		p2 = sum( (K2 .<= K1[current_point_1]) & (RSQ2 .<= RSQ1[current_point_1]) )
		
		d = d + (p1 - p2)^2
	end

	# Evaluate 2-D ECDF for both data sets, at points in data set 2.
	#for current_point_2 = 1:n2
	#	p1 = sum( (K1 .<= K2[current_point_2]) & (RSQ1 .<= RSQ2[current_point_2]) )
	#	p2 = sum( (K2 .<= K2[current_point_2]) & (RSQ2 .<= RSQ2[current_point_2]) )
	#	
	#	d = d + (p1 - p2)^2
	#end
	
	d = d / convert(Float64, n1*n2)
	return d
end