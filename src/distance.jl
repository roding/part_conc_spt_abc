function distance(K1::Array{Int64, 1}, DE1::Array{Float64, 1}, K2::Array{Int64, 1}, DE2::Array{Float64, 1})
	n1::Int64 = length(K1)
	n2::Int64 = length(K2)
	
	d::Float64 = 0.0
	
	#sort!(K1)
	#sort!(K2)
	#sort!(DE1)
	#sort!(DE2)
	
#	if n1 <= n2
#		d = d + sum((K1 - K2[1:n1]).^2) + sum(K2[n1+1:end].^2)
#	else
#		d = d + sum((K1[1:n2] - K2).^2) + sum(K1[n2+1:end].^2)
#	end
#	if n1 <= n2
#		d = d + sum((DE1 - DE2[1:n1]).^2) + sum(DE2[n1+1:end].^2)
#	else
#		d = d + sum((DE1[1:n2] - DE2).^2) + sum(DE1[n2+1:end].^2)
#	end

#	# Distance for D inference only.
#	if length(DE1) > 0 && length(DE2) > 0
#		d = abs( mean(DE1) - mean(DE2) )
#	else
#		#println(length(DE1))
#		#println(length(DE2))
#		d = Inf
#	end

	# Distance for c inference only.
#	if length(K1) > 0 && length(K2) > 0
#		d = abs( sum(K1) - sum(K2) )
#	else
#		d = Inf
#	end
	
	# Distance for az inference only.
#	if length(K1) > 0 && length(K2) > 0
#		d = abs( mean(K1) - mean(K2) )
#	else
#		d = Inf
#	end
	
	if length(K1) > 0 && length(K2) > 0 # Then also D1 and D2 are non-empty.
		d = d + abs( mean(DE1) - mean(DE2) ) / ( 0.5 * (mean(DE1) + mean(DE2)) )
		d = d + abs( std(DE1) - std(DE2) ) / ( 0.5 * (std(DE1) + std(DE2)) )
		d = d + abs( mean(K1) - mean(K2) ) / ( 0.5 * (mean(K1) + mean(K2)) )
		d = d + abs( sum(K1) - sum(K2) ) / ( 0.5 * (sum(K1) + sum(K2)) )
	else
		d = Inf
	end

	return d
end