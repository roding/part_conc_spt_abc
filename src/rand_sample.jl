function rand_weighted_integer(cumulative_weights::Array{Float64, 1})
	val::Float64 = rand()
	i::Int64 = 0
	while cumulative_weights[i] < val
		i += 1
	end
	return i
end