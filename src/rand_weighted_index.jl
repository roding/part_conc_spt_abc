function rand_weighted_index(cumulative_weights::Array{Float64, 1})
	val::Float64 = rand()
	idx::Int64 = 1
	while cumulative_weights[idx] < val
		idx += 1
	end
	return idx
end