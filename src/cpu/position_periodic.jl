function position_periodic(x::Float64, L::Float64)
	if x > L
		x = x - L
	elseif x < 0.0
		x = x + L
	end

	return x
end
