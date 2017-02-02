function distance(K1::Array{Int64, 1}, DE1::Array{Float64, 1}, K2::Array{Int64, 1}, DE2::Array{Float64, 1})
	
	d::Float64 = 0.0

	if length(K1) > 0 && length(K2) > 0 # Then also D1 and D2 are non-empty.
		kmax = convert(Float64, max(maximum(K1), maximum(K2)))
		(~, n1) = hist(K1, 0.5:1:kmax+0.5)
		(~, n2) = hist(K2, 0.5:1:kmax+0.5)
		d = d + sum( (n1 - n2).^2 )
		(~, n1) = hist(DE1, 0.0:0.1:12.5)
		(~, n2) = hist(DE2, 0.0:0.1:12.5)
		d = d + sum( (n1 - n2).^2 )
	else
		d = Inf
	end

	return d
end