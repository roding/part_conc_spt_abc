function distance2(K1::Array{Int64, 1}, RSQ1::Array{Float64, 1}, K2::Array{Int64, 1}, RSQ2::Array{Float64, 1})
	n1::Int64 = length(K1)
	n2::Int64 = length(K2)
	
	d::Float64 = 0.0
	
	sort!(K1)
	sort!(K2)
	sort!(RSQ1)
	sort!(RSQ2)
	
	if n1 <= n2
		d = d + sum((K1 - K2[1:n1]).^2) + sum(K2[n1+1:end].^2)
	else
		d = d + sum((K1[1:n2] - K2).^2) + sum(K1[n2+1:end].^2)
	end
	if n1 <= n2
		d = d + sum((RSQ1 - RSQ2[1:n1]).^2) + sum(RSQ2[n1+1:end].^2)
	else
		d = d + sum((RSQ1[1:n2] - RSQ2).^2) + sum(RSQ1[n2+1:end].^2)
	end
	
	return d
end