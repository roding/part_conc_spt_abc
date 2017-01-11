# Simple modulo function for imposing periodic boundary conditions.

function periodic(a::Float64, b::Float64)
	if a > b
		a = a - b
	elseif a < 0.0
		a = a + b
	end
	
	return a
end