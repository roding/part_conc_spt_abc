# Simple modulo function for imposing periodic boundary conditions.

function periodic(x::Float64, L::Float64)
	if x > L
		x = x - L
	elseif x < 0.0
		x = x + L
	end
	
	return x
end
x = 2.0
L = 1.5
x = periodic(x, L)
println(x)