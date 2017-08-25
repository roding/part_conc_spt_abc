function normpdf(x::Float64, mu::Float64, sigma::Float64)
	return 1.0 / ( sqrt(2.0 * pi) * sigma ) * exp( - 0.5 * (x - mu)^2 / sigma^2 )
end	