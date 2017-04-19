function normpdf(x, mu, sigma)
	return 1.0 / ( sqrt(2.0 * pi) * sigma ) * exp( - 0.5 * (x - mu)^2 / sigma^2 )
end	