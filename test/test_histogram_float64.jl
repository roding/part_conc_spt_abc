workspace()

include("../src/histogram_float64.jl")

function test_histogram_float64()
	N::Int64 = 100000000
	x::Array{Float64, 1} = 10.0 * rand(N)
	
	x_bin_edges::FloatRange{Float64} = 0.0:1.0:10.0
	
	n::Array{Int64, 1} = histogram_float64(x, x_bin_edges)
	
	println(n)
	
	nothing
end

test_histogram_float64()