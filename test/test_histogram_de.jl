workspace()

include("../src/histogram_de.jl")

function test_histogram_de()
	N::Int64 = 100000000
	DE::Array{Float64, 1} = 10.0 * rand(N)
	
	de_bin_edges::FloatRange{Float64} = 0.0:1.0:10.0
	
	n_DE::Array{Int64, 1} = histogram_float64(DE, de_bin_edges)
	
	println(n)
	
	nothing
end

test_histogram_de()